#!/usr/bin/env luajit
-- lapse.gnuplot implemented in lua, using the gnuplot library
-- this way I can do more stuff with the data
require 'ext'
local gnuplot = require 'gnuplot'

local n = 1000
local I = function(i) return i end
local C = function(c) return function() return c end end
local U = function(f) return range(n):map(f) end
local data = function(t) return table.map(t, U) end

local function inc()
	M = m(R)	--  matter within planet radius, in m

	--  lapse 
	--  lapse without considering hydrostatic pressure (i.e. if m(r) = m(radius), i.e. if radius=0, i.e. a black hole)
	-- schwarzschild_alpha(r) = sqrt(1 - 2*m(r)/r)
	--  lapse with considerations (MTW box 23.2 eqn 6)
	schwarzschild_alpha = function(r) 
		return r > R 
			and math.sqrt(1 - 2*M/r) 
			or (1.5*math.sqrt(1 - 2*M/R) - .5*math.sqrt(1 - 2*M*r*r/(R*R*R))) 
	end

	--  derivative of lapse wrt radius
	schwarzschild_dr_alpha = function(r)
		return r > R 
			and (M / (r * r * math.sqrt(1 - 2*M/r))) 
			or M*r/(R*R*R*math.sqrt(1 - 2*M*r*r/(R*R*R)))
	end

	--  metric, based on g_tt = -alpha^2 + beta^2 and assuming beta=0
	schwarzschild_g_tt = function(r)
		return -schwarzschild_alpha(r)^2
	end

	--  schwarzschild gravity = Gamma_rtt = -alpha d/dr alpha ... times c^2 for m/s^2
	schwarzschild_gravity = function(r) 
		return -schwarzschild_alpha(r) * schwarzschild_dr_alpha(r) * c^2
	end

	--  newton gravity ... times c^2 for m/s^2
	newton_gravity = function(r) 
		return -m(r) / r^2 * c^2
	end

	--  now for kerr.
	--  I'll just compute the gravity
	--  there is no solution yet solved for the equations of structure inside a star using the Kerr metric
	angularVelocity = 2 * math.pi / (60 * 60 * 24) / c	--  angular velocity, in m^-1
	inertia = 2 / 5 * M * R^2						--  moment of inertia about a sphere, in m^3
	angularMomentum = inertia * angularVelocity			--  angular momentum in m^2
	earth_a = angularMomentum / M						--  in m

	Delta = function(r,a)
		return r^2 - 2*m(r) * r + a^2
	end
	Sigma = function(r,a,theta)
		return r^2 + a^2 * math.cos(theta)^2
	end
	kerr_gravity = function(r,a,theta) 
		return -2*m(r) * Delta(r,a) * (r^2 - a^2 * math.cos(theta)^2) / (2 * Sigma(r,a,theta)^3) * c^2
	end
end

sphereVolume = function(r) return 4/3*math.pi*r*r*r end
sphereSurface = function(r) return 4*math.pi*r*r end

c = 299792458						--  ... m/s = 1
G = 6.67384e-11						--  ... m^3 / (kg s^2) = 1

km_in_m = 1000	-- 1 km = 1000 m
cm_in_m = 1/100
kg_in_g = 1000
g_in_kg = 1/kg_in_g
kg_in_m = G / c / c
Pa_in_m = kg_in_m / c / c 

theta = math.pi/2						--  polar angle for Kerr metric

--  set radius and density to temporary values:
R = 1
density = .1
m = function(x) return density * sphereVolume(math.min(R,x)) end			--  matter within radius, in m
dm_dr = function(x) return R < x and 0 or (density * sphereSurface(x)) end	--  derivative of schwarzschild radius function wrt radius
inc()

--  show how schwarzschild alpha, d/dr alpha, g_tt, and Gamma^r_tt are all related:
local x = (I-.5)/n * R*10
gnuplot{
	output = 'images/schwarzschild_eos.png',
	style = 'data lines',
	data = data{x,
		schwarzschild_alpha:o(x),
		schwarzschild_dr_alpha:o(x),
		schwarzschild_g_tt:o(x),
		-schwarzschild_alpha:o(x) * schwarzschild_dr_alpha:o(x),
	},
	{using='1:2', title='schwarzschild alpha'},
	{using='1:3', title='schwarzschild d_r alpha'},
	{using='1:4', title='g_t_t'},
	{using='1:5', title='schwarzschild gravity = -Gamma_r_t_t = -alpha d_r alpha'},
}

--  set radius and density to earth parameters:
R = 6.37101e+6						--  matter radius, in m
mass = 5.9736e+24 * G / c / c		--  total mass within radius, in m
density = mass / sphereVolume(R)	--  average density, in m^-2
m = function(x) return density * sphereVolume(math.min(R,x)) end			--  matter within radius, in m
dm_dr = function(x) return R < x and 0 or (density * sphereSurface(x)) end	--  derivative of schwarzschild radius function wrt radius
inc()

kerr_gravity_earth_rotation = function(x) return kerr_gravity(x, earth_a, theta) end
kerr_gravity_no_rotation = function(x) return kerr_gravity(x, 0, theta) end

--  show gravity only, in m/s^2
local x = (I-.5)/n * R*10
gnuplot{
	output = 'images/gravity models.png',
	log = 'y',
	style = 'data lines',
	ylabel = 'm/s^2',
	data = data{x,
		-newton_gravity:o(x),
		-schwarzschild_gravity:o(x),
		-kerr_gravity_earth_rotation:o(x),
	},
	{using='1:2', title='newton gravity'},
	{using='1:3', title='schwarzschild gravity'},
	{using='1:4', title='kerr gravity'},
}

--  show differences in gravity 
local x = (I-.5)/n * R*2
gnuplot{
	output = 'images/gravity differences.png',
	style = 'data lines',
	ylabel = 'm/s^2',
	data = data{x,
		math.abs:o(schwarzschild_gravity:o(x)) - math.abs:o(newton_gravity:o(x)),
		math.abs:o(schwarzschild_gravity:o(x)) - math.abs:o(kerr_gravity_earth_rotation:o(x)),
		math.abs:o(newton_gravity:o(x)) - math.abs:o(kerr_gravity_earth_rotation:o(x)),
	},
	{using='1:2', title='|schwarzschild|-|newton|'},
	{using='1:3', title='|schwarzschild|-|kerr|'},
	{using='1:4', title='|newton|-|kerr|'},
}

--  rotation-less Kerr isn't the same as Schwarzschild ... It weaker than rotation Kerr.  It is within 2e-10 of rotation Kerr ... 2e-11 at the Earth's surface
--  note that Schwarzschild (which is rotation-less) is stronger than rotation-Kerr and rotation-less Kerr by about 1.4e-8

local x = (I-.5)/n * R*2
gnuplot{
	output = 'images/kerr with vs. without rotation.png',
	style = 'data lines',
	ylabel = 'm/s^2',
	log = 'y',
	data = data{x, math.abs:o(kerr_gravity_earth_rotation:o(x)) - math.abs:o(kerr_gravity_no_rotation:o(x))},
	{using='1:2', title='|kerr w/rotation|-|kerr w/o rotation|'},
}

local x = (I-.5)/n * R + R
gnuplot{
	output = 'images/kerr fast rotation vs. kerr rotation.png',
	style = 'data lines',
	ylabel = 'm/s^2',
	data = data{x, 
		function(i) return math.abs(kerr_gravity(x(i),520000*earth_a,theta)) - math.abs(kerr_gravity(x(i),earth_a,theta)) end,
	},
	{using='1:2', title='|kerr fast rotation|-|kerr w/rotation| ... how to double gravity only using rotation'},
}

-- constant density vs PREM
local premdata = file['prem.txt']:trim():split'\n':filter(function(l)
	return l:sub(1,1) ~= '#'
end):map(function(l)
	return l:split'%s+':map(tonumber:o(I))
end)
-- reverse
premdata = premdata:map(function(_,i)
	return premdata[#premdata-i+1]
end)
local prem = function(i) return premdata[i] end
-- convert to meters
local prem_r = premdata:map(function(row) return row[2] * km_in_m end)	-- km -> m
local prem_rho = premdata:map(function(row) return row[5] * g_in_kg * kg_in_m / cm_in_m^3 end)		-- g/cm^3 -> m^-2
local prem_P = premdata:map(function(row) return row[8] * 1e+9 * Pa_in_m end)	-- GPa -> m^-2

-- now let's try to integrate prem_rho through the sphere ...
local function integral(f, x1, x2)
	local n = 2000 * math.ceil(x2 - x1) / 6e+6
	local sum = 0
	local h = (x2 - x1) / n
	for i=0,n do
		local x = i/n * (x2 - x1) + x1
		sum = sum + f(x) * (
			(i == 0 or i == n) and 1
			or ((i%2==1) and 4 or 2)
		)
	end
	return sum * h / 3
end
local function rho(r) 
	for i=1,#premdata-1 do
		if prem_r[i+1] > r then
			local f = (r - prem_r[i]) / (prem_r[i+1] - prem_r[i])
			return (1-f) * prem_rho[i] + f * prem_rho[i+1]
		end
	end
	return 0
end

m = function(R) return integral(function(r) return 4 * math.pi * rho(r) * r^2 end, 0, R) end
local epsilon = 1e-6
dm_dr = function(R) return (m(R+epsilon) - m(R-epsilon)) / (2 * epsilon) end
local x = (I-.5)/n * R*2
gnuplot{
	output = 'images/prem mass.png',
	style = 'data lines',
	data = data{x, m:o(x)},
	{using='1:2', title='m(r)'},
}
gnuplot{
	output = 'images/prem d_dr mass.png',
	style = 'data lines',
	data = data{x, dm_dr:o(x)},
	{using='1:2', title='dm/dr'},
}

inc()
local x = (I-.5)/n * R*10
gnuplot{
	output = 'images/prem gravity models.png',
	log = 'y',
	style = 'data lines',
	ylabel = 'm/s^2',
	data = data{x,
		-newton_gravity:o(x),
		-schwarzschild_gravity:o(x),
		-kerr_gravity_earth_rotation:o(x),
	},
	{using='1:2', title='newton gravity'},
	{using='1:3', title='schwarzschild gravity'},
	{using='1:4', title='kerr gravity'},
}

--  show differences in gravity 
local x = (I-.5)/n * R*2
gnuplot{
	output = 'images/prem gravity differences.png',
	style = 'data lines',
	ylabel = 'm/s^2',
	data = data{x,
		math.abs:o(schwarzschild_gravity:o(x)) - math.abs:o(newton_gravity:o(x)),
		math.abs:o(schwarzschild_gravity:o(x)) - math.abs:o(kerr_gravity_earth_rotation:o(x)),
		math.abs:o(newton_gravity:o(x)) - math.abs:o(kerr_gravity_earth_rotation:o(x)),
	},
	{using='1:2', title='|schwarzschild|-|newton|'},
	{using='1:3', title='|schwarzschild|-|kerr|'},
	{using='1:4', title='|newton|-|kerr|'},
}
