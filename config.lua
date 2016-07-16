--[[
natural units ...
1 = c m/s = 299792458 m/s
	1 s = c m
	1 s = 299792458 m
1 = G m^3 / (kg s^2) = 6.67384e-11 m^3 / (kg s^2)
    kg = G m^3 / s^2 = G / c^2 m
	kg = 7.4256484500929e-28 m
1 = kB m^2 kg / (K s^2) = 1.3806488e-23 m^2 kg / (K s^2)
	K = kB kg m^2 / s^2 = kB / c^2 kg = kB G / c^4 m
	K = 1.1407124948367e-67 m
joules: J = kg m^2 / s^2
electronvolts: 1 eV = 1.6e-19 J
Gauss: 1 Gauss^2 = g / (cm s^2) = .1 kg / (m s^2)
	1 Gauss^2 = .1 G/c^4 1/m^2
	Gauss = sqrt(.1 G/c^2) 1/m

I'm going to use meters as my units ...

Radius of Earth = 6.37101e+6 m
Mass of Earth = 5.9736e+24 kg
--]]

c = 299792458	-- m/s 
G = 6.67384e-11	-- m^3 / (kg s^2)
kB = 1.3806488e-23	-- m^2 kg / (K s^2)

-- [[ earth
radius = 6.37101e+6	-- m
mass = 5.9736e+24 * G / c / c	-- m
--]]
--[[ sun
radius = 6.960e+8	-- m
mass = 1.9891e+30 * G / c / c	-- m
--]]

volume = 4/3*math.pi*radius^3	-- m^3
--earth volume: 1.0832120174985e+21 m^3
density = mass / volume	-- 1/m^2
--earth density: 4.0950296770075e-24 1/m^2
schwarzschildRadius = 2 * mass	--Schwarzschild radius: 8.87157 mm, which is accurate

--earth magnetic field at surface: .25-.26 gauss
magneticField = .45 * math.sqrt(.1 * G) / c	-- 1/m

local boundsRadius = 2
xmin = {-boundsRadius*radius,-boundsRadius*radius,-boundsRadius*radius}
xmax = {boundsRadius*radius,boundsRadius*radius,boundsRadius*radius}

size = 16
maxiter = 0

solver = 'conjgrad'
--solver = 'conjres'
--solver = 'gmres'

-- metric prims
--[[ flat
metricPrims = function(x,y,z)
	return 
		0,	-- ln(alpha)
		0,0,0,	-- beta^i
		1,0,0,1,0,1	-- gamma_ xx xy xz yy yz zz
end
--]]
-- [[ Schwarzschild stellar model.  MTW 23.2 box 6
metricPrims = function(x,y,z)
	local r = math.sqrt(x^2+y^2+z^2)
	local matterRadius = math.min(r, radius)
	local volumeOfMatterRadius = 4/3*math.pi*matterRadius^3
	local m = density * volumeOfMatterRadius	-- m^3		
	local ln_alpha = math.log(r > radius 
		and math.sqrt(1 - 2*mass/r)
		or (1.5 * math.sqrt(1 - 2*mass/radius) - .5 * math.sqrt(1 - 2*mass*r*r/(radius*radius*radius)))
	)
	local gamma_xx = 1 + x/r * x/r * 2*m/(r - 2*m)
	local gamma_yy = 1 + y/r * y/r * 2*m/(r - 2*m)
	local gamma_zz = 1 + z/r * z/r * 2*m/(r - 2*m)
	local gamma_xy = x/r * y/r * 2*m/(r - 2*m)
	local gamma_xz = x/r * z/r * 2*m/(r - 2*m)
	local gamma_yz = y/r * z/r * 2*m/(r - 2*m)
	return
		ln_alpha,	-- ln(alpha)
		0,0,0,	-- beta^i
		gamma_xx,gamma_xy,gamma_xz,gamma_yy,gamma_yz,gamma_zz	-- gamma_ij
end
--]]

-- stress energy prim - planet configuration
stressEnergyPrims = function(x,y,z)
	local r = math.sqrt(x^2+y^2+z^2)
	return 
		(r < radius and density or 0),	-- density
		0,	-- specific internal energy 
		0,	-- pressure
		0,0,0,	-- velocity
		0,0,0,	-- electric field
		0,0,0	-- magnetic field
end

analyticalGravity = function(x,y,z)
	local r = math.sqrt(x^2+y^2+z^2)
	local matterRadius = math.min(r, radius)
	local volumeOfMatterRadius = 4/3*math.pi*matterRadius^3
	local m = density * volumeOfMatterRadius	-- m^3
	--now that I'm using the correct alpha equation, my dm/dr term is causing the analytical gravity calculation to be off ...
	--local dm_dr = r > radius ? 0 : density * 4 * M_PI * matterRadius * matterRadius
	-- ... maybe it shouldn't be there to begin with?
	local dm_dr = 0
	local GammaUr_tt = (2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r)
		* c * c		--	+9 at earth surface, without matter derivatives
	--acceleration is -Gamma^r_tt along the radial direction (i.e. upwards from the surface), or Gamma^r_tt downward into the surface
	return GammaUr_tt
end

outputFilename = 'out.txt'
