#!/usr/bin/env luajit

-- same as magneticfield-cl-krylov except only solving the electric potential

local template = require 'template'
local ffi = require 'ffi'

local c = 299792458
local G = 6.67408e-11
local ke = 8.9875517873681764e+9	-- = 1 / (4 pi epsilon0)
local hBar = 1.05457180013e-34
local kB = 1.3806488e-23
local e = 6.2415093414e+18

print'creating env...'
local n = 64
local env = require 'cl.obj.env'{size={n,n,n}}

-- init

local rhoCPU = ffi.new('real[?]', env.domain.volume)
local phiCPU = ffi.new('real[?]', env.domain.volume)

for i=0,env.domain.volume-1 do
	rhoCPU[i] = 0
end

local vec3d = require 'ffi.vec.vec3d'
local xmin = vec3d(-1, -1, -1)
local xmax = vec3d(1, 1, 1)
local dx = (xmax - xmin) / vec3d(env.domain.size:unpack())
	
do	--lazy rasterization
	local q = 1							-- Coulombs (C)
	q = q * math.sqrt(ke * G) / (c * c)	-- ...times Columbs to meters (m)
	q = q * dx:volume()					-- ... per meter cubed (1/m^2)
	
	local divs = 8 * n
	for i=0,divs-1 do
		local frac = i / divs
		local th = 2 * math.pi * frac
		local x = math.floor(.5 * tonumber(env.domain.size.x) + .25 * n * math.cos(th))
		local y = math.floor(.5 * tonumber(env.domain.size.y) + .25 * n * math.sin(th))
		local z = math.floor(.5 * tonumber(env.domain.size.z)) 
	
		local index = x + env.domain.size.x * (y + env.domain.size.y * z)
		rhoCPU[index] = q
	end
end

for i=0,env.domain.volume-1 do
	phiCPU[i] = -rhoCPU[i]
end

local rho = env:buffer{name='rho', data=rhoCPU}
local phi = env:buffer{name='phi', data=phiCPU}

local A = env:kernel{
	argsOut = {rho},
	argsIn = {phi},
	body = template([[	
#if 1	//set boundary to zero?	
	if (i.x == 0 || i.x >= size.x-1 ||
		i.y == 0 || i.y >= size.y-1 ||
		i.z == 0 || i.z >= size.z-1)
	{
		rho[index] = 0;	//boundary conditions
		return;
	} 
#endif
	real sum = -phi[index] * <?= 2. * (
		1. / (dx.x * dx.x)
		+ 1. / (dx.y * dx.y)
		+ 1. / (dx.z * dx.z)
	)?>;
	<? for k=0,dim-1 do ?>{
		int4 iL = i;
		iL.s<?=k?> = max(0, iL.s<?=k?>-1);
		int indexL = indexForInt4(iL);

		int4 iR = i;
		iR.s<?=k?> = min(<?=tonumber(size:ptr()[k]-1)?>, iR.s<?=k?>+1);
		int indexR = indexForInt4(iR);

		sum += (phi[indexR] + phi[indexL]) * <?=1 / (dx:ptr()[k] * dx:ptr()[k]) ?>;
	}<? end ?>
	
	rho[index] = sum;
]], {
	dx = dx,
	dim = env.domain.dim,
	size = env.domain.size,
	clnumber = require 'cl.obj.number',
})}

local solver = 
--require 'solver.cl.conjgrad'
require 'solver.cl.conjres'	-- took 0.768707s to solve within 1e-7
--require 'solver.cl.gmres'
{
	env = env,
	A = A,
	b = rho,
	x = phi,
	epsilon = 1e-7,
	errorCallback = function(err,iter)
		io.stderr:write(tostring(err)..'\t'..tostring(iter)..'\n')
		assert(err == err)
	end,
}

local beginTime = os.clock()
solver()
local endTime = os.clock()
print('took '..(endTime - beginTime)..' seconds')

print'writing results...'
rho:toCPU(rhoCPU)
phi:toCPU(phiCPU)
local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z rho phi\n')
local index = 0
for k=0,tonumber(env.domain.size.z)-1 do
	for j=0,tonumber(env.domain.size.y)-1 do
		for i=0,tonumber(env.domain.size.x)-1 do
			file:write(i,'\t',j,'\t',k,
				'\t',rhoCPU[index],
				'\t',phiCPU[index],'\n')
			index = index + 1
		end
	end
end
file:close()
print'done!'
