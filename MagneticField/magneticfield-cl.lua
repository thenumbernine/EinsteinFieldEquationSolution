#!/usr/bin/env luajit
-- TODO put this in EinsteinFieldEquation/MagneticField or something?  this doesn't use LinearSolvers just yet
-- TODO use a linear solver?
-- ye olde poisson problem
-- A x = y, y = rho, x = phi, A = del^2

local ffi = require 'ffi'

local c = 299792458
local G = 6.67408e-11
local ke = 8.9875517873681764e+9	-- = 1 / (4 pi epsilon0)
local hBar = 1.05457180013e-34
local kB = 1.3806488e-23
local e = 6.2415093414e+18

local size = {64,64,64}
local env = require 'cl.obj.env'{size=size}

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
	
	local divs = 8 * size[1]
	for i=0,divs-1 do
		local frac = i / divs
		local th = 2 * math.pi * frac
		local x = math.floor(.5 * tonumber(env.domain.size.x) + .25 * size[1] * math.cos(th))
		local y = math.floor(.5 * tonumber(env.domain.size.y) + .25 * size[2] * math.sin(th))
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

-- [=[
-- jacobi
local iter = env:kernel{
	domain = env.domain,
	argsOut = {phi},
	argsIn = {rho},
	body=require 'template'([[
#if 1	//set boundary to zero?	
	if (i.x == 0 || i.x >= size.x-1 ||
		i.y == 0 || i.y >= size.y-1 ||
		i.z == 0 || i.z >= size.z-1)
	{
		phi[index] = 0.;
		return;
	}
#endif
	real skewSum = 0;
	<? for i=0,dim-1 do ?>{
		int4 iL = i;
		iL.s<?=i?> = max(0, iL.s<?=i?>-1);
		int indexL = indexForInt4(iL);

		int4 iR = i;
		iR.s<?=i?> = min(<?=tonumber(size:ptr()[i]-1)?>, iR.s<?=i?>+1);
		int indexR = indexForInt4(iR);

		skewSum += (phi[indexL] + phi[indexR]) * <?=1 / (dx:ptr()[i] * dx:ptr()[i])?>;
	}<? end ?>

<? 
local diag = 0
for i=0,dim-1 do
	diag = diag + 1. / (dx:ptr()[i] * dx:ptr()[i])
end
local invDiag = 1 / (-2 * diag)
?>	phi[index] = (rho[index] - skewSum) * <?=invDiag?>;
]], {
	clnumber = require 'cl.obj.number',
	dim = env.domain.dim,
	size = env.domain.size,
	dx = dx,
})}

for i=1,100 do
	iter()
end
--]=]

print'writing results...'
local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z rho phi\n')
local rhoCPU = rho:toCPU()
local phiCPU = phi:toCPU()
local index = 0
for k=0,size[3]-1 do
	for j=0,size[2]-1 do
		for i=0,size[1]-1 do
			file:write(i,'\t',j,'\t',k,'\t',rhoCPU[index],'\t',phiCPU[index],'\n')
			index = index + 1
		end
	end
end
file:close()
print'done!'
