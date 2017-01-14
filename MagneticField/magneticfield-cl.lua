#!/usr/bin/env luajit
-- Poisson problem solved with Jacobi method
-- A x = y, y = JU, x = AU, A = del^2

local ffi = require 'ffi'

local c = 299792458
local G = 6.67408e-11
local ke = 8.9875517873681764e+9	-- = 1 / (4 pi epsilon0)
local hBar = 1.05457180013e-34
local kB = 1.3806488e-23
local e = 6.2415093414e+18

local n = 16	--64
local size = {16,16,16}
local env = require 'cl.obj.env'{size=size, verbose=true}

local JU_CPU = ffi.new('real4[?]', env.domain.volume)
local AU_CPU = ffi.new('real4[?]', env.domain.volume)

for i=0,env.domain.volume-1 do
	for j=0,3 do
		JU_CPU[i].s[j] = 0
	end
end

local vec3d = require 'ffi.vec.vec3d'
local xmin = vec3d(-1, -1, -1)
local xmax = vec3d(1, 1, 1)
local dx = (xmax - xmin) / vec3d(env.domain.size:unpack())
	
do	--lazy rasterization
	local q = 1							-- Coulombs (C)
	q = q * math.sqrt(ke * G) / (c * c)	-- ...times Columbs to meters (m)
	q = q * dx:volume()					-- ... per meter cubed (1/m^2)
	local I = q	-- amps (C / s) = 1/(m^2 s)
	I = I * c	-- 1/m^3
	local divs = 8 * n
	for i=0,divs-1 do
		local frac = i / divs
		local th = 2 * math.pi * frac
		local x = math.floor(.5 * size[1] + .25 * size[1] * math.cos(th))
		local y = math.floor(.5 * size[2] + .25 * size[2] * math.sin(th))
		local z = math.floor(.5 * size[3]) 
		local index = x + env.domain.size.x * (y + env.domain.size.y * z)
		assert(index >= 0 and index < env.domain.volume)
		JU_CPU[index].s0 = q
		JU_CPU[index].s1 = -I * math.sin(th)
		JU_CPU[index].s2 = I * math.cos(th)
		JU_CPU[index].s3 = 0
	end
end

for i=0,env.domain.volume-1 do
	for j=0,3 do
		AU_CPU[i].s[j] = -JU_CPU[i].s[j]
	end
end

local JU = env:buffer{name='JU', type='real4', data=JU_CPU}
local AU = env:buffer{name='AU', type='real4', data=AU_CPU}

-- [=[
-- jacobi
local iter = env:kernel{
	domain = env.domain,
	argsOut = {AU},
	argsIn = {JU},
	body=require 'template'([[
#if 1	//set boundary to zero?	
	if (i.x == 0 || i.x >= size.x-1 ||
		i.y == 0 || i.y >= size.y-1 ||
		i.z == 0 || i.z >= size.z-1)
	{
		AU[index] = 0.;
		return;
	}
#endif
	real4 skewSum = 0;
	<? for i=0,dim-1 do ?>{
		int4 iL = i;
		iL.s<?=i?> = max(0, iL.s<?=i?>-1);
		int indexL = indexForInt4(iL);

		int4 iR = i;
		iR.s<?=i?> = min(<?=tonumber(size:ptr()[i]-1)?>, iR.s<?=i?>+1);
		int indexR = indexForInt4(iR);

		skewSum += (AU[indexL] + AU[indexR]) * <?=1 / (dx:ptr()[i] * dx:ptr()[i])?>;
	}<? end ?>

<? 
local diag = 0
for i=0,dim-1 do
	diag = diag + 1. / (dx:ptr()[i] * dx:ptr()[i])
end
local invDiag = 1 / (-2 * diag)
?>	AU[index] = (JU[index] - skewSum) * <?=invDiag?>;
]], {
	clnumber = require 'cl.obj.number',
	dim = env.domain.dim,
	size = env.domain.size,
	dx = dx,
})}

local sumReduce = env:reduce{
	size = env.domain.volume * 4,
	type = 'real',
	op = function(x,y) return x..' + '..y end,
}

-- why does keeping the type real4 make a difference?!?!?!?!
-- and why does it fail when I only allocate room for 1 result in the reduce object?
-- [=[
local square = env:kernel{
	argsOut = {{name='result', type='real', buf=sumReduce.buffer}},
	argsIn = {{name='x', type='real', buf=AU.buf}},
	domain = require 'cl.obj.domain'{env=env, size=env.domain.volume * 4},
	body = [[	real v = x[index];	result[index] = v * v;	]],
}
--]=]
--[=[
local square = env:kernel{
	argsOut = {{name='result', type='real4', buf=sumReduce.buffer}},
	argsIn = {{name='x', type='real4', buf=AU.buf}},
	domain = env.domain.volume,
	body = [[	real4 v = x[index];	result[index] = v * v;	]],
}
--]=]
square:compile()
square.kernel:setArg(1, JU.buf)
square()
square.kernel:setArg(1, AU.buf)
local JULen = math.sqrt(sumReduce())	

-- TODO use a norm to find epsilon
print('# iter err')
for i=1,100 do
	iter()
	square()
	local AULen = math.sqrt(sumReduce())
	local err = AULen / JULen
	print(i, err)
	if err < 1e-5 then break end
end
--]=]

print'writing results...'
local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z J^t J^x J^y J^z A^t A^x A^y A^z\n')
JU:toCPU(JU_CPU)
AU:toCPU(AU_CPU)
local index = 0
for k=0,size[3]-1 do
	for j=0,size[2]-1 do
		for i=0,size[1]-1 do
			file:write(i,'\t',j,'\t',k,
				'\t',JU_CPU[index].s0,
				'\t',JU_CPU[index].s1,
				'\t',JU_CPU[index].s2,
				'\t',JU_CPU[index].s3,
				'\t',AU_CPU[index].s0,
				'\t',AU_CPU[index].s1,
				'\t',AU_CPU[index].s2,
				'\t',AU_CPU[index].s3,
				'\n')
			index = index + 1
		end
	end
end
file:close()
print'done!'
