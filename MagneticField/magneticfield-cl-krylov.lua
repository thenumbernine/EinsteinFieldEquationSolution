#!/usr/bin/env luajit

local template = require 'template'
local ffi = require 'ffi'

local c = 299792458
local G = 6.67408e-11
local ke = 8.9875517873681764e+9	-- = 1 / (4 pi epsilon0)
local hBar = 1.05457180013e-34
local kB = 1.3806488e-23
local e = 6.2415093414e+18

local n = 64
local env = require 'cl.obj.env'{size={n,n,n}}

-- init

local JU_CPU = ffi.new('real4[?]', env.domain.volume)
local AU_CPU = ffi.new('real4[?]', env.domain.volume)

for i=0,env.domain.volume do
	for j=0,3 do
		JU_CPU[i].s[j] = 0
	end
end

local vec3d = require 'ffi.vec.vec3d'
local xmin = vec3d(-1, -1, -1)
local xmax = vec3d(1, 1, 1)
local dx = (xmax - xmin) / env.domain.size
	
do	--lazy rasterization
	local divs = 8 * n
	for i=0,divs-1 do
		local frac = i / divs
		local th = 2 * math.pi * frac
		local x = math.floor(.5 * tonumber(env.domain.size.x) + .25 * n * math.cos(th))
		local y = math.floor(.5 * tonumber(env.domain.size.y) + .25 * n * math.sin(th))
		local z = math.floor(.5 * tonumber(env.domain.size.z)) 
	
		local q = 1							-- Coulombs (C)
		q = q * math.sqrt(ke * G) / (c * c)	-- ...times Columbs to meters (m)
		q = q * dx:volume()					-- ... per meter cubed (1/m^2)

		local index = x+env.domain.size.x*(y+env.domain.size.y+z)
		JU_CPU[index].s0 = q
		
		local I = q		-- amps (C / s) => 1/(m^2 s)
		I = I * c		-- (1/m^3)
		
		JU_CPU[index].s1 = -I * math.sin(th)
		JU_CPU[index].s2 = I * math.cos(th)
		JU_CPU[index].s3 = 0
	end
end

for i=0,env.domain.volume do
	for j=0,3 do
		AU_CPU[i].s[j] = -JU_CPU[i].s[j]
	end
end

local JU = env:buffer{name='JU', type='real4', data=JU_CPU}
local AU = env:buffer{name='AU', type='real4', data=AU_CPU}

local A = env:kernel{
	argsOut = {JU},
	argsIn = {AU},
	body = template([[	
#if 0	//set boundary to zero?	
	if (i.x == 0 || i.x >= size.x-1 ||
		i.y == 0 || i.y >= size.y-1 ||
		i.z == 0 || i.z >= size.z-1)
	{
		JU[index] = 0;	//boundary conditions
		return;
	} 
#endif

	const real4 dx = (real4)(<?=clnumber(dx.x)?>, <?=clnumber(dx.y)?>, <?=clnumber(dx.z)?>, 1);
	
	real4 sum = - AU[index] * 2. * (
		1. / (dx.x * dx.x)
		+ 1. / (dx.y * dx.y)
		+ 1. / (dx.z * dx.z)
	);
	<? for k=0,dim-1 do ?>{
		int4 iL = i;
		iL.s<?=k?> = max(0, iL.s<?=k?>-1);
		int indexL = indexForInt4(iL);

		int4 iR = i;
		iR.s<?=k?> = min(<?=tonumber(size:ptr()[k]-1)?>, iL.s<?=k?>+1);
		int indexR = indexForInt4(iR);

		sum += (AU[indexR] + AU[indexL]) / (dx.s<?=k?> * dx.s<?=k?>);
	}<? end ?>

	JU[index] = sum;
]], {
	dx = dx,
	clnumber = require 'cl.obj.number',
})}

local solver = 
--require 'LinearSolvers.cl.conjgrad'
require 'LinearSolvers.cl.conjres'	-- took 0.768707s to solve within 1e-7
{
	env = env,
	A = A,
	b = JU,
	x = AU,
	epsilon = 1e-7,

	-- hmm, for buffers of vector types,
	-- I need sizes not equal to the env volume (i.e. factors of it)
	-- and that means overriding 'n' here (for the buffer allocation)
	-- and that also means overriding 'size' for the buffer allocation
	type = 'real',
	size = env.volume * 4,
}

local beginTime = os.clock()
solver()
local endTime = os.clock()
print('took '..(endTime - beginTime)..' seconds')

print'writing results...'
JU:toCPU(JU_CPU)
AU:toCPU(AU_CPU)
local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z JU AU\n')
local index = 0
for i=0,tonumber(env.domain.size.x)-1 do
	for j=0,tonumber(env.domain.size.y)-1 do
		for k=0,tonumber(env.domain.size.z)-1 do
			file:write(i,'\t',j,'\t',k,
				'\t',JU_CPU[index],
				'\t',AU_CPU[index],'\n')
			index = index + 1
		end
	end
end
file:close()
print'done!'
