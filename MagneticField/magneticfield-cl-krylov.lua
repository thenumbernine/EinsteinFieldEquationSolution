#!/usr/bin/env luajit

local template = require 'template'

local env = require 'cl.obj.env'{size={64,64,64}}

local rho = env:buffer{name='rho'}
local phi = env:buffer{name='phi'}

-- init
env:kernel{
	argsOut = {phi, rho},
	body = [[ 
	const real M = 1e+6;
	rho[index] = abs(i.x - size.x/2) < 1 
				&& abs(i.y - size.y/2) < 1
				&& abs(i.z - size.z/2) < 1
				? (
					i.x >= size.x/2 ? -M : M
				) : 0;
	phi[index] = -rho[index];
]],
}()

local A = env:kernel{
	argsOut = {{name='y', buf=true}}, 
	argsIn = {{name='x', buf=true}}, 
	body = template([[	
	if (i.x == 0 || i.x >= size.x-1 ||
		i.y == 0 || i.y >= size.y-1 ||
		i.z == 0 || i.z >= size.z-1)
	{
		y[index] = 0;	//boundary conditions
	} else {
		const real4 dx = (real4)(<?=clnumber(1/tonumber(size.x))?>, <?=clnumber(1/tonumber(size.y))?>, <?=clnumber(1/tonumber(size.z))?>, 1);
		y[index] = (x[index + stepsize.x] + x[index - stepsize.x]) / (dx.x * dx.x)
					+ (x[index + stepsize.y] + x[index - stepsize.y]) / (dx.y * dx.y)
					+ (x[index + stepsize.z] + x[index - stepsize.z]) / (dx.z * dx.z)
					- x[index] * 2. * (
						1. / (dx.x * dx.x)
						+ 1. / (dx.y * dx.y)
						+ 1. / (dx.z * dx.z)
					);
	}
]], {
	dim = env.dim,
	size = env.size,
	clnumber = require 'cl.obj.number',
})}

local solver = 
--require 'LinearSolvers.cl.conjgrad'
require 'LinearSolvers.cl.conjres'
{
	env = env,
	A = A,
	b = rho,
	x = phi,
	errorCallback = function(err, iter)
		print('err', err, 'iter', iter)
	end,
}

solver()

print'writing results...'
local rhoCPU = rho:toCPU()
local phiCPU = phi:toCPU()
local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z rho phi\n')
local index = 0
for i=0,tonumber(env.size.x)-1 do
	for j=0,tonumber(env.size.y)-1 do
		for k=0,tonumber(env.size.z)-1 do
			file:write(i,'\t',j,'\t',k,
				'\t',rhoCPU[index],
				'\t',phiCPU[index],'\n')
			index = index + 1
		end
	end
end
file:close()
print'done!'
