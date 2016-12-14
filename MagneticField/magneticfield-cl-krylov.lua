#!/usr/bin/env luajit

local ffi = require 'ffi'
local template = require 'template'

--[[
just like in LinearSolvers.ConjugateGradient, except with favor towards reusing object (for the sake of CL buffers)
	A = function(y,x) for x ro, y rw
		applies the linear function of A
		reads from x vector, writes to y vector
	b = object to hold 'b' vector
	x = object to hold 'x' vector
	MInv = (optional) function(y,x) for x ro and y rw vectors. preconditioner
	errorCallback (optional) returns true to stop
	epsilon (optional)
	maxiter (optional)
	
	new = function() returns and create a new vector
	free = function(v) frees vector
	copy = function(dst, src) copies contents of src into dst

	dot = function(a,b) returns a number of the inner product of a and b
	mulAdd = function(y,a,b,c) y = a + b * c, for y,a,b vectors and c scalar
--]]
local function conjgradCL(args)
	local A = assert(args.A)
	local b = assert(args.b)
	local x = assert(args.x)
	
	local MInv = args.MInv
	local errorCallback = args.errorCallback
	local epsilon = args.epsilon or 1e-7
	local maxiter = args.maxiter or 1000

	local copy = assert(args.copy)
	local new = assert(args.new)
	local free = args.free
	local dot = assert(args.dot)
	local mulAdd = assert(args.mulAdd)

	local r = new()
	local p = new()
	local Ap = new()
	local MInvR = MInv and new() or r

	local bNorm = dot(b,b)
	if bNorm == 0 then bNorm = 1 end

	A(r, x)
	mulAdd(r, b, r, -1)
	
	if MInv then MInv(MInvR, r) end
	local rDotMInvR = dot(r, MInvR)

	repeat
		local err = dot(r, r) / bNorm
		if errorCallback and errorCallback(err, 0) then break end
		if err < epsilon then break end
		
		copy(p, MInvR)
		for iter=1,maxiter do
			A(Ap, p)
			local alpha = rDotMInvR / dot(p, Ap)
			mulAdd(x, x, p, alpha)
			mulAdd(r, r, Ap, -alpha)
			
			local err = dot(r, r) / bNorm
			if errorCallback and errorCallback(err, iter) then break end
			if err < epsilon then break end
			
			if MInv then MInv(MInvR, r) end
			local nRDotMInvR = dot(r, MInvR)
			local beta = nRDotMInvR / rDotMInvR
			mulAdd(p, MInvR, p, beta)
			
			rDotMInvR = nRDotMInvR
		end
	until true	-- run once, use break to jump out. my stupid CS education has scarred me from ever using goto's again.

	if free then 
		free(r) 
		free(p)
		free(Ap)
		if MInv then free(MInvR) end
	end

	return x
end

local function conjresCL(args)
	local A = assert(args.A)
	local b = assert(args.b)
	local x = assert(args.x)

	local MInv = args.MInv
	local errorCallback = args.errorCallback
	local epsilon = args.epsilon or 1e-7
	local maxiter = args.maxiter or 1000

	local copy = assert(args.copy)
	local new = assert(args.new)
	local free = args.new
	local dot = assert(args.dot)
	local mulAdd = assert(args.mulAdd)

	local r = new()
	local p = new()
	local Ap = new()
	local Ar = new()
	local MInvAp = MInv and new() or Ap

	local bNorm = dot(b,b)
	if bNorm == 0 then bNorm = 1 end

	A(r, x)
	mulAdd(r, b, r, -1)
	
	if MInv then MInv(r, r) end
	local rDotMInvR = dot(r, r)

	repeat
		local err = dot(r, r) / bNorm
		if errorCallback and errorCallback(err, 0) then break end
		if err < epsilon then break end

		A(Ar, r)
		local rAr = dot(r, Ar)
		copy(p, r)
		A(Ap, p)
		for iter=1,maxiter do
			if MInv then MInv(MInvAp, Ap) end
			local alpha = rAr / dot(Ap, MInvAp)
			mulAdd(x, x, p, alpha)
			mulAdd(r, r, MInvAp, -alpha)
		
			local err = dot(r, r) / bNorm
			if errorCallback and errorCallback(err, iter) then break end
			if err < epsilon then break end
		
			A(Ar, r)
			local nrAr = dot(r, Ar)
			local beta = nrAr

			rAr = nrAr
			mulAdd(p, r, p, beta)
			mulAdd(Ap, Ar, Ap, beta)
		end
	
	until true -- just run once / use for break jumps

	if free then
		free(r)
		free(p)
		free(Ap)
		free(Ar)
		if MInv then free(MInvAp) end
	end

	return x
end

-- TODO put this in EinsteinFieldEquation/MagneticField or something?  this doesn't use LinearSolvers just yet
-- ye olde poisson problem
-- A x = y, y = rho, x = phi, A = del^2

local env = require 'cl.obj.env'{size={64,64,64}}

local rho = env:buffer{name='rho'}
local phi = env:buffer{name='phi'}

local program = env:program{
	code = template([[
constant real4 dx = (real4)(<?=clnumber(1/tonumber(size.x))?>, <?=clnumber(1/tonumber(size.y))?>, <?=clnumber(1/tonumber(size.z))?>, 1);
]], {
	size = env.size,
	clnumber = require 'cl.obj.number',
}),
}

local init = program:kernel{
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
}

-- jacobi
local A = program:kernel{
	argsOut = {{name='A_phi', buf=true}}, 
	argsIn = {{name='phi', buf=true}}, 
	body = template([[	
	if (i.x == 0 || i.x >= size.x-1 ||
		i.y == 0 || i.y >= size.y-1 ||
		i.z == 0 || i.z >= size.z-1)
	{
		A_phi[index] = 0;	//boundary conditions
	} else {
		A_phi[index] = (phi[index + stepsize.x] + phi[index - stepsize.x]) / (dx.x * dx.x)
					+ (phi[index + stepsize.y] + phi[index - stepsize.y]) / (dx.y * dx.y)
					+ (phi[index + stepsize.z] + phi[index - stepsize.z]) / (dx.z * dx.z)
					- phi[index] * 2. * (
						1. / (dx.x * dx.x)
						+ 1. / (dx.y * dx.y)
						+ 1. / (dx.z * dx.z)
					);
	}
]], {
	gridDim = env.gridDim,
	gridSize = env.size,
})}

local dotBuf = env:buffer{name='dotBuf'}
local mul = program:kernel{
	argsOut = {{name='y', buf=true}},
	argsIn = {{name='a', buf=true}, {name='b', buf=true}},
	body = [[	y[index] = a[index] * b[index];]],
}

local mulAdd = program:kernel{
	argsOut = {{name='y', buf=true}},
	argsIn = {{name='a', buf=true}, {name='b', buf=true}, {name='s', type='real'}},
	body = [[	y[index] = a[index] + b[index] * s;]],
}

program:compile()

local dot = env:reduce{
	buffer = dotBuf.buf,
	op = function(x,y) return x..' + '..y end,
}

init()

conjgradCL
--conjresCL	-- not working yet
{
	A = function(y, x) A(y.buf, x.buf) end,
	b = rho,
	x = phi,
	maxiter = env.volume,
	new = function() return env:buffer() end,
	
	-- hmm, this release, coupled with the __gc's release, makes things crash ...
	-- at least I know the __gc is cleaning up correctly
	--free = function(buffer) buffer.buf:release() end,
	
	copy = function(dst, src) 
		env.cmds:enqueueCopyBuffer{
			src = src.buf,
			dst = dst.buf,
			-- TODO? store the size in the cl.buffer?  or the cl.obj.buffer?
			size = env.volume * ffi.sizeof(dst.type),
		}
	end,
	dot = function(a,b)
		mul(dotBuf.buf, a.buf, b.buf)
		return dot()
	end,
	mulAdd = function(y,a,b,s)
		mulAdd(y.buf, a.buf, b.buf, ffi.new('real[1]', s))
	end,
	errorCallback = function(err, iter)
		print('err', err, 'iter', iter)
	end,
}

print'writing results...'
local rhoCPU = rho:toCPU()
local phiCPU = phi:toCPU()
local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z rho phi\n')
local index = 0
for i=0,tonumber(env.size.x)-1 do
	for j=0,tonumber(env.size.y)-1 do
		for k=0,tonumber(env.size.z)-1 do
			local line = tostring(i)
				..'\t'..tostring(j)
				..'\t'..tostring(k)
				..'\t'..tostring(rhoCPU[index])
				..'\t'..tostring(phiCPU[index])
				..'\n'
			file:write(line)
			index = index + 1
		end
	end
end
file:close()
print'done!'
