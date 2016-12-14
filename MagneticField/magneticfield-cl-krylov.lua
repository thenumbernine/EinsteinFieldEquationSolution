#!/usr/bin/env luajit

local ffi = require 'ffi'
local math = require 'ext.math'
local range = require 'ext.range'

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
	local free = args.new
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

-- TODO put this in EinsteinFieldEquation/MagneticField or something?  this doesn't use LinearSolvers just yet
-- ye olde poisson problem
-- A x = y, y = rho, x = phi, A = del^2

local size = {2,2,2}
local env = require 'cl.obj.env'{size=size}

local rho = env:buffer{name='rho'}
local phi = env:buffer{name='phi'}

local program = env:program()

local init = program:kernel{
	argsOut = {phi, rho},
	body = [[ 
	rho[index] = i.x==size.x/2 && i.y==size.y/2 && i.z==size.z/2 ? -1 : 0;
	phi[index] = -rho[index];
]],
}

-- jacobi
local A = program:kernel{
	argsOut = {phi}, 
	argsIn = {rho}, 
	body = require 'template'([[
	const real4 dx = (real4)(1,1,1,1);
	
	real skewSum = 0;
	<? for i=0,gridDim-1 do ?>{	
		int4 iL = i;
		iL.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int indexL = indexForInt4(iL);
		real phiL = phi[indexL];
		
		int4 iR = i;
		iR.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int indexR = indexForInt4(iR);
		real phiR = phi[indexR];
		
		skewSum += (phiR + phiL) / (dx.s<?=i?> * dx.s<?=i?>);
	}<? end ?>

	const real diag = -2. * (0
<? for i=0,gridDim-1 do ?>
		+ 1. / (dx.s<?=i?> * dx.s<?=i?>)
<? end ?>
	);

	phi[index] = (rho[index] - skewSum) / diag;
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

local function gpustr(b)
	local c = b:toCPU()
	return range(env.volume):map(function(i) return c[i-1] end):map(tostring):concat'\n\t'
end

conjgradCL{
	A = function(y, x)
print(debug.traceback())
print('A')
print('x', gpustr(x))
		A(y.buf, x.buf)
print('y', gpustr(y))
	end,
	b = rho,
	x = phi,
	new = function() return env:buffer() end,
	free = function(buffer) buffer.buf:free() end,
	copy = function(dst, src) 
		env.cmds:enqueueCopyBuffer{
			src = src.buf,
			dst = dst.buf,
			-- TODO? store the size in the cl.buffer?  or the cl.obj.buffer?
			size = env.volume * ffi.sizeof(dst.type),
		}
	end,
	dot = function(a,b)
print(debug.traceback())
print('a', gpustr(a))
print('b', gpustr(b))
print('mul',dotBuf.buf,a.buf,b.buf)
		mul(
			assert(dotBuf.buf),
			assert(a.buf),
			assert(b.buf))
print('dotBuf', gpustr(dotBuf))
		local result = dot()
print('dot.result',dot.result)
print('dot.result[0]',dot.result[0])
print('dot() result',result)
		result = assert(tonumber(result))
print('tonumber(result)',result)
		assert(math.isfinite(result))
		return result
	end,
	mulAdd = function(y,a,b,s)
print(debug.traceback())
print('mulAdd',y.buf,a.buf,b.buf,s)
print('a',gpustr(a))
print('b',gpustr(b))
print('s',s)
		mulAdd(
			assert(y.buf),
			assert(a.buf),
			assert(b.buf), 
			ffi.new('real[1]', s))
print('y',gpustr(y))
	end,
	errorCallback = function(err, iter)
		print('err', err, 'iter', iter)
	end,
}

print('#x y z rho phi')
local rhoCPU = rho:toCPU()
local phiCPU = phi:toCPU()
local index = 0
for i=0,size[1]-1 do
	for j=0,size[2]-1 do
		for k=0,size[3]-1 do
			print(i,j,k,rhoCPU[index],phiCPU[index])
			index = index + 1
		end
	end
end
