#!/usr/bin/env luajit

local ffi = require 'ffi'

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

end

-- TODO put this in EinsteinFieldEquation/MagneticField or something?  this doesn't use LinearSolvers just yet
-- ye olde poisson problem
-- A x = y, y = rho, x = phi, A = del^2

local env = require 'cl.obj.env'{size={8,8,8}}

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
	argsOut = {{name='A_phi', buf=true}}, 
	argsIn = {{name='phi', buf=true}}, 
	body = require 'template'([[
	const real4 dx = (real4)(1,1,1,1);
	
	real sum = 0;
	<? for i=0,gridDim-1 do ?>{	
		int4 iL = i;
		iL.s<?=i?> = max(i.s<?=i?> - 1, 0);
		int indexL = indexForInt4(iL);
		
		int4 iR = i;
		iR.s<?=i?> = min(i.s<?=i?> + 1, size.s<?=i?> - 1);
		int indexR = indexForInt4(iR);
		
		sum += (phi[indexL] + phi[indexR]) / (dx.s<?=i?> * dx.s<?=i?>);
	}<? end ?>

	sum -= 2. * (0
<? for i=0,gridDim-1 do ?>
		+ 1. / (dx.s<?=i?> * dx.s<?=i?>)
<? end ?>
	) * phi[index];

	A_phi[index] = sum;
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

conjgradCL{
	A = function(y, x)
		A(y.buf, x.buf)
	end,
	b = rho,
	x = phi,
	maxiter = env.volume,
	new = function() return env:buffer() end,
	free = function(buffer) buffer.buf:release() end,
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

local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z rho phi\n')
local rhoCPU = rho:toCPU()
local phiCPU = phi:toCPU()
local index = 0
for i=0,tonumber(env.size.x)-1 do
	for j=0,tonumber(env.size.y)-1 do
		for k=0,tonumber(env.size.z)-1 do
			file:write(i,'\t',j,'\t',k,'\t',rhoCPU[index],'\t',phiCPU[index],'\n')
			index = index + 1
		end
	end
end
file:close()
