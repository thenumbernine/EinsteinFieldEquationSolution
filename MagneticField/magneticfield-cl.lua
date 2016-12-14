#!/usr/bin/env luajit
-- TODO put this in EinsteinFieldEquation/MagneticField or something?  this doesn't use LinearSolvers just yet
-- TODO use a linear solver?
-- ye olde poisson problem
-- A x = y, y = rho, x = phi, A = del^2

local size = {64,64,64}
local env = require 'cl.obj.env'{size=size}

local rho = env:buffer{name='rho'}
local phi = env:buffer{name='phi'}

env:kernel{argsOut={phi, rho}, body=[[ 
	rho[index] = i.x==size.x/2 && i.y==size.y/2 && i.z==size.z/2 ? -1 : 0;
	phi[index] = -rho[index];
]]}()

-- jacobi
local iter = env:kernel{argsOut={phi}, argsIn={rho}, body=require 'template'([[
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

for i=1,100 do
	iter()
end

local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z rho phi\n')
local rhoCPU = rho:toCPU()
local phiCPU = phi:toCPU()
local index = 0
for i=0,size[1]-1 do
	for j=0,size[2]-1 do
		for k=0,size[3]-1 do
			file:write(i,'\t',j,'\t',k,'\t',rhoCPU[index],'\t',phiCPU[index],'\n')
			index = index + 1
		end
	end
end
file:close()
