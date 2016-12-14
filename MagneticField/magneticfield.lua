#!/usr/bin/env luajit
local ni = ... and tonumber(...) or 64
local conjgrad = require 'LinearSolvers.ConjugateGradient'
local conjres = require 'LinearSolvers.ConjugateResidual'
local gmres = require 'LinearSolvers.GeneralizedMinimalResidual'
local matrix = require 'matrix'

local n = {ni,ni,ni}
local h2 = 1/ni^2

local M = 1e+6
local rho = matrix.lambda(n, function(i,j,k)
	return (math.abs(i - (n[1]+1)/2) < 1
		and math.abs(j - (n[2]+1)/2) < 1
		and math.abs(k - (n[3]+1)/2) < 1)
		and (
			i >= (n[1]+1)/2 and -M or M
		) or 0
end)

local phi = 
conjgrad
--conjres
--gmres
{
	x = -rho,
	b = rho,
	A = function(phi)
		return matrix.lambda(n, function(i,j,k)
			if i==1 or i==n[1]
			or j==1 or j==n[2]
			or k==1 or k==n[3]
			then
				return 0
			else
				return (phi[i+1][j][k]
					+ phi[i-1][j][k]
					+ phi[i][j+1][k]
					+ phi[i][j-1][k]
					+ phi[i][j][k+1]
					+ phi[i][j][k-1]
					- 6 * phi[i][j][k]) / h2
			end
		end)
	end,
	clone = matrix,
	dot = matrix.dot,
	errorCallback = function(err,iter)
		io.stderr:write(tostring(err)..'\t'..tostring(iter)..'\n')
	end,
	epsilon = 1e-5,
	maxiter = ni^3,
	restart = 100,	--gmres-only
}

local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z rho phi\n')
for i=1,n[1] do
	for j=1,n[2] do
		for k=1,n[3] do
			file:write(i-1,'\t',j-1,'\t',k-1,'\t',rho[i][j][k],'\t',phi[i][j][k],'\n')
		end
	end
end
file:close()
