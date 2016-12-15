#!/usr/bin/env luajit
local ni = ... and tonumber(...) or 64
local conjgrad = require 'LinearSolvers.ConjugateGradient'
local conjres = require 'LinearSolvers.ConjugateResidual'
local gmres = require 'LinearSolvers.GeneralizedMinimalResidual'
local matrix = require 'matrix'
local vec3 = require 'vec.vec3'
local table = require 'ext.table'

local n = {ni,ni,ni}
local h2 = 1/ni^2

local c = 299792458
local G = 6.67408e-11
local ke = 8.9875517873681764e+9	-- = 1 / (4 pi epsilon0)
local hBar = 1.05457180013e-34
local kB = 1.3806488e-23
local e = 6.2415093414e+18

local xmin = vec3(-1, -1, -1)
local xmax = vec3(1, 1, 1)
local dx = vec3()
for i=1,3 do
	dx[i] = (xmax[i] - xmin[i]) / n[i]
end
	
local rho = matrix.zeros(table.unpack(n))
do	--lazy rasterization
	local q = 1							-- Coulombs (C)
	q = q * math.sqrt(ke * G) / (c * c)	-- ...times Columbs to meters (m)
	q = q * dx[1] * dx[2] * dx[3]					-- ... per meter cubed (1/m^2)
	local divs = 8 * ni
	for i=0,divs-1 do
		local frac = i / divs
		local th = 2 * math.pi * frac
		local x = math.floor(.5 * n[1] + .25 * ni * math.cos(th))
		local y = math.floor(.5 * n[2] + .25 * ni * math.sin(th))
		local z = math.floor(.5 * n[3]) 

		rho[x+1][y+1][z+1] = q
	end
end

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
		print(err,iter)
	end,
	epsilon = 1e-5,
	maxiter = ni^3,
	restart = 100,	--gmres-only
}

local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z rho phi\n')
for k=1,n[3] do
	for j=1,n[2] do
		for i=1,n[1] do
			file:write(i-1,'\t',j-1,'\t',k-1,'\t',rho[i][j][k],'\t',phi[i][j][k],'\n')
		end
	end
end
file:close()
