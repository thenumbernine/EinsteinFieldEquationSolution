#!/usr/bin/env luajit
local ni = ... and tonumber(...) or 16 --64
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
	
local JU = matrix.zeros(n[1], n[2], n[3], 4)
do	--lazy rasterization
	local q = 1							-- Coulombs (C)
	q = q * math.sqrt(ke * G) / (c * c)	-- ...times Columbs to meters (m)
	q = q * dx[1] * dx[2] * dx[3]					-- ... per meter cubed (1/m^2)
	local I = q	-- amps (C / s) = 1/(m^2 s)
	I = I * c	-- 1/m^3
	local divs = 8 * ni
	for i=0,divs-1 do
		local frac = i / divs
		local th = 2 * math.pi * frac
		local x = math.floor(.5 * n[1] + .25 * ni * math.cos(th))
		local y = math.floor(.5 * n[2] + .25 * ni * math.sin(th))
		local z = math.floor(.5 * n[3]) 
		JU[x+1][y+1][z+1][1] = q
		JU[x+1][y+1][z+1][2] = -I * math.sin(th)
		JU[x+1][y+1][z+1][3] = I * math.cos(th)
		JU[x+1][y+1][z+1][4] = 0
	end
end

local AU = 
--conjgrad
--conjres
gmres
{
	x = -JU,	-- initial AU
	b = JU,
	A = function(AU)
		return matrix.lambda(JU:size(), function(i,j,k,l)
			if i==1 or i==n[1]
			or j==1 or j==n[2]
			or k==1 or k==n[3]
			then
				return 0
			else
				return (AU[i+1][j][k][l]
						+ AU[i-1][j][k][l]
						+ AU[i][j+1][k][l]
						+ AU[i][j-1][k][l]
						+ AU[i][j][k+1][l]
						+ AU[i][j][k-1][l]
						- 6 * AU[i][j][k][l]) / h2
			end
		end)
	end,
	clone = matrix,
	dot = matrix.dot,
	errorCallback = function(err, iter, x, rSq, bSq)
		-- err varies from algorithm to algorithm ... hmm ... maybe it shouldn't ... 
		print(iter, err)
	end,
	epsilon = 1e-5,
	maxiter = 4 * ni^3,
	restart = 100,	--gmres-only
}

local file = assert(io.open('out.txt', 'wb'))
file:write('#x y z J^t J^x J^y J^z A^t A^x A^y A^z\n')
for k=1,n[3] do
	for j=1,n[2] do
		for i=1,n[1] do
			for l=1,4 do
				file:write(i-1,'\t',j-1,'\t',k-1,
					'\t',JU[i][j][k][1],
					'\t',JU[i][j][k][2],
					'\t',JU[i][j][k][3],
					'\t',JU[i][j][k][4],
					'\t',AU[i][j][k][1],
					'\t',AU[i][j][k][2],
					'\t',AU[i][j][k][3],
					'\t',AU[i][j][k][4],
					'\n')
			end
		end
	end
end
file:close()
