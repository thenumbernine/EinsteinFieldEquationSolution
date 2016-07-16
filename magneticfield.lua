#!/usr/bin/env luajit
local ni = ... and tonumber(...) or 8
local jfnk = require 'jfnk'
local matrix = require 'matrix'

local n = {ni,ni,ni}
local h2 = 1/ni^2

local M = 1e+6
local rho = matrix.lambda(n, function(i,j,k)
	return (math.abs(i - (n[1]+1)/2) < 1
		and math.abs(j - (n[2]+1)/2) < 1
		and math.abs(k - (n[3]+1)/2) < 1)
		and (
			i >= (n[1]+1)/2 and M or -M
		) or 0
end)

local phi = jfnk{
	x = rho,
	f = function(phi)
		return matrix.lambda(n, function(i,j,k)
			if i==1 or i==n[1]
			or j==1 or j==n[2]
			or k==1 or k==n[3]
			then
				return 0
			else
				if (math.abs(i - (n[1]+1)/2) < 1
					and math.abs(j - (n[2]+1)/2) < 1
					and math.abs(k - (n[3]+1)/2) < 1)
				then
					return rho[i][j][k]
				else
					return (phi[i+1][j][k]
						+ phi[i-1][j][k]
						+ phi[i][j+1][k]
						+ phi[i][j-1][k]
						+ phi[i][j][k+1]
						+ phi[i][j][k-1]
						- 6 * phi[i][j][k]) / (-h2 * 6 * math.pi)
				end
			end
		end)
	end,
	dot = function(a,b)
		local sum = 0
		for i=1,n[1] do
			for j=1,n[2] do
				for k=1,n[3] do
					sum = sum + a[i][j][k] * b[i][j][k]
				end
			end
		end
		return sum
	end,
}

print('#x y z rho phi')
for i=1,n[1] do
	for j=1,n[2] do
		for k=1,n[3] do
			print(i-1,j-1,k-1,rho[i][j][k],phi[i][j][k])
		end
	end
end
