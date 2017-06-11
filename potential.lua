#!/usr/bin/env luajit
require 'ext'
local gnuplot = require 'gnuplot'

local n = 4000
local I = function(i) return i end
local C = function(c) return function() return c end end
local U = function(f) return range(n):map(f) end
local data = function(t) return table.map(t, U) end

local c = 299792458						--  ... m/s = 1
local G = 6.67384e-11						--  ... m^3 / (kg s^2) = 1
local kg_in_m = G / c / c

local theta = math.pi/2						--  polar angle for Kerr metric

local R = 100
local M = 1
local E = .5	-- energy of the particle
local l = 4.1
local function schwarzschild_V(r,l)
	return math.sqrt((1 - 2 * M / r) * (1 + l^2 / r^2))
end
local function schwarzschild_V_l(r)
	return schwarzschild_V(r,l)
end

local r = (I-.5)/n * R*10
gnuplot{
	output = 'images/schwarzschild potential.png',
	style = 'data lines',
	data = data{r, schwarzschild_V_l:o(r)},
	xrange = {1, R},
	yrange = {.9, 1.05},
	log = 'x',
	xlabel = 'log(r/M)',
	ylabel = 'V(r)',
	{using='1:2', title='schwarzschild potential'},
}

local m = 100
local lmax = 5
local ls = (I-.5)/m * lmax
gnuplot{
	output = 'images/schwarzschild potential 2.png',
	style = 'data lines',
	xlabel = 'r/M',
	ylabel = 'l/M',
	zlabel = 'V(r)',
	contour = true,
	xrange = {1, R},
	yrange = {0, lmax},
	zrange = {.8, 1.2},
	griddata = {
		x=range(n):map(r), 
		y=range(m):map(ls), 
		range(n):map(function(i)
			return range(m):map(function(j)
				return schwarzschild_V(r(i),ls(j))
			end)
		end),
	},
	{splot=true, using='1:2:3', title='schwarzschild potential'},
}
