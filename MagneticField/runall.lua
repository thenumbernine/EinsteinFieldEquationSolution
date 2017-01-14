#!/usr/bin/env lua
local template = require 'template'
for _,name in ipairs{'conjgrad', 'conjres', 'gmres'} do
	print('creating '..template('out.cpu.lua.<?=name?>.txt', {name=name}))
	os.execute(template('./magneticfield.lua 16 <?=name?> > out.cpu.lua.<?=name?>.txt', {name=name}))
end
for _,name in ipairs{'conjgrad', 'conjres', 'gmres'} do
	print('creating '..template('out.gpu.lua.<?=name?>.txt', {name=name}))
	os.execute(template('./magneticfield-cl-krylov.lua 16 <?=name?> > out.gpu.lua.<?=name?>.txt', {name=name}))
end
