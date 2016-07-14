#!/usr/bin/env luajit
require 'ext'
local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'ffi.OpenGL'
local sdl = require 'ffi.sdl'
local ig = require 'ffi.imgui'
local ImGuiApp = require 'imguiapp'
local Mouse = require 'gui.mouse'
local quat = require 'vec.quat'
local vec4 = require 'vec.vec4'
local vec4d = require 'ffi.vec.vec4d'
local glreport = require 'gl.report'
local GradientTex = require 'gl.gradienttex'
local Tex2D = require 'gl.tex2d'
local Tex3D = require 'gl.tex3d'
local Program = require 'gl.program'

local App = class(ImGuiApp)

local mouse = Mouse()
local viewAngle = quat()
local viewDist = 2
	
local texForCol = {}
local hsvTex
local volumeShader

-- which column to render?
local col = ffi.new('int[1]', 4)
local colmax	-- max # of columns
local colnames

local clipEnabled = ffi.new('bool[1]', true)
local rotateClip = ffi.new('int[1]', 0)

local clipInfos = range(4):map(function(i)
	local plane = vec4(0,0,0,0)
	plane[math.min(i,3)] = -1
	return {
		enabled = ffi.new('bool[1]', i==3),
		plane = plane,
	}
end)

local alpha = ffi.new('float[1]', 1.5e-1)
local alphaGamma = ffi.new('float[1]', 1)

function App:initGL()
	App.super.initGL(self)
	self.pts = table()
	self.min = table()
	self.max = table()
	
	for l in io.lines'out.txt' do
		l = l:trim()
		if #l > 0 then
			if l:sub(1,1) == '#' then
				if not colnames then colnames = l:sub(2):trim():split'%s+' end
			else
				-- TODO pick out column titles from first line that starts with '#'
				w = l:split'%s+'
				if not colmax then
					colmax = #w
				else 
					assert(#w == colmax)
				end
				local pt = w:map(function(x) 
					return tonumber(x) or error("expected a number but got "..x)
				end)
				self.pts:insert(pt)
				for i=1,#pt do
					self.min[i] = math.min(self.min[i] or pt[i], pt[i])
					self.max[i] = math.max(self.max[i] or pt[i], pt[i])
				end
			end
		end
	end

--[[
local half = math.floor(self.max[1]/2)
for i=half,self.max[1] do
	local index = i + (self.max[1]+1) * (half + (self.max[2]+1) * half) + 1 
	local pt = self.pts[index]
	print(pt[1], pt[7])
end
--]]

	-- this assumes our points are sequential (and that they are power-of-two?)
	local data = ffi.new('unsigned char[?]', #self.pts*4)
	for col=4,colmax do
		for i,pt in ipairs(self.pts) do
			local f = (pt[col] - self.min[col]) / (self.max[col] - self.min[col])
			data[0+4*(i-1)] = 255*math.clamp(f,0,1)
			data[1+4*(i-1)] = 127
			data[2+4*(i-1)] = 255*math.clamp(1-f,0,1)
			data[3+4*(i-1)] = 255
		end
		local tex = Tex3D{
			width = self.max[1]+1,
			height = self.max[2]+1,
			depth = self.max[3]+1,
			internalFormat = gl.GL_RGBA,
			format = gl.GL_RGBA,
			type = gl.GL_UNSIGNED_BYTE,
			data = data,
			minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
			magFilter = gl.GL_LINEAR,
			wrap = {
				s = gl.GL_CLAMP_TO_EDGE,
				t = gl.GL_CLAMP_TO_EDGE,
				r = gl.GL_CLAMP_TO_EDGE,
			},
			generateMipmap = true,
		}
		texForCol[col] = tex
	end 

	local hsvWidth = 256
	hsvTex = GradientTex(hsvWidth,
--[[ rainbow or heatmap or whatever
		{
			{0,0,0,0},
			{1,0,0,1/6},
			{1,1,0,2/6},
			{0,1,1,3/6},
			{0,0,1,4/6},
			{1,0,1,5/6},
			{0,0,0,6/6},
		},
--]]
-- [[ sunset pic from https://blog.graphiq.com/finding-the-right-color-palettes-for-data-visualizations-fcd4e707a283#.inyxk2q43
		table{
			vec3(22,31,86),
			vec3(34,54,152),
			vec3(87,49,108),
			vec3(156,48,72),
			vec3(220,60,57),
			vec3(254,96,50),
			vec3(255,188,46),
			vec3(255,255,55),
		}:map(function(c,i)
			return table(c/255):append{1}
		end),
--]]
		false)
	-- change to 2D so imgui can use it
	local data = ffi.new('unsigned char[?]', hsvWidth*4)
	gl.glGetTexImage(gl.GL_TEXTURE_1D, 0, gl.GL_RGBA, gl.GL_UNSIGNED_BYTE, data)
	hsvTex:unbind()
	hsvTex:delete()
	hsvTex = Tex2D{
		internalFormat = gl.GL_RGBA,
		width = hsvWidth,
		height = 1,
		format = gl.GL_RGBA,
		type = gl.GL_UNSIGNED_BYTE,
		data = data,
		minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
		magFilter = gl.GL_LINEAR,
		wrap = {
			s = gl.GL_CLAMP_TO_EDGE,
			t = gl.GL_REPEAT,
		},
		generateMipmap = true,
	}

	volumeShader = Program{
		vertexCode = [[
varying vec3 pos;
void main() {
	pos = gl_Vertex.xyz;
	gl_Position = ftransform();
}
]],
		fragmentCode = [[
varying vec3 pos;
uniform sampler3D volTex;
uniform sampler2D hsvTex;
uniform vec3 normal;
uniform float alpha;
uniform float alphaGamma;
void main() {
	
	float value = texture3D(volTex, pos).r;
	vec4 voxelColor = vec4(texture2D(hsvTex, vec2(value, .5)).rgb, pow(alpha, alphaGamma));
	
	//calculate normal in screen coordinates
	vec4 n = gl_ModelViewProjectionMatrix * vec4(normal, 0.);
	//determine length of line through slice at its angle
	voxelColor.a /= -n.w;
	
	gl_FragColor = vec4(voxelColor.rgb, voxelColor.a * alpha);
}
]],
		uniforms = {
			'volTex',
			'hsvTex',
			'normal',
			'alpha',
			'alphaGamma',
		},
	}
	volumeShader:use()
	gl.glUniform1i(volumeShader.uniforms.volTex, 0)
	gl.glUniform1i(volumeShader.uniforms.hsvTex, 1)
	volumeShader:useNone()

	glreport'here'

	print('min',self.min:unpack())
	print('max',self.max:unpack())

	gl.glEnable(gl.GL_DEPTH_TEST)
end

local leftShiftDown
local rightShiftDown 
local imguiCapturing
function App:event(event, eventPtr)
	App.super.event(self, event, eventPtr)
	imguiCapturing = ig.igGetIO()[0].WantCaptureKeyboard
	if imguiCapturing then return end
	
	if event.type == sdl.SDL_MOUSEBUTTONDOWN then
		if event.button.button == sdl.SDL_BUTTON_WHEELUP then
			orbitTargetDistance = orbitTargetDistance * orbitZoomFactor
		elseif event.button.button == sdl.SDL_BUTTON_WHEELDOWN then
			orbitTargetDistance = orbitTargetDistance / orbitZoomFactor
		end
	elseif event.type == sdl.SDL_KEYDOWN or event.type == sdl.SDL_KEYUP then
		if event.key.keysym.sym == sdl.SDLK_LSHIFT then
			leftShiftDown = event.type == sdl.SDL_KEYDOWN
		elseif event.key.keysym.sym == sdl.SDLK_RSHIFT then
			rightShiftDown = event.type == sdl.SDL_KEYDOWN
		elseif event.key.keysym.sym == sdl.SDLK_UP and event.type == sdl.SDL_KEYUP then
			col[0] = math.max(4, col[0] - 1)
		elseif event.key.keysym.sym == sdl.SDLK_DOWN and event.type == sdl.SDL_KEYDOWN then
			col[0] = math.min(colmax, col[0] + 1)
		end
	end
end

function App:update()
	if not imguiCapturing then 
		mouse:update()
	end
	if mouse.leftDragging then
		if leftShiftDown or rightShiftDown then
			if rotateClip[0] == 0 then
				viewDist = viewDist * math.exp(10 * mouse.deltaPos[2])
			else
				local clipPlane = clipInfos[rotateClip[0]].plane
				clipPlane[4] = clipPlane[4] - mouse.deltaPos[2]
			end
		else
			local magn = mouse.deltaPos:length() * 1000
			if magn > 0 then
				local normDelta = mouse.deltaPos / magn
				local r = quat():fromAngleAxis(-normDelta[2], normDelta[1], 0, -magn)
				if rotateClip[0] == 0 then
					viewAngle = (viewAngle * r):normalize()
				else
					local clipPlane = clipInfos[rotateClip[0]].plane
					local clipNormal = (viewAngle * r * viewAngle:conjugate()):conjugate():rotate(vec3(clipPlane:unpack()))
					for i=1,3 do
						clipPlane[i] = clipNormal[i]
					end
				end
			end
		end
	end
	
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	local w, h = self:size()
	local ar = w / h
	local znear, zfar = .1, 100
	gl.glFrustum(-ar*znear, ar*znear, -znear, znear, znear, zfar)

	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	
	gl.glTranslatef(0, 0, -viewDist)

	local aa = viewAngle:toAngleAxis()
	gl.glRotated(-aa[4], aa[1], aa[2], aa[3])

	for i,clipInfo in ipairs(clipInfos) do
		gl.glClipPlane(gl.GL_CLIP_PLANE0+i-1, vec4d(clipInfo.plane:unpack()):ptr())
		if clipInfo.enabled[0] then 
			gl.glEnable(gl.GL_CLIP_PLANE0+i-1)
		end
	end

	gl.glTranslatef(-.5, -.5, -.5)

	local tex = texForCol[col[0]]

	volumeShader:use()
	tex:bind(0)
	hsvTex:bind(1)
	gl.glUniform1f(volumeShader.uniforms.alpha, alpha[0])
	gl.glUniform1f(volumeShader.uniforms.alphaGamma, alphaGamma[0])

	gl.glEnable(gl.GL_TEXTURE_GEN_S)
	gl.glEnable(gl.GL_TEXTURE_GEN_T)
	gl.glEnable(gl.GL_TEXTURE_GEN_R)
	gl.glTexGeni(gl.GL_S, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGeni(gl.GL_T, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGeni(gl.GL_R, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGendv(gl.GL_S, gl.GL_OBJECT_PLANE, vec4d(1,0,0,0):ptr())
	gl.glTexGendv(gl.GL_T, gl.GL_OBJECT_PLANE, vec4d(0,1,0,0):ptr())
	gl.glTexGendv(gl.GL_R, gl.GL_OBJECT_PLANE, vec4d(0,0,1,0):ptr())

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glEnable(gl.GL_BLEND)

	--[[ points
	gl.glPointSize(2)
	gl.glBegin(gl.GL_POINTS)
	for _,pt in ipairs(self.pts) do
		gl.glVertex3d( 
			(pt[1] - .5)/(self.max[1] + 1),
			(pt[2] - .5)/(self.max[2] + 1),
			(pt[3] - .5)/(self.max[3] + 1))
	end
	gl.glEnd()
	--]]
	-- [[ slices
	local n = 255
	local fwd = -viewAngle:zAxis()
	local fwddir = select(2, table(fwd):map(math.abs):sup())
	local quad = {{0,0},{1,0},{1,1},{0,1}}
	local jmin, jmax, jdir
	if fwd[fwddir] < 0 then
		jmin, jmax, jdir = 0, n, 1
	else
		jmin, jmax, jdir = n, 0, -1
	end
	gl.glUniform3f(volumeShader.uniforms.normal, fwddir==1 and jdir or 0, fwddir==2 and jdir or 0, fwddir==3 and jdir or 0)
	
	gl.glBegin(gl.GL_QUADS)
	for j=jmin,jmax,jdir do
		local f = j/n
		for _,vtx in ipairs(quad) do
			if fwddir == 1 then
				gl.glVertex3f(f, vtx[1], vtx[2])
			elseif fwddir == 2 then
				gl.glVertex3f(vtx[1], f, vtx[2])
			elseif fwddir == 3 then
				gl.glVertex3f(vtx[1], vtx[2], f)
			end
		end
	end
	gl.glEnd()
	--]]
	
	gl.glDisable(gl.GL_BLEND)

	gl.glDisable(gl.GL_TEXTURE_GEN_S)
	gl.glDisable(gl.GL_TEXTURE_GEN_T)
	gl.glDisable(gl.GL_TEXTURE_GEN_R)

	hsvTex:unbind(1)
	tex:unbind(0)
	volumeShader:useNone()

	for i,clipInfo in ipairs(clipInfos) do
		gl.glDisable(gl.GL_CLIP_PLANE0+i-1)
	end
	
	glreport'here'

	App.super.update(self)
end

function App:updateGUI()
	col[0] = col[0] - 4
	ig.igCombo('column', col, colnames:sub(4))
	col[0] = col[0] + 4
	ig.igText(self.min[col[0]]..' to '..self.max[col[0]])
	
	local gradImageSize = ig.ImVec2(128, 32)
	ig.igImage(
		ffi.cast('void*',ffi.cast('intptr_t',hsvTex.id)),
		gradImageSize)
	local gradScreenPos = ig.igGetCursorScreenPos()
	local mousePos = ig.igGetMousePos()
	local cursorX = mousePos.x - gradScreenPos.x
	local cursorY = gradScreenPos.y - mousePos.y
	if cursorX >= 0 and cursorX <= gradImageSize.x
	and cursorY >= 0 and cursorY <= gradImageSize.y
	then
		local frac = cursorX / gradImageSize.x
		ig.igBeginTooltip()
		ig.igText(tostring( self.min[col[0]] * (1-frac) + self.max[col[0]] * frac ))
		ig.igEndTooltip()
	end
	
	ig.igSliderFloat('alpha', alpha, 0, 1, '%.3e', 10)
	ig.igSliderFloat('gamma', alphaGamma, 0, 1000, '%.3e', 10)
	ig.igRadioButton("rotate camera", rotateClip, 0)
	for i,clipInfo in ipairs(clipInfos) do
		ig.igPushIdStr('clip '..i)
		ig.igCheckbox('clip', clipInfo.enabled)
		ig.igSameLine()
		ig.igRadioButton('rotate', rotateClip, i)
		ig.igPopId()
	end
end

local app = App()
app:run()
