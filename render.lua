#!/usr/bin/env luajit
--[[
this should be a stand-alone tool
--]]
require 'ext'
local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'gl'
local sdl = require 'sdl'
local ig = require 'imgui'
local ImGuiApp = require 'imgui.app'
local Mouse = require 'glapp.mouse'
local quat = require 'vec.quat'
local vec3 = require 'vec.vec3'
local vec4 = require 'vec.vec4'
local vec4d = require 'vec-ffi.vec4d'
local glreport = require 'gl.report'
local glcall = require 'gl.call'
local GradientTex = require 'gl.gradienttex'
local Tex2D = require 'gl.tex2d'
local Tex3D = require 'gl.tex3d'
local Program = require 'gl.program'

local usePoints = false
if _G.usePoints ~= nil then usePoints = _G.usePoints end
local useSlices = true
if _G.useSlices ~= nil then useSlices = _G.useSlices end

local filename, col = ...
filename = filename or 'out.txt'

local mouse = Mouse()
local viewAngle = quat()
local viewDist = 2
	
local texForCol = {}
local hsvTex
local volumeShader

-- which column to render?
local col = tonumber(col) or 4
local colmax	-- max # of columns
local colnames

local rotateClip = 0

local function makeDefaultPlane(i)
	local plane = vec4(0,0,0,0)
	plane[math.min(i,3)] = -1
	return plane
end

local clipInfos = range(4):map(function(i)
	return {
		enabled = i==3,
		plane = makeDefaultPlane(i),
	}
end)

-- _G for imgui table access
alpha = 1.5e-1
alphaGamma = 1
showGradTrace = false
showCurlTrace = false
flipGradient = false

local App = ImGuiApp:subclass()
App.viewUseGLMatrixMode= true

function App:initGL()
	App.super.initGL(self)
	self.pts = table()
	self.min = table()
	self.max = table()
	
	for l in io.lines(filename) do
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
	print(pt[1], pt[6])	-- alpha
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
	hsvTex:toCPU(data)
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
uniform bool clipEnabled[4];
uniform bool flipGradient;
void main() {

	vec4 worldPos = gl_ModelViewMatrix * vec4(pos,1.);
	for (int i = 0; i < 4; ++i) {
		if (clipEnabled[i] && dot(worldPos, gl_ClipPlane[i]) < 0.) discard;
	}

	float value = texture3D(volTex, pos).r;
	if (flipGradient) value = 1. - value;
	vec3 texColor = texture2D(hsvTex, vec2(value, .5)).rgb;
	float voxelAlpha = pow(alpha, alphaGamma);
	vec4 voxelColor = vec4(texColor, voxelAlpha);
	
	//calculate normal in screen coordinates
	vec4 n = gl_ModelViewProjectionMatrix * vec4(normal, 0.);
	//determine length of line through slice at its angle
	voxelColor.a /= -n.w;
	
	gl_FragColor = vec4(voxelColor.rgb, voxelColor.a * alpha);
}
]],
		uniforms = {
			volTex = 0,
			hsvTex = 1,
		},
	}:useNone()

	glreport'here'

	for i,name in ipairs(colnames) do
		print(name, self.min[i], self.max[i])
	end

	gl.glEnable(gl.GL_DEPTH_TEST)
end

local leftShiftDown
local rightShiftDown
local imguiCapturing
function App:event(eventPtr)
	App.super.event(self, eventPtr)
	imguiCapturing = ig.igGetIO()[0].WantCaptureKeyboard
	if imguiCapturing then return end
	
	if eventPtr[0].type == sdl.SDL_EVENT_KEY_DOWN
	or eventPtr[0].type == sdl.SDL_EVENT_KEY_UP
	then
		if eventPtr[0].key.key == sdl.SDLK_LSHIFT then
			leftShiftDown = eventPtr[0].type == sdl.SDL_EVENT_KEY_DOWN
		elseif eventPtr[0].key.key == sdl.SDLK_RSHIFT then
			rightShiftDown = eventPtr[0].type == sdl.SDL_EVENT_KEY_DOWN
		elseif eventPtr[0].key.key == sdl.SDLK_UP
		and eventPtr[0].type == sdl.SDL_EVENT_KEY_UP
		then
			col = math.max(4, col - 1)
		elseif eventPtr[0].key.key == sdl.SDLK_DOWN
		and eventPtr[0].type == sdl.SDL_EVENT_KEY_DOWN
		then
			col = math.min(colmax, col + 1)
		end
	end
end

function App:update()
	if not imguiCapturing then
		mouse:update()
	end
	if mouse.leftDragging then
		if leftShiftDown or rightShiftDown then
			if rotateClip == 0 then
				viewDist = viewDist * math.exp(10 * mouse.deltaPos.y)
			else
				local clipPlane = clipInfos[rotateClip].plane
				clipPlane[4] = clipPlane[4] - mouse.deltaPos.y
			end
		else
			local magn = mouse.deltaPos:length() * 1000
			if magn > 0 then
				local normDelta = mouse.deltaPos / magn
				local r = quat():fromAngleAxis(-normDelta.y, normDelta.x, 0, -magn)
				if rotateClip == 0 then
					viewAngle = (viewAngle * r):normalize()
				else
					local clipPlane = clipInfos[rotateClip].plane
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
		gl.glClipPlane(gl.GL_CLIP_PLANE0+i-1, vec4d(clipInfo.plane:unpack()).s)
	end

	gl.glTranslatef(-.5, -.5, -.5)

	local tex = texForCol[col]

	volumeShader:use()
	tex:bind(0)
	hsvTex:bind(1)
	gl.glUniform1f(volumeShader.uniforms.alpha.loc, alpha)
	gl.glUniform1f(volumeShader.uniforms.alphaGamma.loc, alphaGamma)
	gl.glUniform1iv(volumeShader.uniforms['clipEnabled[0]'].loc, 4,
		ffi.new('int[4]', clipInfos:map(function(info) return info.enabled end)))
	gl.glUniform1i(volumeShader.uniforms.flipGradient.loc, flipGradient)

	gl.glEnable(gl.GL_TEXTURE_GEN_S)
	gl.glEnable(gl.GL_TEXTURE_GEN_T)
	gl.glEnable(gl.GL_TEXTURE_GEN_R)
	gl.glTexGeni(gl.GL_S, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGeni(gl.GL_T, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGeni(gl.GL_R, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGendv(gl.GL_S, gl.GL_OBJECT_PLANE, vec4d(1,0,0,0).s)
	gl.glTexGendv(gl.GL_T, gl.GL_OBJECT_PLANE, vec4d(0,1,0,0).s)
	gl.glTexGendv(gl.GL_R, gl.GL_OBJECT_PLANE, vec4d(0,0,1,0).s)

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glEnable(gl.GL_BLEND)

	if usePoints then
		gl.glPointSize(2)
		gl.glBegin(gl.GL_POINTS)
		for _,pt in ipairs(self.pts) do
			gl.glVertex3d(
				(pt[1] - self.min[1])/(self.max[1] - self.max[1]),
				(pt[2] - self.min[2])/(self.max[2] - self.max[2]),
				(pt[3] - self.min[3])/(self.max[3] - self.max[3]))
		end
		gl.glEnd()
	end
	if useSlices then
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
		gl.glUniform3f(volumeShader.uniforms.normal.loc, fwddir==1 and jdir or 0, fwddir==2 and jdir or 0, fwddir==3 and jdir or 0)
		
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
	end
	
	gl.glDisable(gl.GL_BLEND)

	gl.glDisable(gl.GL_TEXTURE_GEN_S)
	gl.glDisable(gl.GL_TEXTURE_GEN_T)
	gl.glDisable(gl.GL_TEXTURE_GEN_R)

	hsvTex:unbind(1)
	tex:unbind(0)
	volumeShader:useNone()

	if showGradTrace or showCurlTrace then
		for i,clipInfo in ipairs(clipInfos) do
-- intel/ubuntu was having trouble when the clip plane included the viewport
-- so I moved the clipping code to the shader
			if clipInfo.enabled then
				gl.glEnable(gl.GL_CLIP_PLANE0+i-1)
			end
		end
	
		
		gl.glDisable(gl.GL_DEPTH_TEST)

		self.gradLists = self.gradLists or {}
		self.gradLists[col] = self.gradLists[col] or {}
		
		self.curlLists = self.curlLists or {}
		self.curlLists[col] = self.curlLists[col] or {}
		
		for curl=0,1 do
			if (curl==0 and showGradTrace) or (curl==1 and showCurlTrace) then
				glcall(curl==0 and self.gradLists[col] or self.curlLists[col], function()
					-- TODO
					-- 1) space out seeds
					-- 2) integrate rk4 or something the lines
					local divs = 16
					for i=1,self.max[1]-1,(self.max[1]+1)/divs do
						for j=1,self.max[2]-1,(self.max[2]+1)/divs do
							for k=1,self.max[3]-1,(self.max[3]+1)/divs do
								local function trace(x, dir)
									gl.glBegin(gl.GL_LINE_STRIP)
									local lastdx
									local lastT
									while true do
										local function intx(x)
											local i = vec3(math.floor(x[1]), math.floor(x[2]), math.floor(x[3]))
											if i[1] < 1 or i[1] > self.max[1]-2
											or i[2] < 1 or i[2] > self.max[2]-2
											or i[3] < 1 or i[3] > self.max[3]-2
											then
												return
											end
											return i
										end
										local dx_ds = function(x)
											local i = intx(x)
											if not i then return end
											local f = x - i
											
											local indexm = 1 + i[1] + (self.max[1] + 1) * (i[2] + (self.max[2] + 1) * i[3])
											assert(indexm >= 1 and indexm <= #self.pts, "failed for i="..i.." x="..x)
											local dxm = vec3(
												self.pts[indexm + 1][col] - self.pts[indexm - 1][col],
												self.pts[indexm + self.max[1] + 1][col] - self.pts[indexm - self.max[1] - 1][col],
												self.pts[indexm+(self.max[1]+1)*(self.max[2]+1)][col]-self.pts[indexm-(self.max[1]+1)*(self.max[2]+1)][col])
											local indexp = 1 + i[1]+1 + (self.max[1] + 1) * (i[2]+1 + (self.max[2] + 1) * (i[3]+1))
											assert(indexp >= 1 and indexp <= #self.pts, "failed for i="..i.." x="..x)
											local dxp = vec3(
												self.pts[indexp + 1][col] - self.pts[indexp - 1][col],
												self.pts[indexp + self.max[1] + 1][col] - self.pts[indexp - self.max[1] - 1][col],
												self.pts[indexp+(self.max[1]+1)*(self.max[2]+1)][col]-self.pts[indexp-(self.max[1]+1)*(self.max[2]+1)][col])
									
											-- one vector field is along the gradient ... the B field ...
											-- the E field goes around ... normal to the plane of curvature of the B field? Frenet frame
											-- the S vector is E x B

											local dx = vec3(
												dxm[1] * (1 - f[1]) + dxp[1] * f[1],
												dxm[2] * (1 - f[2]) + dxp[2] * f[2],
												dxm[3] * (1 - f[3]) + dxp[3] * f[3])
											local len = dx:length()
											if len < 1e-20 then return end
											return dx / len
										end
										
										local h = 1/(self.max[1]+1)
										
										if curl==1 then
											local T_of_x = dx_ds
											dx_ds = function(x)
												local T = T_of_x(x)
												if not T then return lastT end
												local Tp = T_of_x(x + T*h)
												if not Tp then return lastT end
												local Tm = T_of_x(x - T*h)
												if not Tm then return lastT end
												local N = Tm:cross(Tp)
												local len = N:length()
												if len < 1e-20 then return lastT end
												return N / len
											end
										end

										-- [[ euler
										local dx = dx_ds(x)
										if not dx then break end
										if curl==1 then lastT = dx end
										dx = dx * h
										--]]
										--[[ rk4
										local k1 = dx_ds(x)
										if not k1 then break end
										local k2 = dx_ds(x + k1 * (h/2))
										if not k2 then break end
										local k3 = dx_ds(x + k2 * (h/2))
										if not k3 then break end
										local k4 = dx_ds(x + k3 * h)
										if not k4 then break end
										local dx = k1*(h/6) + k2*(h/3) + k3*(h/3) + k4*(h/6)
										--]]
										if lastdx and dx:dot(lastdx) < 1e-20 then break end
										lastdx = dx
										x = x + dx * dir
										if not intx(x) then break end
										if not math.isfinite(x[1]) or not math.isfinite(x[2]) or not math.isfinite(x[3]) then break end

										gl.glVertex3d(
											(x[1]+.5)/(self.max[1]+1),
											(x[2]+.5)/(self.max[2]+1),
											(x[3]+.5)/(self.max[3]+1))
									end
									gl.glEnd()
								end

								local pt = vec3(i+.5,j+.5,k+.5)
								trace(pt, 1)
								trace(pt, -1)
							end
						end
					end
				end)
			end
		end
		gl.glEnable(gl.GL_DEPTH_TEST)
	
		for i,clipInfo in ipairs(clipInfos) do
			gl.glDisable(gl.GL_CLIP_PLANE0+i-1)
		end
	end

	glreport'here'

	App.super.update(self)
end

print'App:updateGUI() FIXME'
function App:updateGUI()
do return end -- TODO FIXME lots in here is broken
	col = col - 4
	ig.luatableCombo('column', _G, 'col', colnames and colnames:sub(4) or range(4,colmax))
	col = col + 4
	ig.igText(('%.3e to %.3e'):format(self.min[col], self.max[col]))
	
	local gradImageSize = ig.ImVec2(128, 32)
	ig.igImage(
		ffi.cast('void*',ffi.cast('intptr_t',hsvTex.id)),
		gradImageSize,
		ig.ImVec2(flipGradient and 1 or 0,0),
		ig.ImVec2(flipGradient and 0 or 1,1))
	local gradScreenPos = ig.igGetCursorScreenPos()
	local mousePos = ig.igGetMousePos()
	local cursorX = mousePos.x - gradScreenPos.x
	local cursorY = gradScreenPos.y - mousePos.y
	if cursorX >= 0 and cursorX <= gradImageSize.x
	and cursorY >= 0 and cursorY <= gradImageSize.y
	then
		local frac = cursorX / gradImageSize.x
		ig.igBeginTooltip()
		ig.igText(tostring( self.min[col] * (1-frac) + self.max[col] * frac ))
		ig.igEndTooltip()
	end
	
	ig.luatableSliderFloat('alpha', _G, 'alpha', 0, 1, '%.3e', 10)
	ig.luatableSliderFloat('gamma', _G, 'alphaGamma', 0, 1000, '%.3e', 10)
	ig.luatableCheckbox('flip gradient', _G, 'flipGradient')
	ig.luatableRadioButton("rotate camera", _G, 'rotateClip', 0)
	for i,clipInfo in ipairs(clipInfos) do
		ig.igPushID_Str('clip '..i)
		ig.luatableCheckbox('clip', clipInfo, 'enabled')
		ig.igSameLine()
		ig.luatableRadioButton('rotate', _G, 'rotateClip', i)
		ig.igSameLine()
		if ig.igButton('reset') then
			clipInfo.plane = makeDefaultPlane(i)
		end
		ig.igPopId()
	end
	ig.luatableCheckbox('show gradient trace', _G, 'showGradTrace')
	ig.luatableCheckbox('show curl trace', _G, 'showCurlTrace')
end

local app = App()
app:run()
