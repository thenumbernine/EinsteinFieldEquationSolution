#!/usr/bin/env sh
#dist/osx/debug/test.app/Contents/MacOS/test discreteLaplacian > out.txt
lua -lmake && dist/linux/debug/MagneticField 2> converge.txt && ../render.lua 
