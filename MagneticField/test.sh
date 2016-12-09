#!/usr/bin/env sh
#dist/osx/debug/test.app/Contents/MacOS/test discreteLaplacian > out.txt
lua -lmake && dist/linux/debug/MagneticField > out.txt 2> converge.txt && gnuplot -p -e "set term wx; set style data lines; splot 'out.txt' using 1:2:3; pause -1"
