#!/opt/local/bin/gnuplot -p

set terminal wx
a = 5

b(x,y,z) = a**2 - x**2 - y**2 - z**2

discr(x,y,z) = b(x,y,z)**2 - 4.*(-a**2 * z**2)
# discr < 0 for a^2 (1 + 4 z^2) < x^2 + y^2 + z^2

rSqPlus(x,y,z) = (-b(x,y,z) + sqrt(discr(x,y,z)))/2.
rSqMinus(x,y,z) = (-b(x,y,z) + sqrt(discr(x,y,z)))/2.

set isosamples 50

# 1D plots at fixed z ...
z = 0
#plot [5:6] sqrt(rSqMinus(x,0,z)) title 'r-'
#plot [5:6] sqrt(rSqPlus(x,0,z)) title 'r+'
#plot [5:6] sqrt(rSqPlus(x,0,z)) title 'r+', sqrt(rSqMinus(x,0,z)) title 'r-'
#plot [0:10] discr(x,0,1)
splot [0:20][-10:10] sqrt(rSqPlus(x,y,z)) title 'r+', sqrt(rSqMinus(x,y,z)) title 'r-'

pause -1
