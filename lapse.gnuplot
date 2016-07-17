#!/opt/local/bin/gnuplot -persist

sphereVolume(r) = 4./3.*pi*r*r*r
sphereSurface(r) = 4.*pi*r*r
min(x,y) = (x < y) ? x : y

c = 299792458.
G = 6.67384e-11
# earth parameters:
#R = 6.37101e+6						# matter radius
#mass = 5.9736e+24 * G / c / c			# total mass within radius
#density = mass / sphereVolume(R)	# average density 
# temp:
R = 1.
density = .1

# matter within radius
m(x) = density * sphereVolume(min(R,x))
# r^1 gives a constant line within the planet surface
# r^2 gives a line from the surface alpha to alpha=1 at r=0
# r^3 gives a parabola, so zero derivative at r=0
#m(x) = .1 * min(R, x)**3.
# derivative of schwarzschild radius function wrt radius
dm_dr(x) = R < x ? 0 : (density * sphereSurface(x))

# total matter
M = m(R)

# lapse 
# lapse without considering hydrostatic pressure
#alpha(r) = sqrt(1. - 2.*m(r)/r)
# lapse with considerations (MTW box 23.2 eqn 6)
alpha(r) = r > R ? sqrt(1. - 2.*M/r) : (1.5*sqrt(1. - 2.*M/R) - .5*sqrt(1. - 2.*M*r*r/(R*R*R)))

# derivative of lapse wrt radius
dalpha_dr(r) = r > R ? (M / (r * r * sqrt(1. - 2.*M/r))) : M*r/(R*R*R*sqrt(1. - 2.*M*r*r/(R*R*R)))

# metric, based on g_tt = -alpha^2 + beta^2 and assuming beta=0
g_tt(r) = -alpha(r)**2

dx = R/1000.
plot [0:R*10.] alpha(x) title 'alpha', dalpha_dr(x) title 'dalpha/dr', (alpha(x+dx) - alpha(x-dx)) / (2.*dx) title 'dalpha/dr discrete', g_tt(x) title 'g_t_t'
#plot [0:R*10.] log(alpha(x))
