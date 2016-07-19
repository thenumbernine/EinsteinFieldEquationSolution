#!/opt/local/bin/gnuplot -p

sphereVolume(r) = 4./3.*pi*r*r*r
sphereSurface(r) = 4.*pi*r*r
min(x,y) = (x < y) ? x : y

c = 299792458.
G = 6.67384e-11
# earth parameters:
R = 6.37101e+6						# matter radius
mass = 5.9736e+24 * G / c / c		# total mass within radius
density = mass / sphereVolume(R)	# average density 
# temp:
#R = 1.
#density = .1

# matter within radius
m(x) = density * sphereVolume(min(R,x))
# r^1 gives a constant line within the planet surface
# r^2 gives a line from the surface alpha to alpha=1 at r=0
# r^3 gives a parabola, so zero derivative at r=0
#m(x) = .1 * min(R, x)**3.
# derivative of schwarzschild radius function wrt radius
dm_dr(x) = R < x ? 0. : (density * sphereSurface(x))

# total matter
M = m(R)

# lapse 
# lapse without considering hydrostatic pressure (i.e. if m(r) = m(radius), i.e. if radius=0, i.e. a black hole)
#schwarzschild_alpha(r) = sqrt(1. - 2.*m(r)/r)
# lapse with considerations (MTW box 23.2 eqn 6)
schwarzschild_alpha(r) = r > R ? sqrt(1. - 2.*M/r) : (1.5*sqrt(1. - 2.*M/R) - .5*sqrt(1. - 2.*M*r*r/(R*R*R)))

# derivative of lapse wrt radius
schwarzschild_dr_alpha(r) = r > R ? (M / (r * r * sqrt(1. - 2.*M/r))) : M*r/(R*R*R*sqrt(1. - 2.*M*r*r/(R*R*R)))

# metric, based on g_tt = -alpha^2 + beta^2 and assuming beta=0
schwarzschild_g_tt(r) = -schwarzschild_alpha(r)**2

# now for kerr.
# but now the theta and phi coordinates matter ... 
# Alcubierre "Introduction to 3+1 Numerical Relativity" eqn 3.4.33 thru 3.4.35
#alcubierre_alpha = sqrt(1. / (1. + 2. * H(r)))

#plot [0:R*10.] schwarzschild_alpha(x) title 'schwarzschild alpha', schwarzschild_dr_alpha(x) title 'schwarzschild d_r alpha', schwarzschild_g_tt(x) title 'g_t_t', -schwarzschild_alpha(x) * schwarzschild_dr_alpha(x) title 'schwarzschild gravity = -Gamma_r_t_t = -alpha d_r alpha'
# gravity-only, in m/s^2
plot [0:R*10.] -schwarzschild_alpha(x) * schwarzschild_dr_alpha(x) * c * c title 'schwarzschild gravity = -Gamma_r_t_t = -alpha d_r alpha in m/s^2'

pause -1
