#!/opt/local/bin/gnuplot -p

sphereVolume(r) = 4./3.*pi*r*r*r
sphereSurface(r) = 4.*pi*r*r
min(x,y) = (x < y) ? x : y

c = 299792458.						# ... m/s = 1
G = 6.67384e-11						# ... m^3 / (kg s^2) = 1
# earth parameters:
R = 6.37101e+6						# matter radius, in m
mass = 5.9736e+24 * G / c / c		# total mass within radius, in m
density = mass / sphereVolume(R)	# average density, in m^-2
# temporary:
#R = 1.
#density = .1

m(x) = density * sphereVolume(min(R,x))		# matter within radius, in m
# r^1 gives a constant line within the planet surface
# r^2 gives a line from the surface alpha to alpha=1 at r=0
# r^3 gives a parabola, so zero derivative at r=0
#m(x) = .1 * min(R, x)**3.
# derivative of schwarzschild radius function wrt radius
dm_dr(x) = R < x ? 0. : (density * sphereSurface(x))

M = m(R)	# matter within planet radius, in m

# lapse 
# lapse without considering hydrostatic pressure (i.e. if m(r) = m(radius), i.e. if radius=0, i.e. a black hole)
#schwarzschild_alpha(r) = sqrt(1. - 2.*m(r)/r)
# lapse with considerations (MTW box 23.2 eqn 6)
schwarzschild_alpha(r) = r > R ? sqrt(1. - 2.*M/r) : (1.5*sqrt(1. - 2.*M/R) - .5*sqrt(1. - 2.*M*r*r/(R*R*R)))

# derivative of lapse wrt radius
schwarzschild_dr_alpha(r) = r > R ? (M / (r * r * sqrt(1. - 2.*M/r))) : M*r/(R*R*R*sqrt(1. - 2.*M*r*r/(R*R*R)))

# metric, based on g_tt = -alpha^2 + beta^2 and assuming beta=0
schwarzschild_g_tt(r) = -schwarzschild_alpha(r)**2

# schwarzschild gravity = Gamma_rtt = -alpha d/dr alpha ... times c^2 for m/s^2
schwarzschild_gravity(r) = -schwarzschild_alpha(r) * schwarzschild_dr_alpha(r) * c**2

# newton gravity ... times c^2 for m/s^2
newton_gravity(r) = -m(r) / r**2 * c**2

# now for kerr.
# I'll just compute the gravity
# there is no solution yet solved for the equations of structure inside a star using the Kerr metric
theta = pi/2.										# polar angle
angularVelocity = 2. * pi / (60. * 60. * 24.) / c	# angular velocity, in m^-1
inertia = 2. / 5. * M * R**2						# moment of inertia about a sphere, in m^3
angularMomentum = inertia * angularVelocity			# angular momentum in m^2
earth_a = angularMomentum / M						# in m

Delta(r,a) = r**2 - 2.*m(r) * r + a**2
Sigma(r,a) = r**2 + a**2 * cos(theta)**2
kerr_gravity(r,a) = -2.*m(r) * Delta(r,a) * (r**2 - a**2 * cos(theta)**2) / (2 * Sigma(r,a)**3) * c**2

# show how schwarzschild alpha, d/dr alpha, g_tt, and Gamma^r_tt are all related:
#plot [0:R*10.] schwarzschild_alpha(x) title 'schwarzschild alpha', schwarzschild_dr_alpha(x) title 'schwarzschild d_r alpha', schwarzschild_g_tt(x) title 'g_t_t', -schwarzschild_alpha(x) * schwarzschild_dr_alpha(x) title 'schwarzschild gravity = -Gamma_r_t_t = -alpha d_r alpha'

# show gravity only, in m/s^2
#set log y
#set ylabel 'm/s^2'
#plot [0:R*10.] -newton_gravity(x) title 'newton gravity', -schwarzschild_gravity(x) title 'schwarzschild gravity', -kerr_gravity(x) title 'kerr gravity'

# show differences in gravity 
#set ylabel 'm/s^2'
#plot [0:R*2.] abs(schwarzschild_gravity(x)) - abs(newton_gravity(x)) title '|schwarzschild|-|newton|', \
#				abs(schwarzschild_gravity(x)) - abs(kerr_gravity(x,earth_a)) title '|schwarzschild|-|kerr|', \
#				abs(newton_gravity(x)) - abs(kerr_gravity(x,earth_a)) title '|newton|-|kerr|'

# rotation-less Kerr isn't the same as Schwarzschild ... It weaker than rotation Kerr.  It is within 2e-10 of rotation Kerr ... 2e-11 at the Earth's surface
# note that Schwarzschild (which is rotation-less) is stronger than rotation-Kerr and rotation-less Kerr by about 1.4e-8
set ylabel 'm/s^2'
plot [0:R*2.] abs(kerr_gravity(x,earth_a)) - abs(kerr_gravity(x,0)) title '|kerr w/rotation|-|kerr w/o rotation|'

pause -1
