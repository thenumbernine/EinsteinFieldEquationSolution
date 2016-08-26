#!/opt/local/bin/gnuplot -p

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
