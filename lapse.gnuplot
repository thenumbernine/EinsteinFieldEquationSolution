#!/opt/local/bin/gnuplot -p

set terminal png size 800,600

sphereVolume(r) = 4./3.*pi*r*r*r
sphereSurface(r) = 4.*pi*r*r
min(x,y) = (x < y) ? x : y

c = 299792458.						# ... m/s = 1
G = 6.67384e-11						# ... m^3 / (kg s^2) = 1

# set radius and density to temporary values:
R = 1.
density = .1
load 'lapse.inc.gnuplot'

# show how schwarzschild alpha, d/dr alpha, g_tt, and Gamma^r_tt are all related:
set output 'images/schwarzschild_eos.png'
plot [0:R*10.] schwarzschild_alpha(x) title 'schwarzschild alpha', schwarzschild_dr_alpha(x) title 'schwarzschild d_r alpha', schwarzschild_g_tt(x) title 'g_t_t', -schwarzschild_alpha(x) * schwarzschild_dr_alpha(x) title 'schwarzschild gravity = -Gamma_r_t_t = -alpha d_r alpha'

# set radius and density to earth parameters:
R = 6.37101e+6						# matter radius, in m
mass = 5.9736e+24 * G / c / c		# total mass within radius, in m
density = mass / sphereVolume(R)	# average density, in m^-2
load 'lapse.inc.gnuplot'

set ylabel 'm/s^2'

# show gravity only, in m/s^2
set output 'images/schwarzschild_gravity.png'
set log y
plot [0:R*10.] -newton_gravity(x) title 'newton gravity', -schwarzschild_gravity(x) title 'schwarzschild gravity', -kerr_gravity(x,earth_a) title 'kerr gravity'
unset log y

# show differences in gravity 
set output 'images/gravity_differences.png'
plot [0:R*2.] abs(schwarzschild_gravity(x)) - abs(newton_gravity(x)) title '|schwarzschild|-|newton|', \
				abs(schwarzschild_gravity(x)) - abs(kerr_gravity(x,earth_a)) title '|schwarzschild|-|kerr|', \
				abs(newton_gravity(x)) - abs(kerr_gravity(x,earth_a)) title '|newton|-|kerr|'

# rotation-less Kerr isn't the same as Schwarzschild ... It weaker than rotation Kerr.  It is within 2e-10 of rotation Kerr ... 2e-11 at the Earth's surface
# note that Schwarzschild (which is rotation-less) is stronger than rotation-Kerr and rotation-less Kerr by about 1.4e-8

set output 'images/kerr with vs. without rotation.png'
plot [0:R*2.] abs(kerr_gravity(x,earth_a)) - abs(kerr_gravity(x,0)) title '|kerr w/rotation|-|kerr w/o rotation|'

set output 'images/kerr fast rotation vs. kerr rotation.png'
plot [0:R*2.] abs(kerr_gravity(x,200000*earth_a)) - abs(kerr_gravity(x,earth_a)) title '|kerr fast rotation|-|kerr w/rotation| ... how to double gravity only using rotation'
