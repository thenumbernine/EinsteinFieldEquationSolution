/*
instead of calculating the finite difference expression up front and extracting like coefficients,
how about using an iterative solver and only numerically calculating finite differences?
might get into trouble ensuring the hamiltonian or momentum constraints are fulfilled.
*/
#include "Tensor/Tensor.h"
#include "Tensor/Grid.h"
#include "Tensor/Inverse.h"
#include "Solvers/JFNK.h"
#include "Parallel/Parallel.h"
#include <functional>
#include <chrono>

Parallel::Parallel parallel(8);

void time(const std::string name, std::function<void()> f) {
	std::cerr << name << " ... ";
	std::cerr.flush();
	auto start = std::chrono::high_resolution_clock::now();
	f();
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end - start;
	std::cerr << "(" << diff.count() << "s)" << std::endl;
}

//opposite upper/lower of what is already in Tensor/Inverse.h
namespace Tensor {
template<typename Real_>
struct InverseClass<Tensor<Real_, Symmetric<Upper<3>, Upper<3>>>> {
	typedef Real_ Real;
	typedef Tensor<Real, Upper<3>, Upper<3>> InputType;
	typedef Tensor<Real, Lower<3>, Lower<3>> OutputType;
	OutputType operator()(const InputType &a, const Real &det) const {
		OutputType result;
		//symmetric, so only do Lower triangular
		result(0,0) = det22(a(1,1), a(1,2), a(2,1), a(2,2)) / det;
		result(1,0) = det22(a(1,2), a(1,0), a(2,2), a(2,0)) / det;
		result(1,1) = det22(a(0,0), a(0,2), a(2,0), a(2,2)) / det;
		result(2,0) = det22(a(1,0), a(1,1), a(2,0), a(2,1)) / det;
		result(2,1) = det22(a(0,1), a(0,0), a(2,1), a(2,0)) / det;
		result(2,2) = det22(a(0,0), a(0,1), a(1,0), a(1,1)) / det;
		return result;
	}
};
}

//enable this to use the J vector to calculate the A vector, to calculate the E & B vectors
//disable to specify E & B directly
//#define USE_CHARGE_CURRENT_FOR_EM

using namespace Tensor;

typedef double real;
enum { spatialDim = 3 };
enum { dim = spatialDim+1 };

Tensor<real, Upper<3>> cross(Tensor<real, Upper<3>> a, Tensor<real, Upper<3>> b) {
	Tensor<real, Upper<3>> c;
	c(0) = a(1) * b(2) - a(2) * b(1);
	c(1) = a(2) * b(0) - a(0) * b(2);
	c(2) = a(0) * b(1) - a(1) * b(0);
	return c;
}

//variables used to build the metric 
//10 vars
struct MetricPrims {
	real alpha;
	Tensor<real, Upper<spatialDim>> betaU;
	Tensor<real, Symmetric<Lower<spatialDim>, Lower<spatialDim>>> gammaLL;
	MetricPrims() : alpha(0) {}
};

//variables used to build the stress-energy tensor
struct StressEnergyPrims {
	
	//source terms:
	real rho;	//matter density
	real P;		//pressure ... due to matter.  TODO what about magnetic pressure?
	real eInt;	//specific internal energy
	Tensor<real, Upper<spatialDim>> v;	//3-vel (upper, spatial)

#ifdef USE_CHARGE_CURRENT_FOR_EM
	/*
	one way to reconstruct the E & B fields is by representing the charge and current densities	
	B^i = curl(A^i), E^i = -dA^i/dt - grad(A^0)
	-A^a;u_;u + A^u_;u^;a + R^a_u A^u = 4 pi J^a
	 but A^u_;u = 0 (Lorentz gauge condition)
	for J^0 = charge density and J^i = current density
	so -A^a;u_;u + R^a_u A^u = 4 pi J^a
	therefore div(E^i) = 4 pi J^0 
	and dE^i/dt - curl(B^i) = -4 pi J^i
	...
	can we calculate A (and F and T_EM) from J?
	wouldn't I also need the time derivative (and spatial derivatives, which can be calculated from the spatial slice lattice) 
	in order to compute the E & B fields?
	in fact, can any E & B be computed with a given A?
	can any A be computed by a given E & B?
	*/
	real chargeDensity;
	Tensor<real, Upper<spatialDim>> currentDensity;	//TODO how does this relate to matter density?
#else	//USE_CHARGE_CURRENT_FOR_EM
	//electric & magnetic fields, which are used in the EM stress-energy tensor, and are derived from A, which is the inverse de Rham of J
	//(can any E & B be represented by a valid A, & subsequently J?)
	//these are upper and spatial-only
	Tensor<real, Upper<spatialDim>> E, B;
#endif	//USE_CHARGE_CURRENT_FOR_EM
	/*
	1) specify charge density and current density -- components of J^a
	2) inverse de Rham vector wave operator to solve for A^a via (-A^a;u_;u + R^a_u A^u) / 4 pi = J^a
	3) F_uv = A_v;u - A_u;v
	4) T_uv = 1/(4 pi)(F_u^a F_va - 1/4 g_uv F^ab F_ab)
	*/
};

/*
natural units ...
1 = c m/s = 299792458 m/s
	1 s = c m
	1 s = 299792458 m
1 = G m^3 / (kg s^2) = 6.67384e-11 m^3 / (kg s^2)
    kg = G m^3 / s^2 = G / c^2 m
	kg = 7.4256484500929e-28 m
1 = kB m^2 kg / (K s^2) = 1.3806488e-23 m^2 kg / (K s^2)
	K = kB kg m^2 / s^2 = kB / c^2 kg = kB G / c^4 m
	K = 1.1407124948367e-67 m
joules: J = kg m^2 / s^2
electronvolts: 1 eV = 1.6e-19 J
Gauss: 1 Gauss^2 = g / (cm s^2) = .1 kg / (m s^2)
	1 Gauss^2 = .1 G/c^4 1/m^2
	Gauss = sqrt(.1 G/c^2) 1/m

I'm going to use meters as my units ...

Radius of Earth = 6.37101e+6 m
Mass of Earth = 5.9736e+24 kg
*/
const real c = 299792458;	// m/s 
const real G = 6.67384e-11;	// m^3 / (kg s^2)
//const real kB = 1.3806488e-23;	// m^2 kg / (K s^2)

int main() {

#if 1	//earth
	const real radius = 6.37101e+6;	// m
	const real mass = 5.9736e+24 * G / c / c;	// m
#endif
#if 0	//sun
	const real radius = 6.960e+8;	// m
	const real mass = 1.9891e+30 * G / c / c;	// m
#endif
	const real volume = 4/3*M_PI*radius*radius*radius;	// m^3
	//earth volume: 1.0832120174985e+21 m^3
	const real density = mass / volume;	// 1/m^2
	//earth density: 4.0950296770075e-24 1/m^2
	//const real schwarzschildRadius = 2 * mass;	//8.87157 mm, which is accurate

	//earth magnetic field at surface: .25-.26 gauss
	//const real magneticField = .45 * sqrt(.1 * G) / c;	// 1/m

	const real boundsSizeInRadius = 2;
	//grid coordinate bounds
	Vector<real, spatialDim> 
		xmin(-boundsSizeInRadius*radius, -boundsSizeInRadius*radius, -boundsSizeInRadius*radius),
		xmax(boundsSizeInRadius*radius, boundsSizeInRadius*radius, boundsSizeInRadius*radius);
	
	int sizei = 128;
	Vector<int, spatialDim> sizev(sizei, sizei, sizei);
	int gridVolume = sizev.volume();
	Vector<real, spatialDim> dx = (xmax - xmin) / sizev;

	//grids
	Grid<Vector<real, spatialDim>, spatialDim> xs;
	Grid<StressEnergyPrims, spatialDim> stressEnergyPrimGrid;
	Grid<MetricPrims, spatialDim> metricPrimGrid;
	
	//some helper storage...
	Grid<Tensor<real, Symmetric<Lower<dim>, Lower<dim>>>, spatialDim> gLLs;
	Grid<Tensor<real, Symmetric<Upper<dim>, Upper<dim>>>, spatialDim> gUUs;
	Grid<Tensor<real, Upper<dim>, Symmetric<Lower<dim>, Lower<dim>>>, spatialDim> GammaULLs;

	time("allocating", [&]{
		xs.resize(sizev);
		stressEnergyPrimGrid.resize(sizev);
		metricPrimGrid.resize(sizev);
		gLLs.resize(sizev);
		gUUs.resize(sizev);
		GammaULLs.resize(sizev);
	});

	//specify coordinates
	RangeObj<spatialDim> range = xs.range();
	time("calculating grid", [&]{
		parallel.foreach(range.begin(), range.end(), [&](const Vector<int, spatialDim>& index) {
			Vector<real, spatialDim>& xi = xs(index);
			for (int j = 0; j < spatialDim; ++j) {
				xi(j) = (xmax(j) - xmin(j)) * ((real)index(j) + .5) / (real)sizev(j) + xmin(j);
			}
		});
	});

	std::function<Vector<int,spatialDim>(Vector<int,spatialDim>,int)> prev = [&](Vector<int, spatialDim> v, int i) -> Vector<int, spatialDim>{
		v(i) = std::max<int>(v(i) - 1, 0);
		return v;
	};
	std::function<Vector<int,spatialDim>(Vector<int,spatialDim>,int)> next = [&](Vector<int, spatialDim> v, int i) -> Vector<int, spatialDim>{
		v(i) = std::min<int>(v(i) + 1, sizev(i) - 1);
		return v;
	};

	auto calc_gLLs_and_gUUs = [&](const real* x) {
		//calculate gLL and gUU from metric primitives
		parallel.foreach(range.begin(), range.end(), [&](const Vector<int, spatialDim>& index) {

			const MetricPrims &metricPrims = *((const MetricPrims*)x + Vector<int, spatialDim>::dot(metricPrimGrid.step, index));

			const Tensor<real, Upper<spatialDim>> &betaU = metricPrims.betaU;
			const Tensor<real, Symmetric<Lower<spatialDim>, Lower<spatialDim>>> &gammaLL = metricPrims.gammaLL;

			real alpha = metricPrims.alpha;
			real alphaSq = alpha * alpha;

			Tensor<real, Lower<spatialDim>> betaL;
			for (int i = 0; i < spatialDim; ++i) {
				betaL(i) = 0;
				for (int j = 0; j < spatialDim; ++j) {
					betaL(i) += betaU(j) * gammaLL(i,j);
				}
			}
				
			real betaSq = 0;
			for (int i = 0; i < spatialDim; ++i) {
				betaSq += betaL(i) * betaU(i);
			}

			//compute ADM metrics
			
			Tensor<real, Symmetric<Lower<dim>, Lower<dim>>> &gLL = gLLs(index);
			gLL(0,0) = -alphaSq + betaSq;
			for (int i = 0; i < spatialDim; ++i) {
				gLL(i+1,0) = betaL(i);
				for (int j = 0; j < spatialDim; ++j) {
					gLL(i+1,j+1) = gammaLL(i,j);
				}
			}

			//why doesn't this call work?
			//Tensor<real, Symmetric<Upper<spatialDim>, Upper<spatialDim>>> gammaUU = inverse<Tensor<real, Symmetric<Upper<spatialDim>, Upper<spatialDim>>>>(gammaLL);
			//oh well, here's the body:
			//symmetric, so only do Lower triangular
			Tensor<real, Symmetric<Upper<spatialDim>, Upper<spatialDim>>> gammaUU;
			real det = determinant33<real, Tensor<real, Symmetric<Lower<spatialDim>, Lower<spatialDim>>>>(gammaLL);
			gammaUU(0,0) = det22(gammaLL(1,1), gammaLL(1,2), gammaLL(2,1), gammaLL(2,2)) / det;
			gammaUU(1,0) = det22(gammaLL(1,2), gammaLL(1,0), gammaLL(2,2), gammaLL(2,0)) / det;
			gammaUU(1,1) = det22(gammaLL(0,0), gammaLL(0,2), gammaLL(2,0), gammaLL(2,2)) / det;
			gammaUU(2,0) = det22(gammaLL(1,0), gammaLL(1,1), gammaLL(2,0), gammaLL(2,1)) / det;
			gammaUU(2,1) = det22(gammaLL(0,1), gammaLL(0,0), gammaLL(2,1), gammaLL(2,0)) / det;
			gammaUU(2,2) = det22(gammaLL(0,0), gammaLL(0,1), gammaLL(1,0), gammaLL(1,1)) / det;

			Tensor<real, Symmetric<Upper<dim>, Upper<dim>>> &gUU = gUUs(index);
			gUU(0,0) = -1/alphaSq;
			for (int i = 0; i < spatialDim; ++i) {
				gUU(i+1,0) = betaU(i) / alphaSq;
				for (int j = 0; j < spatialDim; ++j) {
					gUU(i+1,j+1) = gammaUU(i,j) - betaU(i) * betaU(j) / alphaSq;
				}
			}
		});
	};

	//depends on calc_gLLs_and_gUUs(x)
	auto calc_GammaULLs = [&](const real* x) {
		parallel.foreach(range.begin(), range.end(), [&](const Vector<int, spatialDim>& index) {
			//derivatives of the metric in spatial coordinates using finite difference
			Tensor<real, Symmetric<Lower<dim>, Lower<dim>>, Lower<dim>> dgLLL;
			for (int a = 0; a < dim; ++a) {
				for (int b = 0; b < dim; ++b) {	
					dgLLL(a,b,0) = 0;	//TODO time derivative? for now assume steady state of metric.  this isn't necessary though ... spacetime vortex?
					for (int i = 0; i < spatialDim; ++i) {
						dgLLL(a,b,i+1) = (gLLs(next(index,i))(a,b) - gLLs(prev(index,i))(a,b)) / (2 * dx(i));
					}
				}
			}
			
			//connections
			Tensor<real, Lower<dim>, Symmetric<Lower<dim>, Lower<dim>>> GammaLLL;
			for (int a = 0; a < dim; ++a) {
				for (int b = 0; b < dim; ++b) {
					for (int c = 0; c <= b; ++c) {
						GammaLLL(a,b,c) = .5 * (dgLLL(a,b,c) + dgLLL(a,c,b) - dgLLL(b,c,a));
					}
				}
			}
			
			const Tensor<real, Symmetric<Upper<dim>, Upper<dim>>> &gUU = gUUs(index);
			
			Tensor<real, Upper<dim>, Symmetric<Lower<dim>, Lower<dim>>> &GammaULL = GammaULLs(index);
			for (int a = 0; a < dim; ++a) {
				for (int b = 0; b < dim; ++b) {
					for (int c = 0; c <= b; ++c) {
						real sum = 0;
						for (int d = 0; d < dim; ++d) {
							sum += gUU(a,d) * GammaLLL(d,b,c);
						}
						GammaULL(a,b,c) = sum;
					}
				}
			}
		});
	};

	//specify stress-energy primitives
	//the stress-energy primitives combined with the current metric are used to compute the stress-energy tensor 
	time("calculating stress-energy primitives", [&]{
		parallel.foreach(range.begin(), range.end(), [&](const Vector<int, spatialDim>& index) {
			StressEnergyPrims &stressEnergyPrims = stressEnergyPrimGrid(index);
			real r = Vector<real, spatialDim>::length(xs(index));
			stressEnergyPrims.rho = r < radius ? density : 0;	// average density of Earth in m^-2
			stressEnergyPrims.eInt = 0;	//internal energy / temperature of the Earth?
			stressEnergyPrims.P = 0;	//pressure inside the Earth?
			//stressEnergyPrims.E = Tensor<real, Upper<spatialDim>>(0,0,0);	//electric field
			//stressEnergyPrims.B = Tensor<real, Upper<spatialDim>>(0,0,0);	//magnetic field
		});
	});

	//initialize metric primitives
	time("calculating metric primitives", [&]{
		parallel.foreach(range.begin(), range.end(), [&](const Vector<int, spatialDim>& index) {
			MetricPrims& metricPrims = metricPrimGrid(index);
			real r = Vector<real, spatialDim>::length(xs(index));
			Vector<real, spatialDim> xi = xs(index);
			//mass within radius
			real m = 4/3 * M_PI * pow(std::min<real>(r, radius), 3) * density;
			//substitute the schwarzschild R for 2 m(r)
			real R = 2 * m;
#if 0	//completely flat space
			metricPrims.alpha = 1;
			for (int i = 0; i < spatialDim; ++i) {
				metricPrims.betaU(i) = 0;
				for (int j = 0; j <= i; ++j) {
					metricPrims.gammaLL(i,j) = (real)(i == j);
				}
			}
#endif
#if 1	//schwarzschild metric
			/*
			g_ti = beta_i = 0
			g_tt = -alpha^2 + beta^2 = -alpha^2 = -1 + Rs/r <=> alpha = sqrt(1 - Rs/r)
			g_ij = gamma_ij = delta_ij + x^i x^j / r^2 R/(r - R)		<- but x is upper, and you can't lower it without specifying gamma_ij
			 ... which might be why the contravariant spatial metrics of spherical and cartesian look so similar 
			*/
			metricPrims.alpha = sqrt(1 - R/r); 
			for (int i = 0; i < spatialDim; ++i) {
				metricPrims.betaU(i) = 0;
				for (int j = 0; j <= i; ++j) {
					metricPrims.gammaLL(i,j) = (real)(i == j) + xi(i)/r * xi(j)/r * R/(r - R);
					/*
					dr^2's coefficient
					spherical: 1/(1 - R/r) = 1/((r - R)/r) = r/(r - R)
					spherical contravariant: 1 - R/r
					cartesian contravariant: delta_ij - x/r y/r R/r
					hmm, contravariant terms of cartesian vs spherical look more similar than covariant terms do
				
					in the OV metric, dr^2's coefficient is exp(2 Lambda) = 1/(1 - 2 m(r) / r) where m(r) is the enclosing mass
					so the contravariant coefficient would be exp(-2 Lambda) = 1 - 2 m(r) / r
					I'm going to do the lazy thing and guess this converts to delta^ij - 2 m(r) x^i x^j / r^3
					*/
				}
			}
#endif
		});
	});

	time("calculating g_ab and g^ab", [&]{
		calc_gLLs_and_gUUs((real*)metricPrimGrid.v);
	});

	time("calculating Gamma^a_bc", [&]{
		calc_GammaULLs((real*)metricPrimGrid.v);
	});

	time("calculating gravitational force", [&]{
		//this is printing output, so don't do it in parallel		
		for (RangeObj<spatialDim>::iterator iter = range.begin(); iter != range.end(); ++iter) {
			MetricPrims& metricPrims = metricPrimGrid(iter.index);
			real r = Vector<real, spatialDim>::length(xs(iter.index));
			real det = determinant33<real, Tensor<real, Symmetric<Lower<spatialDim>, Lower<spatialDim>>>>(metricPrims.gammaLL);
		
			
			//numerical computational...		
			Tensor<real, Upper<dim>, Symmetric<Lower<dim>, Lower<dim>>> &GammaULL = GammaULLs(iter.index);
			//why is there a 2 here?
			//Gamma^r_tt = Gamma^i_tt dr/dx^i
			//r^2 = x^2 + y^2 + z^2
			//so dr/dx^i = x^i / r
			Vector<real, spatialDim> xi = xs(iter.index);
			real gravity = -(GammaULL(1,0,0) * xi(0) / r
							+ GammaULL(2,0,0) * xi(1) / r
							+ GammaULL(3,0,0) * xi(2) / r)
							* c * c;	//times c twice because of the two timelike components of a^i = Gamma^i_tt
		
			//substitute the schwarzschild R for 2 m(r)
			real matterRadius = std::min<real>(r, radius);
			real volumeOfMatterRadius = 4/3*M_PI*matterRadius*matterRadius*matterRadius;
			real matterEnclosed = density * volumeOfMatterRadius;	// m^3
			real R = 2 * matterEnclosed;
			real analyticalGravity = R * (r - R) / (2 * r * r * r) * c * c;

			std::cout << iter.index 
				<< " rho=" << stressEnergyPrimGrid(iter.index).rho 
				<< " r=" << r 
				<< " det=" << det 
				<< " gravity=" << gravity
				<< " analyticalGravity=" << analyticalGravity
				<< std::endl;
		}
	});

	std::cerr << "done!" << std::endl;
exit(0);

	assert(sizeof(MetricPrims) == sizeof(real) * 10);	//this should be 10 real numbers and nothing else
	Solvers::JFNK<real> jfnk(
		sizeof(MetricPrims) / sizeof(real) * gridVolume,	//n = vector size
		(real*)metricPrimGrid.v,	//x = state of vector
		[&](real* y, const real* x) {	//A = vector function to minimize
			real* ystart = y;

			calc_gLLs_and_gUUs(x);
			calc_GammaULLs(x);

			parallel.foreach(range.begin(), range.end(), [&](const Vector<int, spatialDim>& index) {

				//connection derivative
				Tensor<real, Upper<dim>, Symmetric<Lower<dim>, Lower<dim>>, Lower<dim>> dGammaULLL;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						for (int c = 0; c <= b; ++c) {
							dGammaULLL(a,b,c,0) = 0;	//TODO explicit calculate Gamma^a_bc,t in terms of alpha, beta^i, gamma_ij
							for (int i = 0; i < spatialDim; ++i) {
								dGammaULLL(a,b,c,i+1) = (GammaULLs(next(index,i))(a,b,c) - GammaULLs(prev(index,i))(a,b,c)) / (2 * dx(i));
							}
						}
					}
				}
				
				const Tensor<real, Upper<dim>, Symmetric<Lower<dim>, Lower<dim>>> &GammaULL = GammaULLs(index);

				Tensor<real, Upper<dim>, Lower<dim>, Lower<dim>, Lower<dim>> GammaSqULLL;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						for (int c = 0; c < dim; ++c) {
							for (int d = 0; d < dim; ++d) {
								real sum = 0;
								for (int e = 0; e < dim; ++e) {
									sum += GammaULL(a,e,d) * GammaULL(e,b,c);
								}
								GammaSqULLL(a,b,c,d) = sum;
							}
						}
					}
				}

				Tensor<real, Upper<dim>, Lower<dim>, Lower<dim>, Lower<dim>> RiemannULLL;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						for (int c = 0; c < dim; ++c) {
							for (int d = 0; d < dim; ++d) {
								RiemannULLL(a,b,c,d) = dGammaULLL(a,b,d,c) - dGammaULLL(a,b,c,d) + GammaSqULLL(a,b,d,c) - GammaSqULLL(a,b,c,d);
							}
						}
					}
				}

				Tensor<real, Symmetric<Lower<dim>, Lower<dim>>> RicciLL;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						real sum = 0;
						for (int c = 0; c < dim; ++c) {
							sum += RiemannULLL(c,a,c,b);
						}
						RicciLL(a,b) = sum;
					}
				}
				
				const Tensor<real, Symmetric<Upper<dim>, Upper<dim>>> &gUU = gUUs(index);
				
				real Gaussian = 0;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						Gaussian += gUU(a,b) * RicciLL(a,b);
					}
				}
				
				const Tensor<real, Symmetric<Lower<dim>, Lower<dim>>> &gLL = gLLs(index);

				Tensor<real, Symmetric<Lower<dim>, Lower<dim>>> EinsteinLL;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						EinsteinLL(a,b) = RicciLL(a,b) - .5 * Gaussian * gLL(a,b);
					}
				}

				//now we want to find the zeroes of EinsteinLL(a,b) - 8 pi T(a,b)
				// ... which is 10 zeroes ...
				// ... and we are minimizing the inputs to our metric ...
				// alpha, beta x3, gamma x6
				// ... which is 10 variables
				// tada!			

				const MetricPrims &metricPrims = *((const MetricPrims*)x + Vector<int, spatialDim>::dot(metricPrimGrid.step, index));
				real alpha = metricPrims.alpha;
				real alphaSq = alpha * alpha;
				const Tensor<real, Upper<spatialDim>> &betaU = metricPrims.betaU;
				const Tensor<real, Symmetric<Lower<spatialDim>, Lower<spatialDim>>> &gammaLL = metricPrims.gammaLL;

				//now compute stress-energy based on source terms
				//notice: stress energy depends on gLL (i.e. alpha, betaU, gammaLL), which it is solving for, so this has to be recalculated every iteration
		
				StressEnergyPrims &stressEnergyPrims = stressEnergyPrimGrid(index);

				//electromagnetic stress-energy

#ifdef USE_CHARGE_CURRENT_FOR_EM
				Tensor<real, Upper<dim>> JU;
				JU(0) = stressEnergyPrims.chargeDensity;
				for (int i = 0; i < spatialDim; ++i) {
					JU(i+1) = stressEnergyPrims.currentDensity(i);
				}
				Tensor<real, Upper<dim>> AU = JU;
				/*
				A^a;u = A^a_;v g^uv = (A^a_,v + Gamma^a_wv A^w) g^uv
				A^a;u_;u = A^a;u_,u + Gamma^a_bu A^b;u + Gamma^u_bu A^a;b
						= (A^a_;v g^uv = (A^a_,v + Gamma^a_wv A^w) g^uv)_,u
							+ Gamma^a_bu (A^b_,v + Gamma^b_wv A^w) g^uv
							- Gamma^u_bu (A^a_,v + Gamma^a_wv A^w) g^bv
				((A^a_,b + Gamma^a_cb A^c) + R^a_b A^b) / (4 pi) = J^a
				*/
				JFNK(dim,
					JU.v,
					[&](double* y, const double* x) {
						for (int a = 0; i < dim; ++a) {
							
						}
					}
				);
#else	//USE_CHARGE_CURRENT_FOR_EM
				Tensor<real, Upper<spatialDim>> E = stressEnergyPrims.E;
				Tensor<real, Upper<spatialDim>> B = stressEnergyPrims.B;
#endif

				//electromagnetic stress-energy
				real ESq = 0, BSq = 0;
				for (int i = 0; i < spatialDim; ++i) {
					for (int j = 0; j < spatialDim; ++j) {
						ESq += E(i) * E(j) * gammaLL(i,j);
						BSq += B(i) * B(j) * gammaLL(i,j);
					}
				}
				Tensor<real, Upper<spatialDim>> S = cross(E, B);
				
				Tensor<real, Symmetric<Lower<dim>, Lower<dim>>> T_EM_UU;
				T_EM_UU(0,0) = (ESq + BSq) / alphaSq / (8 * M_PI);
				for (int i = 0; i < spatialDim; ++i) {
					T_EM_UU(i+1,0) = (-betaU(i) * (ESq + BSq) / alphaSq + 2 * S(i) / alpha) / (8 * M_PI);
					for (int j = 0; j <= i; ++j) {
						T_EM_UU(i+1,j+1) = -2 * (E(i) * E(j) + B(i) * B(j) + (S(i) * B(j) + S(j) * B(i)) / alpha) + betaU(i) * betaU(j) * (ESq + BSq) / alphaSq;
						if (i == j) {
							T_EM_UU(i+1,j+1) += ESq + BSq;
						}
						T_EM_UU(i+1,j+1) /= 8 * M_PI;
					}
				}

				Tensor<real, Upper<dim>, Lower<dim>> T_EM_LU;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						real sum = 0;
						for (int w = 0; w < dim; ++w) {
							sum += gLL(a,w) * T_EM_UU(w,b);
						}
						T_EM_LU(a,b) = sum;
					}
				}

				Tensor<real, Symmetric<Lower<dim>, Lower<dim>>> T_EM_LL;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b <= a; ++b) {
						real sum = 0;
						for (int w = 0; w < dim; ++w) {
							sum += T_EM_LU(a,w) * gLL(w,b);
						}
						T_EM_LL(a,b) = sum;
					}
				}

				//matter stress-energy

				Tensor<real, Upper<spatialDim>> &v = stressEnergyPrims.v;

				//Lorentz factor
				real vLenSq = 0;
				for (int i = 0; i < spatialDim; ++i) {
					for (int j = 0; j < spatialDim; ++j) {
						vLenSq += v(i) * v(j) * gammaLL(i,j);
					}
				}
				real W = 1 / sqrt( 1 - sqrt(vLenSq) );

				//4-vel upper
				Tensor<real, Upper<dim>> uU;
				uU(0) = W;
				for (int i = 0; i < spatialDim; ++i) {
					uU(i+1) = W * v(i);
				}

				//4-vel lower
				Tensor<real, Lower<dim>> uL;
				for (int a = 0; a < dim; ++a) {
					uL(a) = 0;
					for (int b = 0; b < dim; ++b) {
						uL(a) += uU(b) * gLL(b,a);
					}
				}

				/*
				Right now I'm using the SRHD T_matter_ab = (rho + rho eInt) u_a u_b + P P_ab
					for P^ab = g^ab + u^a u^b = projection tensor
				TODO viscious matter stress-energy: MTW 22.16d: T^ab = rho u^a u^b + (P - zeta theta) P^ab - 2 eta sigma^ab + q^a u^b + u^a q^b
				T_heat_ab = q^a u^b + u^a q^b 
					q^a = the heat-flux 4-vector
				T_viscous_ab = -2 eta sigma^ab - zeta theta P^ab 
					eta >= 0 = coefficient of dynamic viscosity
					zeta >= 0 = coefficient of bulk viscosity
					sigma^ab = 1/2(u^a_;u P^ub + u^b_;u P^ua) - theta P^ab / 3 = shear
					theta = u^a_;a = expansion
				*/	
				Tensor<real, Symmetric<Lower<dim>, Lower<dim>>> T_matter_LL;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b <= a; ++b) {
						T_matter_LL(a,b) = uL(a) * uL(b) * (stressEnergyPrims.rho * (1 + stressEnergyPrims.eInt) + stressEnergyPrims.P) + gLL(a,b) * stressEnergyPrims.P;
					}
				}
				
				//total stress-energy	
				Tensor<real, Symmetric<Lower<dim>, Lower<dim>>> T_LL = T_EM_LL + T_matter_LL;
			
				/*
				now solve the linear system G_uv = G(g_uv) = 8 pi T_uv for g_uv 
				i.e. A(x) = b, assuming A is linear ...
				but it looks like, because T is based on g, it will really look like G(g_uv) = 8 pi T(g_uv, source terms)
				*/

				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b <= a; ++b) {
						*y = EinsteinLL(a,b) - 8 * M_PI * T_LL(a,b);
						++y;
					}
				}
			});

			int n = y - ystart;
			if (n != gridVolume * 10) {
				throw Common::Exception() << "expected " << (gridVolume * 10) << " but found " << n << " entries";
			}

		},
		1e-7, //newton stop epsilon
		100, //newton maxi ter
		1e-7, //gmres stop epsilon
		gridVolume * 10, //gmres max iter
		gridVolume	//gmres restart iter
	);

	jfnk.stopCallback = [&]()->bool{
		printf("jfnk iter %d alpha %f residual %.16f\n", jfnk.iter, jfnk.alpha, jfnk.residual);
		return false;
	};
	jfnk.gmres.stopCallback = [&]()->bool{
		printf("gmres iter %d residual %.16f\n", jfnk.gmres.iter, jfnk.gmres.residual);
		return false;
	};
	jfnk.solve();
}
