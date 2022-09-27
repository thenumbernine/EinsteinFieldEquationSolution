/*
instead of calculating the finite difference expression up front and extracting like coefficients,
how about using an iterative solver and only numerically calculating finite differences?
might get into trouble ensuring the hamiltonian or momentum constraints are fulfilled.
*/

#include <fstream>
#include <iostream>
#include <cmath>

namespace std {
bool isfinite(__float128 f) {
	double d = f;
	return std::isfinite(d);
}
ostream& operator<<(ostream& o, __float128 f) {
	double d = f;
	return o << d;
}
}


#include "Solver/ConjGrad.h"
#include "Solver/ConjRes.h"
#include "Solver/GMRES.h"
#include "Solver/JFNK.h"
#include "Parallel/Parallel.h"
#include "LuaCxx/State.h"
#include "LuaCxx/Ref.h"
#include "Tensor/Inverse.h"
#include "Tensor/Derivative.h"
#include "Tensor/Grid.h"
#include "Tensor/Tensor.h"
#include "Tensor/vec.h"
#include "Common/Exception.h"
#include "Common/Macros.h"
#include <functional>
#include <chrono>
#include <iomanip>

#define CONVERGE_ALPHA_ONLY
//#define PRINTTIME
#define PRINT_RANGES

int const numThreads = 8;
Parallel::Parallel parallel(numThreads);

int const partialDerivativeOrder = 2;//options are 2, 4, 6, 8

void time(std::string const & name, std::function<void()> f) {
	std::cout << name << " ... ";
	std::cout.flush();
	auto start = std::chrono::high_resolution_clock::now();
	f();
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end - start;
	std::cout << "(" << diff.count() << "s)" << std::endl;
}

//enable this to use the J vector to calculate the A vector, to calculate the E & B vectors
//disable to specify E & B directly
//#define USE_CHARGE_CURRENT_FOR_EM

//using real = float;
using real = double;
//using real = __float128;

static constexpr auto subDim = 3;	//spatial dim
static constexpr auto dim = subDim+1;

template<typename T>
using vecDim = Tensor::v2::_vec<T,dim>;

template<typename T>
using vecSubDim = Tensor::v2::_vec<T,subDim>;

template<typename T>
using symDim = Tensor::v2::_sym<T,dim>;

template<typename T>
using symSubDim = Tensor::v2::_sym<T,subDim>;


//subDim
using real_Sub = vecSubDim<real>;
using real_Sub_Sub = vecSubDim<vecSubDim<real>>;
using real_SymSub = symSubDim<real>;

//dim
using real_Dim = vecDim<real>;

using real_SymDim = symDim<real>;
using real_Dim_Dim = vecDim<vecDim<real>>;

using real_Dim_SymDim = vecDim<symDim<real>>;
using real_SymDim_Dim = symDim<vecDim<real>>;

using real_Dim_Dim_Dim = vecDim<vecDim<vecDim<real>>>;	//Tensor::v2::_tensor<real, dim, dim, dim>;
/* 
TODO generalize: 
real_SymDim_Dim = _tensor<real, _realSym<dim>, _real<dim>>;
real_Dim_Dim_Dim = _tensor<real, _real<dim>, _real<dim>, _real<dim>>;
where each arg is a single-arg template that gets concatenated
*/

using real_Dim_SymDim_Dim = vecDim<symDim<vecDim<real>>>;
using real_SymDim_SymDim = symDim<symDim<real>>;
using real_Dim_Dim_Dim_Dim = vecDim<vecDim<vecDim<vecDim<real>>>>; //Tensor::v2::_tensor<real, dim, dim, dim, dim>;

//mixed subDim & dim
using real_Sub_SymDim = vecSubDim<symDim<real>>;
using real_Sub_Dim_SymDim = vecSubDim<vecDim<symDim<real>>>;
using real_Sub_SymDim_Dim = vecSubDim<symDim<vecDim<real>>>;

using int_Sub = vecSubDim<int>;

template<typename T>
struct gridFromPtr {
	Tensor::Grid<T, subDim> grid;
	gridFromPtr(real* ptr, int_Sub size) {
		grid.v = (T*)ptr;
		grid.size = size;
		grid.step(0) = 1;
		for (int i = 1; i < subDim; ++i) {
			grid.step(i) = grid.step(i-1) * grid.size(i-1);
		}
	}
	gridFromPtr(real const * ptr, int_Sub size) {
		grid.v = (T*)ptr;
		grid.size = size;
		grid.step(0) = 1;
		for (int i = 1; i < subDim; ++i) {
			grid.step(i) = grid.step(i-1) * grid.size(i-1);
		}
	}
	~gridFromPtr() {
		grid.v = {};	//don't deallocate
	}
};

real determinant44(real_SymDim const & m) {
	real sign = 1;
	real sum = 0;
	for (int a = 0; a < 4; ++a) {
		real_Sub_Sub subm;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				subm(i,j) = m(i+1, j + (j>=a));
			}
		}
		sum += sign * m(0,a) * Tensor::v2::determinant(subm);
		sign = -sign;
	}
	return sum;
}

//variables used to build the metric 
//dim * (dim+1) / 2 vars
struct MetricPrims {
	real alphaMinusOne;
	real_Sub betaU;
	real_SymSub hLL;	//h_ij = gamma_ij - delta_ij
	MetricPrims() : alphaMinusOne(0) {}
};

//variables used to build the stress-energy tensor
struct StressEnergyPrims {
	
	//source terms:
	real rho;	//matter density
	real P;		//pressure ... due to matter.  TODO what about magnetic pressure?
	real eInt;	//specific internal energy

	bool useV;
	real_Sub v;	//3-vel (upper, spatial)

	bool useEM;
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
	real_Sub currentDensity;	//TODO how does this relate to matter density?
#else	//USE_CHARGE_CURRENT_FOR_EM
	//electric & magnetic fields, which are used in the EM stress-energy tensor, and are derived from A, which is the inverse de Rham of J
	//(can any E & B be represented by a valid A, & subsequently J?)
	//these are upper and spatial-only
	real_Sub E, B;
#endif	//USE_CHARGE_CURRENT_FOR_EM
	/*
	1) specify charge density and current density -- components of J^a
	2) inverse de Rham vector wave operator to solve for A^a via (-A^a;u_;u + R^a_u A^u) / 4 pi = J^a
	3) F_uv = A_v;u - A_u;v
	4) T_uv = 1/(4 pi)(F_u^a F_va - 1/4 g_uv F^ab F_ab)
	*/

	StressEnergyPrims()
	: rho(0), P(0), eInt(0)	//technically these should all be positive...
	, useV(false)
	, useEM(false)
#ifdef USE_CHARGE_CURRENT_FOR_EM
	,chargeDensity(0)
#endif
	{}
};

/*
natural units ...
1 = c m/s = 299792458 m/s
	1 s = c m
	1 s = 299792458 m
1 = G m^3 / (kg s^2) = 6.67384e-11 m^3 / (kg s^2)
    kg = G m^3 / s^2 = G / c^2 m
	kg = 7.4256484500929e-28 m
1 = ke kg m^3 / (s^2 C^2)
	C = sqrt(ke kg m^3 / s^2) = sqrt(ke G / c^2 m^4 / s^2) = sqrt(ke G m^2 / c^4) = sqrt(ke G) / c^2 m
	C = 8.6172202304995e-18 m
1 = kB m^2 kg / (K s^2) = 1.3806488e-23 m^2 kg / (K s^2)
	K = kB kg m^2 / s^2 = kB / c^2 kg = kB G / c^4 m
	K = 1.1407124948367e-67 m

I'm going to use meters as my units ...

Radius of Earth = 6.37101e+6 m
Mass of Earth = 5.9736e+24 kg
*/
real const c = 299792458;	// m/s 
real const G = 6.67384e-11;	// m^3 / (kg s^2)
real const ke = 8987551787.3682;	//kg m^3 / (s^2 C^2)
//real const kB = 1.3806488e-23;	// m^2 kg / (K s^2)

//grid coordinate bounds
Tensor::Vector<real, subDim> xmin, xmax;
Tensor::Vector<int, subDim> sizev;
int gridVolume;
Tensor::Vector<real, subDim> dx;

//some helper storage...
Tensor::Grid<real_SymDim, subDim> gLLs;
Tensor::Grid<real_SymDim, subDim> gUUs;
Tensor::Grid<real_SymDim, subDim> dt_gLLs;
Tensor::Grid<real_SymDim, subDim> d2t_gLLs;
//Tensor::Grid<real_SymDim, subDim> dt_gUUs;
Tensor::Grid<real_SymDim_Dim, subDim> dgLLLs;
//Tensor::Grid<real_Dim_SymDim, subDim> GammaLLLs;
Tensor::Grid<real_Dim_SymDim, subDim> GammaULLs;

template<typename CellType>
void allocateGrid(Tensor::Grid<CellType, subDim>& grid, std::string name, Tensor::Vector<int, subDim> sizev, size_t& totalSize) {
	size_t size = sizeof(CellType) * sizev.volume();
	totalSize += size;
	std::cout << name << ": " << size << " bytes, running total: " << totalSize << std::endl;
	grid.resize(sizev);
}

static real_SymSub make_delta3LL() {
	real_SymSub delta3LL;
	delta3LL(0,0) = 1.;
	delta3LL(1,1) = 1.;
	delta3LL(2,2) = 1.;
	return delta3LL;
}
real_SymSub delta3LL = make_delta3LL();

/*
calculates contents of gUUs, gLLs
	(incl. first deriv: dt_gLLs)
	(incl. second deriv: dt_gUUs)
x is an array of MetricPrims[gridVolume]
*/
void calc_gLLs_and_gUUs(
	//input
	Tensor::Grid<MetricPrims, subDim> const & metricPrimGrid,
	Tensor::Grid<MetricPrims, subDim> const & dt_metricPrimGrid,	//first deriv
	//output:
	Tensor::Grid<real_SymDim, subDim>& gLLs,
	Tensor::Grid<real_SymDim, subDim>& gUUs,
	Tensor::Grid<real_SymDim, subDim>& dt_gLLs	//first deriv
) {
	//calculate gLL and gUU from metric primitives
	Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
	parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
		MetricPrims const &metricPrims = metricPrimGrid(index);
		real alpha = metricPrims.alphaMinusOne + 1.;
//debugging
//looks like, for the Krylov solvers, we have a problem of A(x) producing zero and A(A(x)) giving us zeros here ... which cause singular basises
//lesson: the problem isn't linear.  don't use Krylov solvers.
assert(alpha != 0);
		real_Sub const & betaU = metricPrims.betaU;
		real_SymSub gammaLL = metricPrims.hLL + delta3LL;
	
		//I can only solve for one of these.  or can I do more?  without solving for d/dt variables, I am solving 10 unknowns for 10 constraints. 
		MetricPrims const & dt_metricPrims = dt_metricPrimGrid(index);

		real alphaSq = alpha * alpha;

		real_Sub betaL;
		for (int i = 0; i < subDim; ++i) {
			betaL(i) = 0;
			for (int j = 0; j < subDim; ++j) {
				betaL(i) += betaU(j) * gammaLL(i,j);
			}
		}
			
		real betaSq = 0;
		for (int i = 0; i < subDim; ++i) {
			betaSq += betaL(i) * betaU(i);
		}

		//compute ADM metrics
	
		//g_ab
		real_SymDim & gLL = gLLs(index);
		gLL(0,0) = -alphaSq + betaSq;
		for (int i = 0; i < subDim; ++i) {
			gLL(i+1,0) = betaL(i);
			for (int j = 0; j < subDim; ++j) {
				gLL(i+1,j+1) = gammaLL(i,j);
			}
		}

		real dt_alpha = dt_metricPrims.alphaMinusOne;
		real_Sub const & dt_betaU = dt_metricPrims.betaU;
		real_SymSub dt_gammaLL = dt_metricPrims.hLL;
		
		//g_ab,t
		real_SymDim & dt_gLL = dt_gLLs(index);
		//g_tt,t = (-alpha^2 + beta^2),t
		//		 = -2 alpha alpha,t + 2 beta^i_,t beta_i + beta^i beta^j gamma_ij,t
		dt_gLL(0,0) = -2. * alpha * dt_alpha;
		for (int i = 0; i < subDim; ++i) {
			dt_gLL(0,0) += 2. * dt_betaU(i) * betaL(i);
			for (int j = 0; j < subDim; ++j) {
				dt_gLL(0,0) += betaU(i) * betaU(j) * dt_gammaLL(i,j);
			}
		}
		//g_ti = beta_i,t = (beta^j gamma_ij),j 
		//		 = beta^j_,t gamma_ij + beta^j gamma_ij,t
		for (int i = 0; i < subDim; ++i) {
			dt_gLL(i+1,0) = 0;
			for (int j = 0; j < subDim; ++j) {
				dt_gLL(i+1,0) += dt_betaU(j) * gammaLL(i,j) + betaU(j) * dt_gammaLL(i,j);
			}
		}
		//g_ij,t = gamma_ij,t
		for (int i = 0; i < subDim; ++i) {
			for (int j = 0; j <= i; ++j) {
				dt_gLL(i+1,j+1) = dt_gammaLL(i,j);
			}
		}
	
		//gamma^ij
		real_SymSub gammaUU = Tensor::v2::inverse(gammaLL);

		//g^ab
		real_SymDim & gUU = gUUs(index);
		gUU(0,0) = -1/alphaSq;
		for (int i = 0; i < subDim; ++i) {
			gUU(i+1,0) = betaU(i) / alphaSq;
			for (int j = 0; j <= i; ++j) {
				gUU(i+1,j+1) = gammaUU(i,j) - betaU(i) * betaU(j) / alphaSq;
			}
		}
//debugging
#ifdef DEBUG
for (int a = 0; a < dim; ++a) {
	for (int b = 0; b <= a; ++b) {
		assert(gUU(a,b) == gUU(a,b));
	}
}
#endif
		//gamma^ij_,t
		//https://math.stackexchange.com/questions/1187861/derivative-of-transpose-of-inverse-of-matrix-with-respect-to-matrix
		//d/dt AInv_kl = dAInv_kl / dA_ij d/dt A_ij
		//= -AInv_ki (d/dt A_ij) AInv_jl
		real_SymSub dt_gammaUU;
		real_Sub_Sub tmp;
		for (int k = 0; k < subDim; ++k) {
			for (int j = 0; j < subDim; ++j) {
				real sum = 0;
				for (int i = 0; i < subDim; ++i) {
					sum -= gammaUU(k,i) * dt_gammaLL(i,j);
				}
				tmp(k,j) = sum;
			}
		}
		for (int k = 0; k < subDim; ++k) {
			for (int l = 0; l <= k; ++l) {	//dt_gammaUU is symmetric
				real sum = 0;
				for (int j = 0; j < subDim; ++j) {
					sum += tmp(k,j) * gammaUU(j,l);
				}
				dt_gammaUU(k,l) = sum;
			}
		}

		/*
		//g^ab_,t
		real_SymDim &dt_gUU = dt_gUUs(index);
		//g^tt_,t = (-1/alpha^2),t = 2 alpha,t / alpha^3
		dt_gUU(0,0) = 2. * dt_alpha / (alpha * alphaSq);
		//g^ti_,t = (beta^i/alpha^2),t = beta^i_,t / alpha^2 - 2 beta^i alpha,t / alpha^3
		for (int i = 0; i < subDim; ++i) {
			dt_gUU(i,0) = (dt_betaU(i) * alpha - 2. * betaU(i) * dt_alpha) / (alpha * alphaSq);
			for (int j = 0; j <= i; ++j) {
				//g^ij_,t = (gamma^ij - beta^i beta^j / alpha^2),t = gamma^ij_,t - beta^i_,t beta^j / alpha^2 - beta^i beta^j_,t / alpha^2 + 2 beta^i beta^j alpha_,t / alpha^3
				dt_gUU(i,j) = dt_gammaUU(i,j) - (dt_betaU(i) * betaU(j) + betaU(i) * dt_betaU(j)) / alphaSq + 2. * betaU(i) * betaU(j) * dt_alpha / (alpha * alphaSq);
			}
		}
		*/
	});
}

/*
calculates contents of GammaULLs
	incl second derivs: GammaLLLs
depends on gLLs, gUUs
	incl first derivs: dt_gLLs
prereq: calc_gLLs_and_gUUs()
*/
void calc_GammaULLs(
	//input:
	Tensor::Grid<real_SymDim, subDim> const & gLLs,
	Tensor::Grid<real_SymDim, subDim> const & gUUs,
	Tensor::Grid<real_SymDim, subDim> const & dt_gLLs,	//first deriv
	//output:
	Tensor::Grid<real_Dim_SymDim, subDim>& GammaULLs
	//Tensor::Grid<real_Dim_SymDim, subDim>& GammaLLLs	//second deriv
) {
	Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
	parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
		//derivatives of the metric in spatial coordinates using finite difference
		//the templated method (1) stores derivative first and (2) only stores spatial
#warning "TODO get this working with ::v2"
#if 1 
		real_Sub_SymDim dgLLL3;
#else		
		real_Sub_SymDim dgLLL3 = Tensor::partialDerivative<partialDerivativeOrder, real, subDim, real_SymDim>(
			index, dx,
			[&](Tensor::Vector<int, subDim> index)
				-> real_SymDim
			{
				for (int i = 0; i < subDim; ++i) {
					index(i) = std::max<int>(0, std::min<int>(sizev(i)-1, index(i)));
				}
				return gLLs(index);
			}
		);
#endif		
		real_SymDim const & dt_gLL = dt_gLLs(index);
		real_SymDim_Dim & dgLLL = dgLLLs(index);
		for (int a = 0; a < dim; ++a) {
			for (int b = 0; b < dim; ++b) {	
				dgLLL(a,b,0) = dt_gLL(a,b);
				for (int i = 0; i < subDim; ++i) {
					dgLLL(a,b,i+1) = dgLLL3(i,a,b);
				}
			}
		}
		
		//connections
		//real_Dim_SymDim & GammaLLL = GammaLLLs(index);
		real_Dim_SymDim GammaLLL;
		for (int a = 0; a < dim; ++a) {
			for (int b = 0; b < dim; ++b) {
				for (int c = 0; c <= b; ++c) {
					GammaLLL(a,b,c) = .5 * (dgLLL(a,b,c) + dgLLL(a,c,b) - dgLLL(b,c,a));
//debugging
assert(GammaLLL(a,b,c) == GammaLLL(a,b,c));
				}
			}
		}
		
		real_SymDim const & gUU = gUUs(index);		
		real_Dim_SymDim & GammaULL = GammaULLs(index);
		for (int a = 0; a < dim; ++a) {
			for (int b = 0; b < dim; ++b) {
				for (int c = 0; c <= b; ++c) {
					real sum = 0;
					for (int d = 0; d < dim; ++d) {
						sum += gUU(a,d) * GammaLLL(d,b,c);
					}
					GammaULL(a,b,c) = sum;
//debugging
assert(GammaULL(a,b,c) == GammaULL(a,b,c));
				}
			}
		}
	});
}

/*
index is the location in the grid
depends on GammaULLs, gLLs, gUUs
	second deriv: dt_gUUs, GammaLLLs
prereq: calc_gLLs_and_gUUs(), calc_GammaULLs()
*/
real_SymDim calc_EinsteinLL(
	//input
	Tensor::Vector<int, subDim> index,
	Tensor::Grid<real_SymDim, subDim> const & gLLs,
	Tensor::Grid<real_SymDim, subDim> const & gUUs,
	Tensor::Grid<real_Dim_SymDim, subDim> const & GammaULLs
) {
	real_SymDim const & gUU = gUUs(index);
	real_SymDim_Dim const & dgLLL = dgLLLs(index);
	real_Dim_SymDim const & GammaULL = GammaULLs(index);
#if 0	//calc first derivative of Gamma^a_bc's
	//connection derivative
	real_Sub_Dim_SymDim dGammaLULL3 = Tensor::partialDerivative<partialDerivativeOrder, real, subDim, real_Dim_SymDim>(
		index, dx, [&](Tensor::Vector<int, subDim> index) -> real_Dim_SymDim {
			for (int i = 0; i < subDim; ++i) {
				index(i) = std::max<int>(0, std::min<int>(sizev(i)-1, index(i)));
			}

//debugging
#ifdef DEBUG
for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
		for (int k = 0; k <= j; ++k) {
			assert(GammaULL(i,j,k) == GammaULL(i,j,k));
		}
	}
}
#endif
			return GammaULLs(index);
		});
	
	//real_SymDim const & dt_gUU = dt_gUUs(index);
	//real_Dim_SymDim const & GammaLLL = GammaLLLs(index);

	real_Dim_SymDim_Dim dGammaULLL;
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b < dim; ++b) {
			for (int c = 0; c <= b; ++c) {
				//TODO
				//Gamma^a_bc,t = (g^ad Gamma_dbc),t = g^ad_,t Gamma_dbc + g^ad Gamma_dbc,t = g^ad_,t Gamma_dbc + 1/2 g^ad Gamma_dbc,t
				//but this is where the 2nd derivative comes in, and that means providing 2 sets of initial condition metric primitives
				real sum = 0;
				//for (int d = 0; d < dim; ++d) {
				//	sum += dt_gUU(a,d) * GammaLLL(d,b,c) + gUU(a,d) * dt_GammaLLL(d,b,c);
				//}
				dGammaULLL(a,b,c,0) = sum;
				//finite difference
				for (int i = 0; i < subDim; ++i) {
					dGammaULLL(a,b,c,i+1) = dGammaLULL3(i,a,b,c);
//debugging				
assert(dGammaULLL(a,b,c,i+1) == dGammaULLL(a,b,c,i+1));
				}
			}
		}
	}
#else	//calc 2nd deriv of g_ab's and 1st deriv of g^ab's instead
/*
Gamma^a_bc,d = 1/2 (g^ae (g_eb,c + g_ec,b - g_bc,e)),d
	= g^ae_,d Gamma_ebc + 1/2 g^ae (g_eb,cd + g_ec,bd - g_bc,ed)
	= -g^ae g_ef,d Gamma^f_bc + 1/2 g^ae (g_eb,cd + g_ec,bd - g_bc,ed)
*/
	//g^ae g_ef,d
	real_Dim_Dim_Dim gdgULL;
	for (int a = 0; a < 4; ++a) {
		for (int f = 0; f < 4; ++f) {
			for (int d = 0; d < 4; ++d) {
				real sum = 0;
				for (int e = 0; e < 4; ++e) {
					sum += gUU(a,e) * dgLLL(e,f,d);
				}
				gdgULL(a,f,d) = sum;
			}
		}
	}

	//g_ab,ci
#warning "TODO get this working with ::v2"
#if 1
	real_Sub_SymDim_Dim d2gLLLL3;
#else	
	real_Sub_SymDim_Dim d2gLLLL3 = Tensor::partialDerivative<partialDerivativeOrder, real, subDim, real_SymDim_Dim>(
		index, dx,
		[&](Tensor::Vector<int, subDim> index)
			-> real_SymDim_Dim
		{
			for (int i = 0; i < subDim; ++i) {
				index(i) = std::max<int>(0, std::min<int>(sizev(i)-1, index(i)));
			}
			return dgLLLs(index);
		}
	);
#endif

	//g_ab,cd
	real_SymDim_SymDim d2gLLLL;
	real_SymDim const & d2t_gLL = d2t_gLLs(index);
	for (int a = 0; a < 4; ++a) {
		for (int b = 0; b <= a; ++b) {
			for (int c = 0; c < 4; ++c) {
				for (int d = 0; d <= c; ++d) {
					if (c == d) {
						//then do a 2nd deriv
						if (c == 0) {
							//c or d is 0 = time?  then we need dt2_gLL ...
							d2gLLLL(a,b,c,d) = d2t_gLL(a,b);
						} else {
							Tensor::Vector<int, subDim> ixp = index;
							ixp(c) = std::min(ixp(c) + 1, sizev(c)-1);
							Tensor::Vector<int, subDim> ixm = index;
							ixm(c) = std::max(ixm(c) - 1, 0);
							
							real_SymDim_Dim const & dgLLL_ix = dgLLLs(index);
							real_SymDim_Dim const & dgLLL_ixp = dgLLLs(ixp);
							real_SymDim_Dim const & dgLLL_ixm = dgLLLs(ixm);
							
							d2gLLLL(a,b,c,d) = 
								(dgLLL_ixp(a,b,c) 
								- dgLLL_ix(a,b,c) * 2.
								+ dgLLL_ixm(a,b,c)) / (dx(c-1) * dx(c-1));
						}
					} else {
						//then do a 1st deriv
						d2gLLLL(a,b,c,d) = d2gLLLL3(d,a,b,c-1);
					}
				}
			}
		}
	}

	//Gamma^a_bcd = -g^ae g_ef,d Gamma^f_bc + 1/2 g^ae (g_eb,cd + g_ec,bd - g_bc,ed)
	real_Dim_SymDim_Dim dGammaULLL;
	for (int a = 0; a < 4; ++a) {
		for (int b = 0; b < 4; ++b) {
			for (int c = 0; c < 4; ++c) {
				for (int d = 0; d < 4; ++d) {
					real sum = 0;
					for (int f = 0; f < 4; ++f) {
						//-g^ae g_ef,d Gamma^f_bc
						sum -= gdgULL(a,f,d) * GammaULL(f,b,c);
						//+ 1/2 g^af (g_fb,cd + g_fc,bd - g_bc,fd)
						sum += .5 * gUU(a,f) * (d2gLLLL(f,b,c,d) + d2gLLLL(f,c,b,d) - d2gLLLL(b,c,f,d));
					}
					dGammaULLL(a,b,c,d) = sum;
				}
			}
		}
	}
#endif
#if 0	//calculate the Riemann, then the Ricci
	real_Dim_Dim_Dim_Dim GammaSqULLL;
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b < dim; ++b) {
			for (int c = 0; c < dim; ++c) {
				for (int d = 0; d < dim; ++d) {
					real sum = 0;
					for (int e = 0; e < dim; ++e) {
						sum += GammaULL(a,e,d) * GammaULL(e,b,c);
					}
					GammaSqULLL(a,b,c,d) = sum;
//debugging
assert(GammaSqULLL(a,b,c,d) == GammaSqULLL(a,b,c,d));				
				}
			}
		}
	}

	/*
	technically the last two are antisymmetric ... but I don't think I have that working yet ... 
	...and if I stored RiemannLLLL then the first two would be antisymmetric as well, and the pairs of the first two and second two would be symmetric

	Symmetry goes away with torsion

	Gamma^a_bc e_a = e_b;c

	Torsion = T^a_bc = Gamma^a_bc - Gamma^a_cb
	...so order of Gamma matters
	...so what is the covariant derivative?
	is it v^a_;b = v^a_,b + Gamma^a_cb v^c ?	<- I'm using this so the derivative index is last
	or v^a_;b = v^a_,b + Gamma^a_bc v^c ?
	
	v^a_;bc - v^a_;cb = Riemann^a_bcd v^d
	(v^a_,b + Gamma^a_db v^d);c - (v^a_,c + Gamma^a_dc v^d);b
	v^a_,bc + 

	*/
	real_Dim_Dim_Dim_Dim RiemannULLL;
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b < dim; ++b) {
			for (int c = 0; c < dim; ++c) {
				for (int d = 0; d < dim; ++d) {
					RiemannULLL(a,b,c,d) = dGammaULLL(a,b,d,c) - dGammaULLL(a,b,c,d) + GammaSqULLL(a,b,d,c) - GammaSqULLL(a,b,c,d);
//debugging
assert(RiemannULLL(a,b,c,d) == RiemannULLL(a,b,c,d));
				}
			}
		}
	}

	real_SymDim RicciLL;
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b < dim; ++b) {
			real sum = 0;
			for (int c = 0; c < dim; ++c) {
				sum += RiemannULLL(c,a,c,b);
			}
			RicciLL(a,b) = sum;
//debugggin
assert(RicciLL(a,b) == RicciLL(a,b));
		}
	}
#else	//just calculate the Ricci
	real_Dim Gamma12L;
	for (int a = 0; a < dim; ++a) {
		real sum = 0;
		for (int b = 0; b < dim; ++b) {
			sum += GammaULL(b,b,a);
		}
		Gamma12L(a) = sum;
	}
	
	//R_ab = Gamma^c_ab,c - Gamma^c_ac,b + Gamma^d_ab Gamma^c_cd - Gamma^d_ac Gamma^c_bd
	real_SymDim RicciLL;
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b < dim; ++b) {
			real sum = 0;
			for (int c = 0; c < dim; ++c) {
				sum += dGammaULLL(c,a,b,c) - dGammaULLL(c,a,c,b) + GammaULL(c,a,b) * Gamma12L(c);
				for (int d = 0; d < dim; ++d) {
					sum -= GammaULL(d,a,c) * GammaULL(c,b,d);
				}
			}
			RicciLL(a,b) = sum;
//debugggin
assert(RicciLL(a,b) == RicciLL(a,b));
		}
	}
#endif
	
	real Gaussian = 0;
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b < dim; ++b) {
			Gaussian += gUU(a,b) * RicciLL(a,b);
		}
	}
//debugging
assert(Gaussian == Gaussian);

	real_SymDim const & gLL = gLLs(index);
	real_SymDim EinsteinLL;
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b < dim; ++b) {
			EinsteinLL(a,b) = RicciLL(a,b) - .5 * Gaussian * gLL(a,b);
//debugging
assert(EinsteinLL(a,b) == EinsteinLL(a,b));
		}
	}

	return EinsteinLL;
}

/*
calls calc_EinsteinLL at each point
stores G_ab 
*/
void calc_EinsteinLLs(
	//input
	Tensor::Grid<real_SymDim, subDim> const & gLLs,
	Tensor::Grid<real_SymDim, subDim> const & gUUs,
	Tensor::Grid<real_Dim_SymDim, subDim> const & GammaULLs,
	//output:
	Tensor::Grid<real_SymDim, subDim> & EinsteinLLs
) {
	Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
	assert(sizeof(real_SymDim) == sizeof(MetricPrims));	//10 reals for both
	parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
		real_SymDim & EinsteinLL = EinsteinLLs(index);
		EinsteinLL = calc_EinsteinLL(index, gLLs, gUUs, GammaULLs);

//debugging
#ifdef DEBUG
for (int a = 0; a < dim; ++a) {
	for (int b = 0; b <= a; ++b) {
		assert(EinsteinLL(a,b) == EinsteinLL(a,b));
	}
}
#endif
	});
}

/*
index is the location in the grid
metricPrims are the metric primitives at that point
returns 8*pi*T_ab for T_ab the stress energy at that point
depends on: stressEnergyPrims, gLLs, calc_gLLs_and_gUUs()
now compute stress-energy based on source terms
notice: stress energy depends on gLL (i.e. alpha, betaU, gammaLL), which it is solving for, so this has to be recalculated every iteration
*/
real_SymDim calc_8piTLL(
	MetricPrims const & metricPrims,
	real_SymDim const & gLL,
	real_SymDim const & gUU,
	StressEnergyPrims const & stressEnergyPrims
) {
	real_SymSub gammaLL = metricPrims.hLL + delta3LL;

	//electromagnetic stress-energy
	real_SymDim T_EM_LL;
	if (stressEnergyPrims.useEM) {

#ifdef USE_CHARGE_CURRENT_FOR_EM
		real_Dim JU;
		JU(0) = stressEnergyPrims.chargeDensity;
		for (int i = 0; i < subDim; ++i) {
			JU(i+1) = stressEnergyPrims.currentDensity(i);
		}
		real_Dim AU = JU;
		/*
		A^a;u = A^a_;v g^uv = (A^a_,v + Gamma^a_wv A^w) g^uv
		A^a;u_;u = A^a;u_,u + Gamma^a_bu A^b;u + Gamma^u_bu A^a;b
				= (A^a_;v g^uv = (A^a_,v + Gamma^a_wv A^w) g^uv)_,u
					+ Gamma^a_bu (A^b_,v + Gamma^b_wv A^w) g^uv
					- Gamma^u_bu (A^a_,v + Gamma^a_wv A^w) g^bv
		((A^a_,b + Gamma^a_cb A^c) + R^a_b A^b) / (4 pi) = J^a
		
		but then we have to run the JFNK to find A across all the grid, so we can use spatial derivatives
		...and what about time derivatives?  assume steady state, or provide the time derivative 
		*/
		JFNK(dim,
			dim, 
			JU.v,
			[&](double * y, double const * x) {
				for (int a = 0; i < dim; ++a) {
					y[0] = 
				}
			}
		);
#else	//USE_CHARGE_CURRENT_FOR_EM
		real_Sub E = stressEnergyPrims.E;
		real_Sub B = stressEnergyPrims.B;
#endif
		
		//should I be doing a full 4x4 determinant?
		//if converging beta then yep
		real sqrtDetG = sqrt(fabs(determinant44(gLL)));
	
		//n_a = t_,a
		real_Dim nL;
		for (int a = 0; a < 4; ++a) {
			nL(a) = a == 0 ? 1 : 0;
		}
	
		//n^a = g^ab n_b
		real_Dim nU;
		for (int a = 0; a < 4; ++a) {
			real sum = 0;
			for (int b = 0; b < 4; ++b) {
				sum += gUU(a,b) * nL(b);
			}
			nU(a) = sum;
		}

		real_Dim BU, EU;
		for (int i = 0; i < 3; ++i) {
			BU(i+1) = B(i);
			EU(i+1) = E(i);
		}
	
		real_Dim EL;
		for (int a = 0; a < 4; ++a) {
			real sum = 0;
			for (int b = 0; b < 4; ++b) {
				sum += gLL(a,b) * EU(b);
			}
			EL(a) = sum;
		}

		//Faraday
		//F_ab = n^d epsilon_adbf g^fg B_g - n_a E_b + E_a n_b
		//	= E_a n_b - n_a E_b - sqrt|g| n^c [abcd] B^d
		//	= E_a n_b - n_a E_b - sqrt|g| n^c [abcd] B^d
		//assuming E and B are specified in Cartesian coordinates as one-forms
		real_Dim_Dim F_LL;
		for (int a = 0; a < dim; ++a) {
			for (int b = 0; b < dim; ++b) {
				F_LL(a,b) = EL(a) * nL(b) - nL(a) * EL(b);
			}
		}

		real tmp;
#define ADD_B_TO_F(a,b,c,d)	tmp = sqrtDetG * (BU(c) * nU(d) - BU(d) * nU(c)); F_LL(a,b) += tmp; F_LL(b,a) -= tmp;
		ADD_B_TO_F(0,1,2,3)	//+ 0 1 2 3, - 0 1 3 2
		ADD_B_TO_F(0,2,3,1)	//+ 0 2 3 1, - 0 2 1 3
		ADD_B_TO_F(0,3,1,2)	//+ 0 3 1 2, - 0 3 2 1	
		ADD_B_TO_F(1,2,0,3)	//+ 1 2 0 3, - 1 2 3 0
		ADD_B_TO_F(1,3,2,0)	//+ 1 3 2 0, - 1 3 0 2
		ADD_B_TO_F(2,3,0,1)	//+ 2 3 0 1, - 2 3 1 0

		real_Dim_Dim F_LU;
		for (int a = 0; a < dim; ++a) {
			for (int b = 0; b < dim; ++b) {
				real sum = 0;
				for (int c = 0; c < dim; ++c) {
					sum += F_LL(a,c) * gUU(c,b);
				}
				F_LU(a,b) = sum;
			}
		}
		
		real_Dim_Dim F_UU;
		for (int a = 0; a < dim; ++a) {
			for (int b = 0; b < dim; ++b) {
				real sum = 0;
				for (int c = 0; c < dim; ++c) {
					sum += gUU(a,c) * F_LU(c,b);
				}
				F_UU(a,b) = sum;
			}
		}

		real FNormSq = 0;
		for (int a = 0; a < dim; ++a) {
			for (int b = 0; b < dim; ++b) {
				FNormSq += F_LL(a,b) * F_UU(a,b);
			}
		}
		
		//T_ab = 1/(4 pi) (F_ac F_b^c - 1/4 g_ab F_cd F^cd)
		for (int a = 0; a < dim; ++a) {
			for (int b = 0; b <= a; ++b) {
				real sum = -gLL(a,b) * FNormSq / 4.;
				for (int c = 0; c < dim; ++c) {
					sum += F_LL(a,c) * F_LU(b,c);
				}
				T_EM_LL(a,b) = sum / (4. * M_PI);
			}
		}
	}

	//matter stress-energy

	real_Dim uL;
	if (stressEnergyPrims.useV) {
		real_Sub const & v = stressEnergyPrims.v;

		//Lorentz factor
		real vLenSq = 0;
		for (int i = 0; i < subDim; ++i) {
			for (int j = 0; j < subDim; ++j) {
				vLenSq += v(i) * v(j) * gammaLL(i,j);
			}
		}
		real W = 1 / sqrt( 1 - sqrt(vLenSq) );

		//4-vel upper
		real_Dim uU;
		uU(0) = W;
		for (int i = 0; i < subDim; ++i) {
			uU(i+1) = W * v(i);
		}

		//4-vel lower
		for (int a = 0; a < dim; ++a) {
			uL(a) = 0;
			for (int b = 0; b < dim; ++b) {
				uL(a) += uU(b) * gLL(b,a);
			}
		}
	} else {
		for (int a = 0; a < dim; ++a) {
			uL(a) = gLL(a,0);
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
	real_SymDim T_matter_LL;
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b <= a; ++b) {
			T_matter_LL(a,b) = uL(a) * uL(b) * (stressEnergyPrims.rho * (1 + stressEnergyPrims.eInt) + stressEnergyPrims.P) + gLL(a,b) * stressEnergyPrims.P;
		}
	}

	//total stress-energy	
	real_SymDim _8piT_LL;
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b <= a; ++b) {
			_8piT_LL(a,b) = (T_EM_LL(a,b) + T_matter_LL(a,b)) * 8. * M_PI;
		}
	}
	return _8piT_LL;
}

/*
x holds a grid of MetricPrims
stores at y a grid of the values (G_ab - 8 pi T_ab)
depends on: calc_gLLs_and_gUUs(), calc_GammaULLs()
*/
void calc_EFE_constraint(
	//input
	Tensor::Grid<MetricPrims, subDim> const & metricPrimGrid,
	Tensor::Grid<StressEnergyPrims, subDim> const & stressEnergyPrimGrid,
	//output
	Tensor::Grid<real_SymDim, subDim>& EFEGrid
) {
	Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
	parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
		
		//for the JFNK solver that doesn't cache the EinsteinLL tensors
		// no need to allocate for both an EinsteinLL grid and a EFEGrid
		real_SymDim EinsteinLL = calc_EinsteinLL(index, gLLs, gUUs, GammaULLs);
		
		//now we want to find the zeroes of EinsteinLL(a,b) - 8 pi T(a,b)
		// ... which is 10 zeroes ...
		// ... and we are minimizing the inputs to our metric ...
		// alpha, beta x3, gamma x6
		// ... which is 10 variables
		// tada!
		real_SymDim _8piT_LL = calc_8piTLL(
			metricPrimGrid(index),
			gLLs(index),
			gUUs(index),
			stressEnergyPrimGrid(index));
		
		/*
		now solve the linear system G_uv = G(g_uv) = 8 pi T_uv for g_uv 
		i.e. A(x) = b, assuming A is linear ...
		but it looks like, because T is based on g, it will really look like G(g_uv) = 8 pi T(g_uv, source terms)
		*/
		
		real_SymDim &EFE = EFEGrid(index);
		for (int a = 0; a < dim; ++a) {
			for (int b = 0; b <= a; ++b) {
				EFE(a,b) = EinsteinLL(a,b) - _8piT_LL(a,b);
			}
		}
	});
}

struct EFESolver {
	int maxiter;
	EFESolver(int maxiter_) : maxiter(maxiter_) {}
	size_t getN() { return sizeof(MetricPrims) / sizeof(real) * gridVolume; }
	virtual void solve(
		//input/output
		Tensor::Grid<MetricPrims, subDim>& metricPrimGrid,
		//input
		Tensor::Grid<MetricPrims, subDim> const & dt_metricPrimGrid,	//first deriv
		Tensor::Grid<StressEnergyPrims, subDim> const & stressEnergyPrimGrid
	) = 0;
};

/*
use a linear solver and treat G_ab = 8 pi T_ab like a linear system A x = b for x = (alpha, beta, gamma), A x = G_ab(x), and b = 8 pi T_ab ... which is also a function of x ...
nothing appears to be moving ...
or it's diverging ...
looks like there's an inherent problem in all the Krylov solvers, because they're based on A^n(x), and beacuse initial-condition flat A(x) gives all zeroes for G_ab(x)
... and as long as 'x' is the primitive variables, the second that x=0 for A(x) we end up with a singular basis, and everything fails.
... so, for constant A(x) = G_ab(x), G_ab(G_ab(x)) is all nans 
*/
struct KrylovSolver : public EFESolver {
	using Super = EFESolver;
	using Super::Super;
	
	Tensor::Grid<real_SymDim, subDim> _8piTLLs;
	std::shared_ptr<Solver::Krylov<real>> krylov;

	KrylovSolver(int maxiter)
	: Super(maxiter)
	, _8piTLLs(sizev)
	{}

	virtual char const * name() = 0;

	/*
	x holds a grid of MetricPrims 
	calls calc_8piTLL() at each point on the grid
	stores results in _8piTLLs
	depends on: calc_gLLs_and_gUUs()
	*/
	void calc_8piTLLs(
		//input
		Tensor::Grid<MetricPrims, subDim> const & metricPrimGrid,
		Tensor::Grid<StressEnergyPrims, subDim> const & stressEnergyPrimGrid,
		//output
		Tensor::Grid<real_SymDim, subDim>& _8piTLLs
	) {
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			_8piTLLs(index) = calc_8piTLL(
				metricPrimGrid(index), 
				gLLs(index), 
				gUUs(index),
				stressEnergyPrimGrid(index));
		});
	}

	void linearFunc(
		real * y,
		real const * x,
		//extra inputs
		Tensor::Grid<MetricPrims, subDim> const & dt_metricPrimGrid	//first deriv
	) {
//debugging
#ifdef DEBUG
for (int i = 0; i < (int)getN(); ++i) {
	assert(x[i] == x[i]);
}
#endif

#ifdef PRINTTIME
		std::cout << "iteration " << jfnk.iter << std::endl;
		time("calculating g_ab and g^ab", [&](){
#endif				
			Tensor::Grid<MetricPrims, subDim> metricPrimGrid(sizev, (MetricPrims*)x);
			calc_gLLs_and_gUUs(
				//input
				metricPrimGrid,
				dt_metricPrimGrid,	//first deriv
				//output
				gLLs, gUUs, dt_gLLs);
#ifdef PRINTTIME
		});
		time("calculating Gamma^a_bc", [&](){
#endif				
			calc_GammaULLs(gLLs, gUUs, dt_gLLs, GammaULLs);
#ifdef PRINTTIME
		});
		time("calculating G_ab", [&]{
#endif
			Tensor::Grid<real_SymDim, subDim> EinsteinLLs(sizev, (real_SymDim*)y);
			calc_EinsteinLLs(
				//input
				gLLs, gUUs, GammaULLs,
				//output
				EinsteinLLs);
//debugging
#ifdef DEBUG
for (int i = 0; i < (int)getN(); ++i) {
	assert(y[i] == y[i]);
}
#endif

#ifdef PRINTTIME
		});
#endif				

		//here's me abusing GMRES.
		//I'm updating the 'b' vector mid-algorithm since it is dependent on the 'x' vector
		//maybe I shouldn't be doing so here, but instead only before every linear solver solve()?
		//that way the 'b' vector is constant during the linear solver solve() ...
#if 0
#ifdef PRINTTIME
		time("calculating T_ab", [&]{
#endif				
			calc_8piTLLs(Tensor::Grid<MetricPrims const, subDim>(sizev, (MetricPrims const *)x), stressEnergyPrimGrid, _8piTLLs);
#ifdef PRINTTIME
		});
#endif
#endif
	}

	virtual std::shared_ptr<Solver::Krylov<real>> makeSolver(
		Tensor::Grid<MetricPrims, subDim> & metricPrimGrid,
		Tensor::Grid<MetricPrims, subDim> const & dt_metricPrimGrid	//first deriv
	) = 0;
	
	virtual void solve(
		//input/output
		Tensor::Grid<MetricPrims, subDim>& metricPrimGrid,
		//input
		Tensor::Grid<MetricPrims, subDim> const & dt_metricPrimGrid,	//first deriv
		Tensor::Grid<StressEnergyPrims, subDim> const & stressEnergyPrimGrid
	) {
		time("calculating T_ab", [&]{
			calc_8piTLLs(metricPrimGrid, stressEnergyPrimGrid, _8piTLLs);
		});
		
		krylov = makeSolver(
			metricPrimGrid, 
			dt_metricPrimGrid	//first deriv
		);
		
		//seems this is stopping too early, so try scaling both x and y ... or at least the normal that is used ...
#if 0
		krylov->MInv = [&](real * y, real const * x) {
			for (int i = 0; i < (int)getN(); ++i) {
				y[i] = x[i] / (8. * M_PI) * c * c / G / 1000.;
			}
		};
#endif		
		krylov->stopCallback = [&]()->bool{
			std::cout << name() << " iter " << krylov->getIter() << " residual " << krylov->getResidual() << std::endl;
			return false;
		};
		time("solving", [&](){
			krylov->solve();
		});
	}
};

struct ConjGradSolver : public KrylovSolver {
	using Super = KrylovSolver;
	using Super::Super;
	
	virtual char const * name() { return "conjgrad"; }	
	
	virtual std::shared_ptr<Solver::Krylov<real>> makeSolver(
		Tensor::Grid<MetricPrims, subDim> & metricPrimGrid,
		Tensor::Grid<MetricPrims, subDim> const & dt_metricPrimGrid	//first deriv
	) {
		return std::make_shared<Solver::ConjGrad<real>>(
			getN(),
			(real *)metricPrimGrid.v,
			(real const *)_8piTLLs.v,
			[&](real* y, real const * x) { linearFunc(y, x, dt_metricPrimGrid); },
			1e-100,	//epsilon
			getN()	//maxiter
		);
	}
};

struct ConjRes : public Solver::ConjRes<real> {
	using Solver::ConjRes<real>::ConjRes;
	virtual real calcResidual(real rNormL2, real bNormL2, real const * r) {
//debugging
#ifdef DEBUG
for (int i = 0; i < (int)n; ++i) {
	assert(r[i] == r[i]);
}
#endif
		std::cout << "ConjRes::calcResidual"
			<< " n=" << n
			<< " iter=" << iter
			<< " rNormL2=" << rNormL2
			<< " bNormL2=" << bNormL2
			<< std::endl;
		
		return rNormL2;
	}
};

struct ConjResSolver : public KrylovSolver {
	using Super = KrylovSolver;
	using Super::Super;
	
	virtual char const * name() { return "conjres"; }

	virtual std::shared_ptr<Solver::Krylov<real>> makeSolver(
		Tensor::Grid<MetricPrims, subDim> & metricPrimGrid,
		Tensor::Grid<MetricPrims, subDim> const & dt_metricPrimGrid	//first deriv
	) {
		return std::make_shared<ConjRes>(
			getN(),
			(real *)metricPrimGrid.v,
			(real const *)_8piTLLs.v,
			[&](real* y, real const * x) { linearFunc(y, x, dt_metricPrimGrid); },
			1e-100,	//epsilon
			getN()	//maxiter
		);
	}
};

struct GMRES : public Solver::GMRES<real> {
	using Solver::GMRES<real>::GMRES;
	virtual real calcResidual(real rNormL2, real bNormL2, real const * r) {
#if 0
		//error is 16, is sqrt(total sum of errors) which is 256, which is 4 * 64
		// 64 is the # of grid elements, 4 is how much error per grid
		// because the inputs have 4 1's and the rest is 0's.
		real real_rNormL2 = Solver::Vector<real>::normL2(n, r);
		std::cout << "GMRES::calcResidual"
			<< " iter=" << iter
			<< " rNormL2=" << rNormL2
			<< " bNormL2=" << bNormL2
			<< " real_rNormL2=" << real_rNormL2 
			<< std::endl;
		int e = 0;
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		std::for_each(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			std::cout << "r[" << index << "] =";
			for (int j = 0; j < 10; ++j, ++e) {
				std::cout << " " << r[e];
			}
			std::cout << std::endl;
		});
		return real_rNormL2;
#endif
		return rNormL2;
	}
};

struct GMRESSolver : public KrylovSolver {
	using Super = KrylovSolver;
	using Super::Super;
	
	virtual char const * name() { return "gmres"; }

	virtual std::shared_ptr<Solver::Krylov<real>> makeSolver(
		Tensor::Grid<MetricPrims, subDim> & metricPrimGrid,
		Tensor::Grid<MetricPrims, subDim> const & dt_metricPrimGrid	//first deriv
	) {
		return std::make_shared<GMRES>(
			getN(),	//n = vector size
			(real *)metricPrimGrid.v,		//x = state vector
			(real const *)_8piTLLs.v,		//b = solution vector
			[&](real* y, real const * x) { linearFunc(y, x, dt_metricPrimGrid); },	//A = linear function to solve x for A(x) = b
			1e-100,			//epsilon
			getN(),			//maxiter ... = n^3 ... 262144
			100				//restart
		);
		//the largest allocation is restart * maxiter, which would be 100 * n^3, 
		// for a 64^3 grid is 2,621,440,000
		// for a 32^3 grid is 327,680,000
	}
};

//scale input to keep solution from being zero 
double const jfnkInputScale = 1;
//double const jfnkInputScale = 8 * M_PI * c * c / G / 1000.;	//m -> g/cm^3

//scale output to keep residual from being zero 
double const jfnkOutputScale = 1; 
//double const jfnkOutputScale = 8 * M_PI * c * c / G / 1000.;	//m -> g/cm^3

struct JFNK : public Solver::JFNK<real> {
	using Super = typename Solver::JFNK<real>;
	using Super::JFNK;
	virtual real calcResidual(real const * r, real alpha) const {
		real residual = Super::calcResidual(r, alpha);
		residual *= jfnkOutputScale;
		
		std::cout << "JFNK::calcResidual"
			<< " n=" << n
			<< " iter=" << iter
			<< " alpha=" << std::setprecision(49) << alpha << std::setprecision(6)
			<< " residual=" << std::setprecision(49) << residual  << std::setprecision(6)
			<< std::endl;
		return residual;
	}
	size_t getN() const { return n; }
};

//use JFNK
//as soon as this passes 'restart' it explodes.
struct JFNKSolver : public EFESolver {
	using Super = EFESolver;

	Tensor::Grid<real_SymDim, subDim> EFEGrid;	

	JFNKSolver(int maxiter)
	: Super(maxiter)
	, EFEGrid(sizev)
	{
	}

	virtual void solve(
		//input/output
		Tensor::Grid<MetricPrims, subDim> & metricPrimGrid,
		//input
		Tensor::Grid<MetricPrims, subDim> const & dt_metricPrimGrid,	//first deriv
		Tensor::Grid<StressEnergyPrims, subDim> const & stressEnergyPrimGrid
	) {
	
		std::ofstream jfnkFile("jfnk.txt");
		jfnkFile << "#iter residual alpha" << std::endl;
	
		std::ofstream gmresFile("gmres.txt");
		gmresFile << "#jfnk_iter gmres_iter residual" << std::endl;
		
		assert(sizeof(MetricPrims) == sizeof(EFEGrid.v[0]));	//this should be 10 real numbers and nothing else
		
#ifdef CONVERGE_ALPHA_ONLY
		std::vector<real> alphaMinusOnes(gridVolume);
		for (int i = 0; i < gridVolume; ++i) {
			alphaMinusOnes[i] = metricPrimGrid.v[i].alphaMinusOne;
			//scale up alphas before feeding them to the EFE constraint, so they are further from zero
			alphaMinusOnes[i] *= jfnkInputScale;
		}
#endif
		
		int const gmresRestart = 100;
		JFNK jfnk(
#ifdef CONVERGE_ALPHA_ONLY
			gridVolume,	//n = vector size
			alphaMinusOnes.data(),	//x = state vector
#else
			getN(),	//n = vector size
			(real*)metricPrimGrid.v,	//x = state vector
#endif			
			[&](real* y, real const * x) {	//A = vector function to minimize

#ifdef CONVERGE_ALPHA_ONLY
				//Tensor::Grid<MetricPrims, subDim> metricPrimGrid(sizev);
				for (int k = 0; k < gridVolume; ++k) {
					metricPrimGrid.v[k].alphaMinusOne = x[k];
					//scale alphaMinusOne's back down now that we're inside the linear function 
					metricPrimGrid.v[k].alphaMinusOne /= jfnkInputScale;
				}
#else
				Tensor::Grid<MetricPrims, subDim> metricPrimGrid(sizev, (MetricPrims*)x); 
#endif

#ifdef PRINTTIME
				std::cout << "iteration " << jfnk.iter << std::endl;
				time("calculating g_ab and g^ab", [&](){
#endif			
				//g_ab = [-1/alpha^2, beta^i/alpha, gamma_ij]
				//g^ab = inv(g_ab)
				calc_gLLs_and_gUUs(
					//input:
					metricPrimGrid,
					dt_metricPrimGrid,	//first deriv
					//output:
					gLLs, gUUs, dt_gLLs);
#ifdef PRINTTIME
				});
#endif

#ifdef PRINTTIME
				time("calculating Gamma^a_bc", [&](){
#endif			
				//Gamma^a_bc = 1/2 g^ad (g_db,c + g_dc,b - g_bc,d)
				calc_GammaULLs(gLLs, gUUs, dt_gLLs, GammaULLs);
#ifdef PRINTTIME
				});
#endif

#ifdef CONVERGE_ALPHA_ONLY
				Tensor::Grid<real_SymDim, subDim> EFEGrid(sizev);
#else				
				Tensor::Grid<real_SymDim, subDim> EFEGrid(sizev, (real_SymDim*)y);
#endif			
				//EFE_ab = G_ab - 8 pi T_ab
				//T_ab = stress energy constraint, whose calculations depend on g_ab and the stress-energy primitives 
				//G_ab = R_ab - 1/2 R g_ab
				//R = g^ab R_ab
				//R_ab = R^c_acb = (pick a more optimized implementation)
				//R^c_acb = Gamma^c_ab,c - Gamma^c_ac,b + Gamma^c_dc Gamma^d_ab - Gamma^c_db Gamma^d_ac

#ifdef PRINTTIME
				time("calculating G_ab = 8 pi T_ab", [&]{
#endif
				calc_EFE_constraint(
					metricPrimGrid,
					stressEnergyPrimGrid,
					EFEGrid);
#ifdef PRINTTIME
				});
#endif

//scale up the EFE constraint here, so the residual gets a better value
#if 1
				for (int k = 0; k < gridVolume; ++k) {
					for (int a = 0; a < dim; ++a) {
						for (int b = 0; b <= a; ++b) {
							EFEGrid.v[k](a,b) *= jfnkOutputScale;
						}
					}
				}
#endif

#ifdef CONVERGE_ALPHA_ONLY
				for (int k = 0; k < gridVolume; ++k) {
					real sum = 0;
					for (int a = 0; a < dim; ++a) {
						for (int b = 0; b <= a; ++b) {
							real d = EFEGrid.v[k](a,b);
							sum += d * d;
						}
					}
					y[k] = sum;
				}
#endif

#if 0 //debug output
std::cout << "efe constraint" << std::endl;
for (int i = 0; i < gridVolume*10; ++i) {
	std::cout << " " << y[i];
}
std::cout << std::endl;
//this is happening after the first gmres iteration, but I don't think it should ...
std::cout << "got jfnk residual == 0" << std::endl;
int e = 0;
Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
std::for_each(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
	std::cout << "index=" << index;
	std::cout << " metricPrims=";
	for (int j = 0; j < 10; ++j, ++e) {
		std::cout << " " << x[e];
	}
	
	std::cout << " G_ab=";
	real_SymDim EinsteinLL = calc_EinsteinLL(index, gLLs, gUUs, GammaULLs);
	real* p = &EinsteinLL(0,0);
	for (int j = 0; j < 10; ++j, ++p) {
		std::cout << " " << *p;
	}

	std::cout << " 8piT_ab=";
	real_SymDim _8piT_LL = calc_8piTLL(Tensor::Grid<MetricPrims, subDim>(sizev, (MetricPrims const *)x)(index), gLLs(index), stressEnergyPrimGrid(index));
	for (int a = 0; a < dim; ++a) {
		for (int b = 0; b <= a; ++b) {
			std::cout << " " << _8piT_LL(a,b);
		}
	}
	
	std::cout << std::endl;
});
#endif	//debug output
			},
			1e-100, 				//newton stop epsilon
			maxiter, 			//newton max iter
			[&](size_t n, real* x, real* b, JFNK::Func A) -> std::shared_ptr<Solver::Krylov<real>> {
				return std::make_shared<GMRES>(
					n, x, b, A,
					1e-100,	 				//gmres stop epsilon
					n, //n*10,				//gmres max iter
					gmresRestart			//gmres restart iter
				);
			}
		);
		jfnk.jacobianEpsilon = 1e-10;
		jfnk.maxAlpha = 1;
		//jfnk.lineSearch = &JFNK::lineSearch_none;
		jfnk.lineSearch = &JFNK::lineSearch_bisect;
#ifdef CONVERGE_ALPHA_ONLY	
		jfnk.lineSearchMaxIter = 50;
#else
		jfnk.lineSearchMaxIter = 20;
#endif
		jfnk.stopCallback = [&]()->bool{
			
#ifdef PRINT_RANGES
			real_SymDim gLL_mins, gLL_maxs;
			real_SymDim gUU_mins, gUU_maxs;
			real_Dim_SymDim GammaULL_mins, GammaULL_maxs;
			real_SymDim EFE_mins, EFE_maxs;
			//TODO store EinsteinLL and _8piT_LL so I can see ranges on those too
			for (int a = 0; a < dim; ++a) {
				for (int b = 0; b <= a; ++b) {
					gLL_maxs(a,b) = -(gLL_mins(a,b) = std::numeric_limits<real>::infinity());
					gUU_maxs(a,b) = -(gUU_mins(a,b) = std::numeric_limits<real>::infinity());
					EFE_maxs(a,b) = -(EFE_mins(a,b) = std::numeric_limits<real>::infinity());
				}
			}
			for (int a = 0; a < dim; ++a) {
				for (int b = 0; b < dim; ++b) {
					for (int c = 0; c <= b; ++c) {
						GammaULL_maxs(a,b,c) = -(GammaULL_mins(a,b,c) = std::numeric_limits<real>::infinity());
					}
				}
			}
			for (int k = 0; k < gridVolume; ++k) {
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b <= a; ++b) {
						real d = gLLs.v[k](a,b);
						gLL_mins(a,b) = std::min(gLL_mins(a,b), d);
						gLL_maxs(a,b) = std::max(gLL_maxs(a,b), d);
						
						d = gUUs.v[k](a,b);
						gUU_mins(a,b) = std::min(gUU_mins(a,b), d);
						gUU_maxs(a,b) = std::max(gUU_maxs(a,b), d);
						
						d = EFEGrid.v[k](a,b);
						EFE_mins(a,b) = std::min(EFE_mins(a,b), d);
						EFE_maxs(a,b) = std::max(EFE_maxs(a,b), d);
					}
				}
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						for (int c = 0; c <= b; ++c) {
							real d = GammaULLs.v[k](a,b,c);
							GammaULL_mins(a,b,c) = std::min(GammaULL_mins(a,b,c), d);
							GammaULL_maxs(a,b,c) = std::max(GammaULL_maxs(a,b,c), d);
						}
					}
				}
			}
			std::cout << "gLL range: " << std::endl;
			std::cout << " mins " << gLL_mins << std::endl;
			std::cout << " maxs " << gLL_maxs << std::endl;
			std::cout << "gUU range: " << std::endl;
			std::cout << " mins " << gUU_mins << std::endl;
			std::cout << " maxs " << gUU_maxs << std::endl;
			std::cout << "GammaULL range: " << std::endl;
			std::cout << " mins " << GammaULL_mins << std::endl;
			std::cout << " maxs " << GammaULL_maxs << std::endl;
			std::cout << "EFE range: " << std::endl;
			std::cout << " mins " << EFE_mins << std::endl;
			std::cout << " maxs " << EFE_maxs << std::endl;
	
#endif

			
			std::cout << "jfnk iter=" << jfnk.getIter() 
				<< " alpha=" << std::setprecision(49) << jfnk.getAlpha() << std::setprecision(6)
				<< " residual=" << std::setprecision(49) << jfnk.getResidual() << std::setprecision(6) 
				<< std::endl;
			
			jfnkFile << jfnk.getIter() 
				<< "\t" << std::setprecision(16) << jfnk.getResidual() << std::setprecision(6)
				<< "\t" << std::setprecision(49) << jfnk.getAlpha() << std::setprecision(6)
				<< std::endl;
			gmresFile << std::endl;
			
			return false;
		};
		std::shared_ptr<GMRES> gmres = std::dynamic_pointer_cast<GMRES>(jfnk.getLinearSolver());
		real lastResidual;
		gmres->stopCallback = [&]()->bool{
			if (gmres->getIter() > (int)jfnk.getN()) {
				if (gmres->getResidual() == lastResidual) {
					std::cout << "gmres stuck -- aborting gmres" << std::endl;
					return true;
				}
			}
			lastResidual = gmres->getResidual();
			
			if (!std::isfinite(gmres->getResidual())) {
				std::cout << "gmres got a non-finite residual -- aborting gmres" << std::endl;
				return true;
			}
			
			//the residual is staying constant ... at 16 even, for a 4*4*4 grid ...
			std::cout << "gmres"
				<< " iter=" << gmres->getIter() 
				<< " residual=" << std::setprecision(49) << gmres->getResidual() << std::setprecision(6)
				<< std::endl;
			
			gmresFile << jfnk.getIter()
				<< "\t" << gmres->getIter()
				<< "\t" << std::setprecision(16) << gmres->getResidual() << std::setprecision(6)
				<< std::endl;

			
			return false;
		};
//I don't think I've ever tested preconditioners ...
#if 0
		gmres->MInv = [&](real* y, real const * x) {
			for (int i = 0; i < (int)jfnk.getN(); ++i) {
				y[i] = x[i] / (8. * M_PI) * c * c / G / 1000.;
			}
		};
#endif
		time("solving", [&](){
			jfnk.solve();
		});

#ifdef CONVERGE_ALPHA_ONLY
		for (int i = 0; i < gridVolume; ++i) {
			metricPrimGrid.v[i].alphaMinusOne = alphaMinusOnes[i];
			//scale back down the alphaMinusOnes now that we're done solving
			metricPrimGrid.v[i].alphaMinusOne /= jfnkInputScale;
		}
#endif

		jfnkFile.close();
		gmresFile.close();
	}
};

struct Body {
	real radius;

	Body(real radius_) : radius(radius_) {}

	virtual void initStressEnergyPrim(
		Tensor::Grid<StressEnergyPrims, subDim> & stressEnergyPrimGrid,
		Tensor::Grid<Tensor::Vector<real, subDim>, subDim> const & xs
	) = 0;
};

struct NullBody : public Body {
	using Body::Body;

	virtual void initStressEnergyPrim(
		Tensor::Grid<StressEnergyPrims, subDim>& stressEnergyPrimGrid,
		Tensor::Grid<Tensor::Vector<real, subDim>, subDim> const & xs
	) { }
};

struct SphericalBody : public Body {
	real mass;
	//derived:
	real volume, density;
	
	SphericalBody(real radius_, real mass_)
	: Body(radius_),
	mass(mass_),
	volume(4./3.*M_PI*radius*radius*radius),	// m^3
	density(mass / volume)	// 1/m^2
	{}

	virtual void initStressEnergyPrim(
		Tensor::Grid<StressEnergyPrims, subDim>& stressEnergyPrimGrid,
		Tensor::Grid<Tensor::Vector<real, subDim>, subDim> const & xs
	) {
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			StressEnergyPrims &stressEnergyPrims = stressEnergyPrimGrid(index);
			real r = xs(index).length();
			stressEnergyPrims.rho = r < radius ? density : 0;	// average density of Earth in m^-2
			stressEnergyPrims.eInt = 0;	//internal energy / temperature of the Earth?
			stressEnergyPrims.P = 0;	//pressure inside the Earth?
			for (int i = 0; i < subDim; ++i) {
				stressEnergyPrims.v(i) = 0;	//3-velocity
				stressEnergyPrims.E(i) = 0;	//electric field
				stressEnergyPrims.B(i) = 0;	//magnetic field
			}
		});
	}
};

//E_i = A_t,i - A_i,t
//B_i = epsilon_i^jk A_k,j

//stress-energy of a uniform electric field in the x-direction
struct EMUniformFieldBody : public Body {
	using Body::Body;
	virtual void initStressEnergyPrim(
		Tensor::Grid<StressEnergyPrims, subDim>& stressEnergyPrimGrid,
		Tensor::Grid<Tensor::Vector<real, subDim>, subDim> const & xs
	) {
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			StressEnergyPrims &stressEnergyPrims = stressEnergyPrimGrid(index);

			stressEnergyPrims.useEM = true; 

			real volts = 100000;			//V
			real dist = 0.1;				//m
			real E = volts / dist;			//V / m
			E = E * sqrt(ke * G) / (c * c);	//m^-1

			stressEnergyPrims.E(0) = E;
			stressEnergyPrims.E(1) = 0;
			stressEnergyPrims.E(2) = 0;
	
			stressEnergyPrims.B(0) = 0;
			stressEnergyPrims.B(1) = 0;
			stressEnergyPrims.B(2) = 0;
		});
	}
};

struct EMLineBody : public Body {
	using Body::Body;

	//torus-shaped-something
	//radius is the big radius of the torus
	virtual void initStressEnergyPrim(
		Tensor::Grid<StressEnergyPrims, subDim>& stressEnergyPrimGrid,
		Tensor::Grid<Tensor::Vector<real, subDim>, subDim> const & xs
	) {
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			StressEnergyPrims &stressEnergyPrims = stressEnergyPrimGrid(index);
			Tensor::Vector<real, subDim> const & xi = xs(index);
			real x = xi(0);
			real y = xi(1);
			real z = xi(2);
			real polar_rSq = x*x + y*y; 
			real polar_r = sqrt(polar_rSq);		//r in polar coordinates	
			real dr = polar_r - radius;	//difference from polar radius to torus big radius
			real r = sqrt(z*z + dr*dr);			//r in torus radial coordinates
			real theta = atan2(z, dr);			//angle around the small radius
			real phi = atan2(x, y);				//angle around the big radius
		
			//F^uv_;v = -4 pi J^u
			// means that the divergence of the EM is the 4-current 
			//the divergence of the exterior of the 4-potential is the 4-current
			//so if the 4-current is a Dirac delta function along the line in space where there is current
			// then the EM tensor is going to be an inverse falloff around it
	
			//4-current
			//t is current density
			//i is charge density = current density * drift velocity
			//stressEnergyPrims.chargeDensity
			//stressEnergyPrims.A
	
			/*
			point on the surface:
				r * cos(phi) * cos(theta)
				r * sin(phi) * cos(theta)
				r * sin(theta)
			*/
		
			stressEnergyPrims.E(0) = -y / polar_rSq;
			stressEnergyPrims.E(1) = x / polar_rSq;
			stressEnergyPrims.E(2) = 0;
	
			stressEnergyPrims.B(0) = cos(theta) / r * cos(phi);
			stressEnergyPrims.B(1) = cos(theta) / r * sin(phi);
			stressEnergyPrims.B(2) = -sin(theta) / r;
		});
	}
};

struct InitCond {
	virtual void initMetricPrims(
		Tensor::Grid<MetricPrims, subDim>& metricPrimGrid,
		Tensor::Grid<Tensor::Vector<real, subDim>, subDim> const & xs
	) = 0;
};

struct SphericalBodyInitCond : public InitCond {
	std::shared_ptr<SphericalBody> body;	
	SphericalBodyInitCond(std::shared_ptr<SphericalBody> body_) : body(body_) {}
};

struct FlatInitCond : public InitCond {
	//substitute the schwarzschild R for 2 m(r)
	virtual void initMetricPrims(
		Tensor::Grid<MetricPrims, subDim>& metricPrimGrid,
		Tensor::Grid<Tensor::Vector<real, subDim>, subDim> const & xs
	) {
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			MetricPrims& metricPrims = metricPrimGrid(index);
			metricPrims.alphaMinusOne = 0;
			for (int i = 0; i < subDim; ++i) {
				metricPrims.betaU(i) = 0;
				for (int j = 0; j <= i; ++j) {
					metricPrims.hLL(i,j) = 0;
				}
			}
		});
	}
};

/*
The stellar Schwarzschild initial conditions have constraint value of zero outside Earth (good)
but inside Earth they give a difference of 2 g/cm^3 ... off from the Earth's density of 5.51 g/cm^3
*/
struct StellarSchwarzschildInitCond : public SphericalBodyInitCond {
	using SphericalBodyInitCond::SphericalBodyInitCond;
	virtual void initMetricPrims(
		Tensor::Grid<MetricPrims, subDim>& metricPrimGrid,
		Tensor::Grid<Tensor::Vector<real, subDim>, subDim> const & xs
	) {
		real radius = body->radius;
		real density = body->density;
		real mass = body->mass;
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			MetricPrims& metricPrims = metricPrimGrid(index);
			Tensor::Vector<real, subDim> const & xi = xs(index);
			real r = xi.length();
			real matterRadius = std::min<real>(r, radius);
			real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
			real m = density * volumeOfMatterRadius;	// m^3		
			
			/*
			g_ti = beta_i = 0
			g_tt = -alpha^2 + beta^2 = -alpha^2 = -1 + Rs/r <=> alpha = sqrt(1 - Rs/r)
			g_ij = gamma_ij = delta_ij + x^i x^j / r^2 2M/(r - 2M)		<- but x is upper, and you can't lower it without specifying gamma_ij
			 ... which might be why the contravariant spatial metrics of spherical and cartesian look so similar 
			*/
			/*
			I'm going by MTW box 23.2 eqn 6 d/dt (proper time) = sqrt(1 - R/r) for r > R
				and ( 3/2 sqrt(1 - 2 M / R) - 1/2 sqrt(1 - 2 M r^2 / R^3) ) for r < R
				for M = total mass 
				and R = planet radius 
			*/
			metricPrims.alphaMinusOne = (
					r > radius 
					? sqrt(1 - 2*mass/r)
					: (1.5 * sqrt(1 - 2*mass/radius) - .5 * sqrt(1 - 2*mass*r*r/(radius*radius*radius)))
				) - 1.;
			
			for (int i = 0; i < subDim; ++i) {
				metricPrims.betaU(i) = 0;
				for (int j = 0; j <= i; ++j) {
					metricPrims.hLL(i,j) = xi(i)/r * xi(j)/r * 2*m/(r - 2*m);
					/*
					dr^2's coefficient
					spherical: 1/(1 - 2M/r) = 1/((r - 2M)/r) = r/(r - 2M)
					spherical contravariant: 1 - 2M/r
					cartesian contravariant: delta_ij - x/r y/r 2M/r
					hmm, contravariant terms of cartesian vs spherical look more similar than covariant terms do
				
					in the OV metric, dr^2's coefficient is exp(2 Lambda) = 1/(1 - 2 m(r) / r) where m(r) is the enclosing mass
					so the contravariant coefficient would be exp(-2 Lambda) = 1 - 2 m(r) / r
					I'm going to do the lazy thing and guess this converts to delta^ij - 2 m(r) x^i x^j / r^3
					*/
				}
			}

#if 0	//rotating about a distance
			/*
			now if we are going to rotate this
			at a distance of L and at an angular frequency of omega
			(not considering relativistic Thomas precession just yet)
			
			this might be a mess, but I'm (1) calculating the change in time as if I were in a frame rotating by L exp(i omega t)
			then (2) re-centering the frame at L exp(i omega t) ... so I can use the original coordinate system
			*/
			real dr_alpha = r > radius 
				? (mass / (r * r * sqrt(1. - 2. * mass / r))) 
				: (mass * r / (radius * radius * radius * sqrt(1. - 2. * mass * r * r / (radius * radius * radius))));
			real dr_m = r > radius ? 0 : (4. * M_PI * r * r * density);
			MetricPrims& dt_metricPrims = dt_metricPrimGrid(index);
			real L = 149.6e+9;	//distance from earth to sun, in m 
			//real omega = 0; //no rotation
			//real omega = 2. * M_PI / (60. * 60. * 24. * 365.25) / c;	//one revolution per year in m^-1 
			//real omega = 1;	//angular velocity of the speed of light
			real omega = c;	//I'm trying to find a difference ...
			real t = 0;	//where the position should be.  t=0 means the body is moved by [L, 0], and its derivatives are along [0, L omega] 
			Tensor::Vector<real,2> dt_xHat(L * omega * sin(omega * t), -L * omega * cos(omega * t));
			dt_metricPrims.alphaMinusOne = dr_alpha * (xi(0)/r * dt_xHat(0) + xi(1)/r * dt_xHat(1)) - 1.;
			for (int i = 0; i < subDim; ++i) {
				dt_metricPrims.betaU(i) = 0;
			}
			for (int i = 0; i < subDim; ++i) {
				for (int j = 0; j < subDim; ++j) {
					real sum = 0;
					for (int k = 0; k < 2; ++k) {
						//gamma_ij = f/g
						//so d/dxHat^k gamma_ij = 
						real dxHat_k_of_gamma_ij = 
						// f' / g
						(
							((i==k)*xi(j) + xi(i)*(j==k)) * 2.*m + xi(i)*xi(j) * 2.*dr_m * xi(k)/r
						) / (r * r * (r - 2 * m))
						// - f g' / g^2
						- (xi(i) * xi(j) * 2 * m) * ( (xi(k) - 2 * dr_m * xi(k)) * r + 2 * xi(k) * (r - 2 * m) )
						/ (r * r * r * r * (r - 2 * m) * (r - 2 * m));
						sum += dxHat_k_of_gamma_ij * dt_xHat(k);
					}
					dt_metricPrims.hLL(i,j) = sum - (real)(i == j);
				}
			}
#endif

#if 0	//work that beta
			/*
			so if we can get the gamma^ij beta_j,t components to equal the Gamma^t_tt components ...
			 voila, gravity goes away.
			I'm approximating this as beta^i_,t ... but it really is beta^j_,t gamma_ij + beta^j gamma_ij,t
			 ... which is the same as what I've got, but I'm setting gamma_ij,t to zero
			*/
			//expanding ...
			MetricPrims& dt_metricPrims = dt_metricPrimGrid(index);
			real_Sub betaL;
			for (int i = 0; i < subDim; ++i) {
				//negate all gravity by throttling the change in space/time coupling of the metric
				real dm_dr = 0;
				betaL(i) = -(2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r) * xi(i)/r;
			}
			real_SymSub gammaUU = Tensor::v2::inverse(metricPrims.hLL + delta3LL);
			for (int i = 0; i < subDim; ++i) {
				real sum = 0;
				for (int j = 0; j < subDim; ++j) {
					sum += gammaUU(i,j) * betaL(j);
				}
				dt_metricPrims.betaU(i) = sum;
			}
#endif
		});
	}
};

struct StellarKerrNewmanInitCond : public SphericalBodyInitCond {
	using SphericalBodyInitCond::SphericalBodyInitCond;
	virtual void initMetricPrims(
		Tensor::Grid<MetricPrims, subDim>& metricPrimGrid,
		Tensor::Grid<Tensor::Vector<real, subDim>, subDim> const & xs
	) {
		real radius = body->radius;
		real mass = body->mass;
		real density = body->density;
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			MetricPrims& metricPrims = metricPrimGrid(index);
			Tensor::Vector<real,subDim> const & xi = xs(index);
			
			real x = xi(0);
			real y = xi(1);
			real z = xi(2);
			
			real angularVelocity = 2. * M_PI / (60. * 60. * 24.) / c;	//angular velocity, in m^-1
			real inertia = 2. / 5. * mass * radius * radius;	//moment of inertia about a sphere, in m^3
			real angularMomentum = inertia * angularVelocity;	//angular momentum in m^2
			real a = angularMomentum / mass;	//m
			
			//real r is the solution of (x*x + y*y) / (r*r + a*a) + z*z / (r*r) = 1 
			// r^4 - (x^2 + y^2 + z^2 - a^2) r^2 - a^2 z^2 = 0
			real RSq_minus_aSq = x*x + y*y + z*z - a*a;
			//so we have two solutions ... which do we use? 
			//from gnuplot it looks like the two of these are the same ...
			real r = sqrt((RSq_minus_aSq + sqrt(RSq_minus_aSq * RSq_minus_aSq + 4.*a*a*z*z)) / 2.);	//use the positive root

			//should I use the Kerr-Schild 'r' coordinate?
			//well, if 'm' is the mass enclosed within the coordinate
			// and that determines 'a', the angular momentum per mass within the coordinate (should it?)
			// then we would have a circular definition
			//real R = sqrt(x*x + y*y + z*z); 
			real matterRadius = std::min<real>(r, radius);
			real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
			real m = density * volumeOfMatterRadius;	// m^3

			real Q = 0;	//charge
			real H = (r*m - Q*Q/2.)/(r*r + a*a*z*z/(r*r));
		
			//3.4.33 through 3.4.35 of Alcubierre "Introduction to 3+1 Numerical Relativity"
			
			/*TODO fix this for the metric within the star
			 in other news, this is an unsolved problem!
			https://arxiv.org/pdf/1503.02172.pdf section 3.11
			https://arxiv.org/pdf/1410.2130.pdf section 4.2 last paragraph
			*/
			//metricPrims.alphaMinusOne = 1./sqrt(1. + 2*H) - 1.;
			metricPrims.alphaMinusOne = sqrt(1. - 2*H/(1+2*H) ) - 1.;
			
			Tensor::Vector<real,subDim> l( (r*x + a*y)/(r*r + a*a), (r*y - a*x)/(r*r + a*a), z/r );
			for (int i = 0; i < subDim; ++i) {
				metricPrims.betaU(i) = 2. * H * l(i) / (1. + 2. * H);
				for (int j = 0; j <= i; ++j) {
					metricPrims.hLL(i,j) = 2 * H * l(i) * l(j); 
				}
			}
		});
	}
};

struct EMUniformFieldInitCond : public FlatInitCond {
	std::shared_ptr<EMUniformFieldBody> body;	
	EMUniformFieldInitCond(std::shared_ptr<EMUniformFieldBody> body_) : body(body_) {}
	//use flat initial metric prims for now
};

struct EMLineInitCond : public FlatInitCond {
	std::shared_ptr<EMLineBody> body;
	EMLineInitCond(std::shared_ptr<EMLineBody> body_) : body(body_) {}
};

int main(int argc, char** argv) {	

	//std::cout << "mass=" << mass << std::endl;		//m
	//std::cout << "radius=" << radius << std::endl;	//m
	//std::cout << "volume=" << volume << std::endl;	//m^3
	//std::cout << "density=" << density << std::endl;	//m^-2

	LuaCxx::State lua;
	lua.loadFile("config.lua");
	
	int maxiter = std::numeric_limits<int>::max();
	if (!lua["maxiter"].isNil()) lua["maxiter"] >> maxiter;	
	std::cout << "maxiter=" << maxiter << std::endl;

	std::string bodyName = "earth";
	if (!lua["body"].isNil()) lua["body"] >> bodyName;
	std::cout << "body=\"" << bodyName << "\"" << std::endl;

	std::string initCondName = "stellar_schwarzschild";
	if (!lua["initCond"].isNil()) lua["initCond"] >> initCondName;
	std::cout << "initCond=\"" << initCondName << "\"" << std::endl;

	std::string solverName = "jfnk";
	if (!lua["solver"].isNil()) lua["solver"] >> solverName;	
	std::cout << "solver=\"" << solverName << "\"" << std::endl;

	sizev = Tensor::Vector<int, subDim>(16, 16, 16);
	if (!lua["size"].isNil()) {
		if (lua["size"].isNumber()) {
			for (int i = 0; i < subDim; ++i) {
				lua["size"] >> sizev(i);
			}
		} else if (lua["size"].isTable()) {
			for (int i = 0; i < subDim; ++i) {
				if (!lua["size"][i+1].isNumber()) throw Common::Exception() << "size[" << (i+1) << "] is not a number";
				lua["size"][i+1] >> sizev(i);
			}
		}
	}
	std::cout << "size=" << sizev << std::endl;

	real bodyRadii = 2;
	if (!lua["bodyRadii"].isNil()) {
		double d = bodyRadii; lua["bodyRadii"] >> d; bodyRadii = d;
	}
	std::cout << "bodyRadii=" << bodyRadii << std::endl;


	std::shared_ptr<Body> body;
	{
		struct {
			char const * name;
			std::function<std::shared_ptr<Body>()> func;
		} bodies[] = {
			{"Null", [&](){
				return std::make_shared<NullBody>(2);
			}},
			
			{"earth", [&](){
				real const earthRadius = 6.37101e+6;	// m
				real const earthMass = 5.9736e+24 * G / c / c;	// m
				//earth volume: 1.0832120174985e+21 m^3
				//earth density: 4.0950296770075e-24 1/m^2 = 5.5147098661212 g/cm^3, which is what Google says.
				//note that G_tt = 8 pi T_tt = 8 pi rho ... for earth = 1.0291932119615e-22 m^-2
				//real const schwarzschildRadius = 2 * mass;	//Schwarzschild radius: 8.87157 mm, which is accurate
				//earth magnetic field at surface: .25-.26 gauss
				//real const earthMagneticField = .45 * sqrt(.1 * G) / c;	// 1/m
				return std::make_shared<SphericalBody>(earthRadius, earthMass);
			}},
		
			{"sun", [&](){
				real const sunRadius = 6.960e+8;	// m
				real const sunMass = 1.9891e+30 * G / c / c;	// m
				return std::make_shared<SphericalBody>(sunRadius, sunMass);
			}},
		
			{"EMUniformField", [&](){
				return std::make_shared<EMUniformFieldBody>(2);
			}},
			
			{"em_line", [&](){
				return std::make_shared<EMLineBody>(2);
			}},
		}, *p;
	
		for (p = bodies; p < endof(bodies); ++p) {
			if (p->name == bodyName) {
std::cout << "creating body " << bodyName << std::endl;				
				body = p->func();
				break;
			}
		}
		if (!body) {
			throw Common::Exception() << "couldn't find body named " << bodyName;
		}
	}


	xmin = Tensor::Vector<real, subDim>(-bodyRadii*body->radius, -bodyRadii*body->radius, -bodyRadii*body->radius),
	xmax = Tensor::Vector<real, subDim>(bodyRadii*body->radius, bodyRadii*body->radius, bodyRadii*body->radius);
	gridVolume = sizev.volume();
	dx = (xmax - xmin) / sizev;

	Tensor::Grid<Tensor::Vector<real, subDim>, subDim> xs;
	Tensor::Grid<MetricPrims, subDim> metricPrimGrid;
	Tensor::Grid<MetricPrims, subDim> dt_metricPrimGrid;	//first deriv
	Tensor::Grid<StressEnergyPrims, subDim> stressEnergyPrimGrid;

	time("allocating", [&]{ 
		std::cout << std::endl;
		size_t totalSize = 0;
#define ALLOCATE_GRID(x)	allocateGrid(x, #x, sizev, totalSize)
		ALLOCATE_GRID(xs);
		ALLOCATE_GRID(metricPrimGrid);
		ALLOCATE_GRID(dt_metricPrimGrid);	//first deriv
		ALLOCATE_GRID(stressEnergyPrimGrid);
		ALLOCATE_GRID(gLLs);
		ALLOCATE_GRID(gUUs);
		ALLOCATE_GRID(dt_gLLs);
		ALLOCATE_GRID(d2t_gLLs);
		//ALLOCATE_GRID(dt_gUUs);
		ALLOCATE_GRID(dgLLLs);
		//ALLOCATE_GRID(GammaLLLs);
		ALLOCATE_GRID(GammaULLs);
#undef ALLOCATE_GRID
	});

	//specify coordinates
	time("calculating grid", [&]{
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			Tensor::Vector<real, subDim>& xi = xs(index);
			for (int j = 0; j < subDim; ++j) {
				xi(j) = (xmax(j) - xmin(j)) * ((real)index(j) + .5) / (real)sizev(j) + xmin(j);
			}
		});
	});

	//specify stress-energy primitives
	//the stress-energy primitives combined with the current metric are used to compute the stress-energy tensor 
	//this is done by choosing the 'body'

	//initialize metric primitives
	time("calculating stress-energy primitives", [&]{
		body->initStressEnergyPrim(stressEnergyPrimGrid, xs);
	});

	//while we're here, set the 'useE' and 'useV' flags, to spare our calculations
	time("determine what stress-energy variables to use", [&]() {
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			StressEnergyPrims& stressEnergyPrims = stressEnergyPrimGrid(index);
			
			stressEnergyPrims.useV = false;
			for (int i = 0; i < subDim; ++i) {
				stressEnergyPrims.useV |= stressEnergyPrims.v(i) != 0;
			}
			
			stressEnergyPrims.useEM = false;
#ifdef USE_CHARGE_CURRENT_FOR_EM
			stressEnergyPrims.useEM |= stressEnergyPrims.chargeDensity != 0;
			for (int i = 0; i < subDim; ++i) {
				stressEnergyPrims.useEM |= stressEnergyPrims.currentDensity(i) != 0;
			}
#else
			for (int i = 0; i < subDim; ++i) {
				stressEnergyPrims.useEM |= stressEnergyPrims.E(i) != 0;
				stressEnergyPrims.useEM |= stressEnergyPrims.B(i) != 0;
			}
#endif
		});
	});

	{
		struct {
			char const * name;
			std::function<std::shared_ptr<InitCond>()> func;
		} initConds[] = {
			{"flat", [&](){ return std::make_shared<FlatInitCond>(); }},
			{"stellar_schwarzschild", [&](){ 
				std::shared_ptr<SphericalBody> sphericalBody = std::dynamic_pointer_cast<SphericalBody>(body);
				assert(sphericalBody);
				return std::make_shared<StellarSchwarzschildInitCond>(sphericalBody);
			}},
			{"stellar_kerr_newman", [&](){ 
				std::shared_ptr<SphericalBody> sphericalBody = std::dynamic_pointer_cast<SphericalBody>(body);
				assert(sphericalBody);
				return std::make_shared<StellarKerrNewmanInitCond>(sphericalBody);
			}},
			{"EMUniformField", [&](){ 
				std::shared_ptr<EMUniformFieldBody> emFieldBody = std::dynamic_pointer_cast<EMUniformFieldBody>(body);
				assert(emFieldBody);
				return std::make_shared<EMUniformFieldInitCond>(emFieldBody);
			}},
			{"em_line", [&](){ 
				std::shared_ptr<EMLineBody> emLineBody = std::dynamic_pointer_cast<EMLineBody>(body);
				assert(emLineBody);
				return std::make_shared<EMLineInitCond>(emLineBody);
			}},
		}, *p;

		std::shared_ptr<InitCond> initCond;
		for (p = initConds; p < endof(initConds); ++p) {
			if (p->name == initCondName) {
				initCond = p->func();
				break;
			}
		}
		if (!initCond) {
			throw Common::Exception() << "couldn't find initial condition named " << initCondName;
		}

		//initialize metric primitives
		time("calculating metric primitives", [&]{
			initCond->initMetricPrims(metricPrimGrid, xs);
		});
	}

	std::shared_ptr<EFESolver> solver;
	{
		struct {
			char const * name;
			std::function<std::shared_ptr<EFESolver>()> func;
		} solvers[] = {
			{"jfnk", [&](){ return std::make_shared<JFNKSolver>(maxiter); }},
			{"gmres", [&](){ return std::make_shared<GMRESSolver>(maxiter); }},
			{"conjres", [&](){ return std::make_shared<ConjResSolver>(maxiter); }},
			{"conjgrad", [&](){ return std::make_shared<ConjGradSolver>(maxiter); }},
		}, *p;
		for (p = solvers; p < endof(solvers); ++p) {
			if (p->name == solverName) {
				solver = p->func();
			}
		}
		if (!solver) {
			throw Common::Exception() << "couldn't find solver named " << solverName;
		}
	}

	if (maxiter > 0) {
		solver->solve(
			metricPrimGrid, 
			dt_metricPrimGrid, 	//first deriv
			stressEnergyPrimGrid);
	}

	//once all is solved for, do some final calculations ...

	time("calculating g_ab and g^ab", [&]{
		calc_gLLs_and_gUUs(
			//input:
			metricPrimGrid,
			dt_metricPrimGrid,	//first deriv
			//output:
			gLLs, gUUs, dt_gLLs);
	});

	time("calculating Gamma^a_bc", [&]{
		calc_GammaULLs(gLLs, gUUs, dt_gLLs, GammaULLs);
	});

	Tensor::Grid<real_SymDim, subDim> EFEGrid(sizev);
	time("calculating EFE constraint", [&]{
		calc_EFE_constraint(metricPrimGrid, stressEnergyPrimGrid, EFEGrid);
	});

	Tensor::Grid<real, subDim> numericalGravity(sizev);
	time("calculating numerical gravitational force", [&]{
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			Tensor::Vector<real, subDim> xi = xs(index);
			real r = xi.length();
			//numerical computational...		
			real_Dim_SymDim & GammaULL = GammaULLs(index);
//the analytical calculations on these are identical, provided Gamma^i_tt is the schwarzschild metric connection
//but the acceleration magnitude method can't show sign
#if 1 // here's change-of-coordinate from G^i_tt to G^r_tt
			//Gamma^r_tt = Gamma^i_tt dr/dx^i
			//r^2 = x^2 + y^2 + z^2
			//so dr/dx^i = x^i / r
			numericalGravity(index) = (GammaULL(1,0,0) * xi(0) / r
									+ GammaULL(2,0,0) * xi(1) / r
									+ GammaULL(3,0,0) * xi(2) / r)
									* c * c;	//times c twice because of the two timelike components of a^i = Gamma^i_tt
#endif
#if 0	//here's taking the acceleration in cartesian and computing the magnitude of the acceleration vector
			numericalGravity(index) = sqrt(GammaULL(1,0,0) * GammaULL(1,0,0)
										+ GammaULL(2,0,0) * GammaULL(2,0,0)
										+ GammaULL(3,0,0) * GammaULL(3,0,0))
										* c * c;	//times c twice because of the two timelike components of a^i = Gamma^i_tt
#endif
		});
	});

	Tensor::Grid<real, subDim> analyticalGravity(sizev);
	std::shared_ptr<SphericalBody> sphericalBody = std::dynamic_pointer_cast<SphericalBody>(body);
	if (sphericalBody) {
		time("calculating analytical gravitational force", [&]{
			Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
			parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
				Tensor::Vector<real, subDim> xi = xs(index);
				real r = xi.length();
				//substitute the schwarzschild R for 2 m(r)
				real matterRadius = std::min<real>(r, sphericalBody->radius);
				real volumeOfMatterRadius = 4./3.*M_PI*matterRadius*matterRadius*matterRadius;
				real m = sphericalBody->density * volumeOfMatterRadius;	// m^3

				//now that I'm using the correct alpha equation, my dm/dr term is causing the analytical gravity calculation to be off ...
				//real dm_dr = r > radius ? 0 : density * 4 * M_PI * matterRadius * matterRadius;
				// ... maybe it shouldn't be there to begin with?
				real dm_dr = 0;
				real GammaUr_tt = (2*m * (r - 2*m) + 2 * dm_dr * r * (2*m - r)) / (2 * r * r * r)
					* c * c;	//+9 at earth surface, without matter derivatives

				//acceleration is -Gamma^r_tt along the radial direction (i.e. upwards from the surface), or Gamma^r_tt downward into the surface
				analyticalGravity(index) = GammaUr_tt;
			});
		});
	}

	{
		struct Col {
			std::string name;
			std::function<real(Tensor::Vector<int,subDim>)> func;
		};
		std::vector<Col> cols = {
			{"ix", [&](Tensor::Vector<int,subDim> index)->real{ return index(0); }},
			{"iy", [&](Tensor::Vector<int,subDim> index)->real{ return index(1); }},
			{"iz", [&](Tensor::Vector<int,subDim> index)->real{ return index(2); }},
			{"rho", [&](Tensor::Vector<int,subDim> index)->real{ return stressEnergyPrimGrid(index).rho; }},
			{"det_h", [&](Tensor::Vector<int,subDim> index)->real{ return Tensor::determinant33<real, real_SymSub>(metricPrimGrid(index).hLL); }},
			{"alpha-1", [&](Tensor::Vector<int,subDim> index)->real{ return metricPrimGrid(index).alphaMinusOne; }},
#if 0	//within 1e-23			
			{"ortho_error", [&](Tensor::Vector<int,subDim> index)->real{
				real_SymDim const & gLL = gLLs(index);
				real_SymDim const & gUU = gUUs(index);
				real err = 0;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						real sum = 0;
						for (int c = 0; c < dim; ++c) {
							sum += gLL(a,c) * gUU(c,b);
						}
						err += fabs(sum - (real)(a == b));
					}
				}
				return err;
			}},
#endif
#if 1		// numerical gravity is double what analytical gravity is ... and flips to negative as it passes the planet surface ...
			{"gravity", [&](Tensor::Vector<int,subDim> index)->real{ return numericalGravity(index); }},
		};
		if (sphericalBody) {
			cols.push_back({"analyticalGravity", [&](Tensor::Vector<int,subDim> index)->real{ return analyticalGravity(index); }});
		}
		std::vector<Col> moreCols = {
#endif
			{"EFE_tt(g/cm^3)", [&](Tensor::Vector<int,subDim> index)->real{
				return EFEGrid(index)(0,0) / (8. * M_PI) * c * c / G / 1000.;	// g/cm^3 ... so in absense of any curvature, the constraint error will now match the density
			}},
			{"EFE_ti", [&](Tensor::Vector<int,subDim> index)->real{ 
				real_SymDim &t = EFEGrid(index);
				return sqrt( t(0,1)*t(0,1) + t(0,2)*t(0,2) + t(0,3)*t(0,3) ) * c;
			}},
			{"EFE_ij", [&](Tensor::Vector<int,subDim> index) -> real {
				real_SymDim &t = EFEGrid(index);
				/* determinant
				return t(1,1) * t(2,2) * t(3,3)
					+ t(1,2) * t(2,3) * t(3,1)
					+ t(1,3) * t(2,1) * t(3,2)
					- t(1,3) * t(2,2) * t(3,1)
					- t(1,1) * t(2,3) * t(3,2)
					- t(1,2) * t(2,1) * t(3,3);
				*/
				// norm
				real sum = 0;
				for (int a = 1; a < dim; ++a) {
					for (int b = 1; b < dim; ++b) {
						sum += t(a,b)*t(a,b);
					}
				}
				return sqrt(sum);
			}},
#if 1
			{"G_ab", [&](Tensor::Vector<int,subDim> index)->real{
				real_SymDim G = calc_EinsteinLL(index, gLLs, gUUs, GammaULLs);
				real sum = 0;
				for (int a = 0; a < dim; ++a) {
					for (int b = 0; b < dim; ++b) {
						sum += G(a,b) * G(a,b);
					}
				}
				return sqrt(sum);
			}},
#endif
		};
		cols.insert(cols.end(), moreCols.begin(), moreCols.end());

		if (!lua["outputFilename"].isNil()) {
			std::string outputFilename;
			lua["outputFilename"] >> outputFilename;

			std::ofstream file(outputFilename);
			if (!file.good()) throw Common::Exception() << "failed to open file " << outputFilename;

			file << "#";
			{
				char const * tab = "";
				for (std::vector<Col>::iterator p = cols.begin(); p != cols.end(); ++p) {
					file << tab << p->name;
					tab = "\t";
				}
			}
			file << std::endl;
			time("outputting", [&]{
				//this is printing output, so don't do it in parallel		
				Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
				for (Tensor::RangeObj<subDim>::iterator iter = range.begin(); iter != range.end(); ++iter) {
					char const * tab = "";
					for (std::vector<Col>::iterator p = cols.begin(); p != cols.end(); ++p) {
						file << tab << std::setprecision(16) << p->func(iter.index) << std::setprecision(6);
						tab = "\t";
					}
					file << std::endl;
				}
			});
			
			file.close();
		}
	}

#if 0
	{
		Tensor::RangeObj<subDim> range(Tensor::Vector<int,subDim>(), sizev);
		
		real EFE_tt_min = std::numeric_limits<real>::infinity();
		real EFE_tt_max = -std::numeric_limits<real>::infinity();
		std::for_each(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			real EFE_tt = EFEGrid(index)(0,0);
			if (EFE_tt < EFE_tt_min) EFE_tt_min = EFE_tt;
			if (EFE_tt > EFE_tt_max) EFE_tt_max = EFE_tt;
		});
		std::cout << "EFE_tt range: " << EFE_tt_min << " to " << EFE_tt_max << std::endl;

		int bins = 256;
		std::vector<real> EFE_tt_distr(bins);
		std::for_each(range.begin(), range.end(), [&](Tensor::Vector<int, subDim> const & index) {
			real EFE_tt = EFEGrid(index)(0,0);
			int bin = (int)((EFE_tt - EFE_tt_min) / (EFE_tt_max - EFE_tt_min) * (real)bins);
			if (bin == bins) --bin;
			++EFE_tt_distr[bin];
		});

		std::cout << "EFE_tt:" << std::endl;
		for (int i = 0; i < bins; ++i) {
			real delta = EFE_tt_max - EFE_tt_min;
			std::cout << (delta * (real)i / (real)bins + EFE_tt_min) << "\t" 
					<< (delta * (real)(i+1) / (real)bins + EFE_tt_min) << "\t"
					<< EFE_tt_distr[i] << std::endl;
		}
	}
#endif

	std::cout << "done!" << std::endl;
}
