#include <iostream>
#include "Parallel/Parallel.h"
#include "Tensor/Grid.h"
#include "Solvers/GMRes.h"

Parallel::Parallel parallel(8);

const int gridDim = 3;
const int dim = 1;

typedef double real;
using namespace Tensor;

/*

J_a = (rho, j_i)

flat space:

F^uv_,v = -4 pi J^u
F_uv = A_v,u - A_u,v
A^v,u_v - A^u,v_v = -4 pi J^u

curved space:

A^a;u_u - A^u_;u^;a - R^a_u A^u = -4 pi J^a

use gauge A^u_;u = 0

A^a;u_u - R^a_u A^u = -4 pi J^a
A^a;u_;u - R^a_u A^u = -4 pi J^a

A^a;u_;u
(A^a_;b g^bu)_;u
((A^a_,b + Gamma^a_cb A^c) g^bu)_;u
(A^a,u + Gamma^au_c A^c)_;u
A^a,u_,u 
	+ (Gamma^au_c A^c)_,u
	+ Gamma^a_bu (A^b,u + Gamma^bu_c A^c)
	+ Gamma^u_bu (A^a,b + Gamma^ab_c A^c)
A^a,u_,u
	+ Gamma^a_cv,u g^vu A^c 
	+ Gamma^a_cv g^vu_,u A^c 
	+ Gamma^a_cv g^vu A^c_,u
	+ Gamma^a_bu A^b,u
	+ Gamma^u_bu A^a,b
	+ Gamma^a_bu Gamma^bu_c A^c
	+ Gamma^u_bu Gamma^ab_c A^c


A^a;u = A^a_;v g^uv = (A^a_,v + Gamma^a_wv A^w) g^uv
A^a;u_;u = A^a;u_,u + Gamma^a_bu A^b;u + Gamma^u_bu A^a;b
		= (A^a_;v g^uv = (A^a_,v + Gamma^a_wv A^w) g^uv)_,u
			+ Gamma^a_bu (A^b_,v + Gamma^b_wv A^w) g^uv
			- Gamma^u_bu (A^a_,v + Gamma^a_wv A^w) g^bv
((A^a_,b + Gamma^a_cb A^c) + R^a_b A^b) / (4 pi) = J^a

*/

int main() {
	size_t n = 50;
	Vector<int, gridDim> size;
	for (int i = 0; i < gridDim; ++i) {
		size(i) = n;
	}
	Grid<Vector<real, dim>, gridDim> JU(size);
	Grid<Vector<real, dim>, gridDim> AU(size);
	real dxSq = 1;

	//solve for AU
	RangeObj<gridDim> range = JU.range();
	parallel.foreach(range.begin(), range.end(), [&](const Vector<int, gridDim>& index) {
		Vector<real,gridDim> dx = ((Vector<real,gridDim>)index + .5) - size/2;
		real linfDist = 0;
		real sign = 1;
		for (int i = 0; i < gridDim; ++i) {
			linfDist = std::max<real>(linfDist, fabs(dx(i)));
			if (dx(i) < 0) sign = -sign;
		}

		const real r = 15;
		JU(index)(0) = linfDist < r ? sign : 0;
		AU(index) = JU(index);
	});

	Solvers::Krylov<real>::Func A = [&](real* JU_, const real* AU_) {
		Grid<Vector<real, dim>, gridDim> JU(size, (Vector<real,dim>*)JU_);
		Grid<const Vector<real, dim>, gridDim> AU(size, (const Vector<real,dim>*)AU_);

		//(ln sqrt|g|),a = Gamma^b_ab

		//Gamma^a_bc,d

		//(ln sqrt|g|),ab = Gamma^c_ac,b

		//R_ab = R^c_acb = Gamma^c_ab,c - Gamma^c_ac,b + Gamma^c_dc Gamma^d_ab - Gamma^c_db Gamma^d_ac

		//R^a_b

		//A^a_,b

		//A^a_;b
		
		//A^a;b

		//A^a;b_,b
		
		//A^a;b_;b

		parallel.foreach(range.begin(), range.end(), [&](const Vector<int, gridDim>& index) {
			Vector<real, dim> sum = AU(index) * (-2. * (real)gridDim);
			for (int k = 0; k < gridDim; ++k) {
				Vector<int, gridDim> ip = index;
				ip(k) = std::min<int>(ip(k)+1, n-1);
				Vector<int, gridDim> im = index;
				im(k) = std::max<int>(im(k)-1, 0);
				sum += AU(ip) + AU(im);
			}
			JU(index) = sum / dxSq;
		});
	};

	int volume = size.volume() * dim;
	Solvers::GMRes<real> solver(volume, (real*)AU.v, (const real*)JU.v, A, 1e-7, volume * 10, volume);
	
	solver.stopCallback = [&]()->bool{
		std::cerr << solver.getResidual() << "\t" << solver.getIter() << std::endl;
		return false;
	};

	solver.solve();

	//hmm, TODO if dim < 4 then xyz else txyz
	const char* xs = {"x", "y", "z", "t"};
	assert(numberof(xs) >= gridDim);
	assert(numberof(xs) >= dim);
	std::cout << "#x";
	for (int i = 1; i < gridDim; ++i) {
		std::cout << "\t" << xs[i];
	}
	for (int i = 0; i < dim; ++i) {
		std::cout << "\tA" << xs[i];
	}
	for (int i = 0; i < dim; ++i) {
		std::cout << "\tJ" << xs[i];
	}
	std::cout << std::endl;

	int lastY = 0;
	std::for_each(range.begin(), range.end(), [&](const Vector<int, gridDim>& index) {
		if (index(1) != lastY) {
			lastY = index(1);
			std::cout << std::endl;
		}
		std::cout << index(0);
		for (int i = 1; i < gridDim; ++i) {
			std::cout << "\t" << index(i);
		}
		for (int i = 0; i < dim; ++i) {
			std::cout << "\t" << AU(index)(i);
		}
		for (int i = 0; i < dim; ++i) {
			std::cout << "\t" << JU(index)(i);
		}
		std::cout << std::endl;
	});
}
