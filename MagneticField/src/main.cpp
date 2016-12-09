#include <iostream>
#include "Parallel/Parallel.h"
#include "Tensor/Grid.h"
#include "Solvers/GMRes.h"
#include "Common/Macros.h"

Parallel::Parallel parallel(8);

const int gridDim = 3;
const int dim = 4;

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
	real dx = 1;
	real dxSq = dx*dx;

	//solve for AU
	RangeObj<gridDim> range = JU.range();
#if 0
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
	});
#endif
	
	//lazy rasterization
	for (int i = 0; i < 8 * (int)n; ++i) {
		real th = 2. * M_PI * (real)i / (real)(4*n);
		int x = .5 * size(0) + .25 * n * cos(th);
		int y = .5 * size(1) + .25 * n * sin(th);
		int z = .5 * size(2);
		JU(x,y,z)(0) = 1;
		JU(x,y,z)(1) = -sin(th);
		JU(x,y,z)(2) = cos(th);
		JU(x,y,z)(3) = 0;
	}

	parallel.foreach(range.begin(), range.end(), [&](const Vector<int, gridDim>& index) {
		AU(index) = -JU(index);
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
	Solvers::GMRes<real> solver(volume, (real*)AU.v, (const real*)JU.v, A, 2e-2, volume, 100);
	
	solver.stopCallback = [&]()->bool{
		std::cerr << solver.getResidual() << "\t" << solver.getIter() << std::endl;
		return false;
	};

	solver.solve();


	//calculate E and B
	Grid<Vector<real, dim-1>, gridDim> E(size), B(size);

	//E^i = -grad phi - d/dt A
	//B^i = curl A
	parallel.foreach(range.begin(), range.end(), [&](const Vector<int, gridDim>& index) {
		Vector<int,dim> dAU[3];
		for (int i = 0; i < gridDim; ++i) {
			Vector<int, gridDim> ip = index;
			ip(i) = std::min<int>(ip(i)+1, n-1);
			Vector<int, gridDim> im = index;
			im(i) = std::max<int>(im(i)-1, 0);
			for (int u = 0; u < gridDim; ++u) {
				dAU[i](u) = (AU(ip)(u) - AU(im)(u)) / (2. * dx);	//... TODO - d/dt AU
			}
		}
		for (int i = 0; i < 3; ++i) {
			int i2 = (i+1)%3;
			int i3 = (i+2)%3;
			E(index)(i) = -dAU[i](0);	//... TODO - d/dt AU
			B(index)(i) = dAU[i2](i3+1) - dAU[i3](i2+1);
		}
	});

	//hmm, TODO if dim < 4 then xyz else txyz
	const char* xs[] = {"t", "x", "y", "z"};
	assert(numberof(xs) >= gridDim);
	assert(numberof(xs) >= dim);
	
	struct Col {
		std::string name;
		std::function<real(Vector<int,gridDim>)> func;
		Col(std::string name_, std::function<real(Vector<int,gridDim>)> func_) : name(name_), func(func_) {}
	};
	std::vector<Col> cols = {
		{"ix", [&](Vector<int,gridDim> index)->real{ return index(0); }},
		{"iy", [&](Vector<int,gridDim> index)->real{ return index(1); }},
		{"iz", [&](Vector<int,gridDim> index)->real{ return index(2); }},
	};
	for (int i = 0; i < dim; ++i) {
		cols.push_back(Col(
			std::string("A^") + std::string(xs[i]),
			[&,i](Vector<int,gridDim> index)->real{ return AU(index)(i); }
		));
	}
	for (int i = 0; i < dim; ++i) {
		cols.push_back(Col(
			std::string("J^") + std::string(xs[i]),
			[&,i](Vector<int,gridDim> index)->real{ return JU(index)(i); }
		));
	}
	for (int i = 0; i < 3; ++i) {
		cols.push_back(Col(
			std::string("E^") + std::string(xs[i+1]),
			[&,i](Vector<int,gridDim> index)->real{ return E(index)(i); }
		));
		cols.push_back(Col(
			std::string("B^") + std::string(xs[i+1]),
			[&,i](Vector<int,gridDim> index)->real{ return B(index)(i); }
		));
	}

	std::string outputFilename = "out.txt";
	FILE* file = fopen(outputFilename.c_str(), "w");
	if (!file) throw Common::Exception() << "failed to open file " << outputFilename;

	fprintf(file, "#");
	{
		const char* tab = "";
		for (std::vector<Col>::iterator p = cols.begin(); p != cols.end(); ++p) {
			fprintf(file, "%s%s", tab, p->name.c_str());
			tab = "\t";
		}
	}
	fprintf(file, "\n");
	//this is printing output, so don't do it in parallel		
	for (RangeObj<gridDim>::iterator iter = range.begin(); iter != range.end(); ++iter) {
		const char* tab = "";
		for (std::vector<Col>::iterator p = cols.begin(); p != cols.end(); ++p) {
			fprintf(file, "%s%.16e", tab, p->func(iter.index));
			tab = "\t";
		}
		fprintf(file, "\n");
	}

	fclose(file);
}
