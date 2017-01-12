#include <iostream>
#include "Parallel/Parallel.h"
#include "Tensor/Grid.h"
#include "Solvers/ConjGrad.h"
#include "Solvers/ConjRes.h"
#include "Solvers/GMRes.h"
#include "Common/Macros.h"

Parallel::Parallel parallel(8);

const int gridDim = 3;
const int stDim = 4;

typedef double real;
using namespace Tensor;

/*
charge density
coulumbs per meter^3
1 coulumb = 6.242e+18 e
how do you convert coulumbs to meters?
by setting the Coulumb constant to 1?  via Planck units?
1 = ~c~ = 299792458 m / s = c m / s
	m / s = 1 / c
	1 s = c m
	1 s = 299792458 m
1 = ~hBar~ = 1.05457180013e-34 kg m^2 / s 
	1 = hBar kg m^2 / s
	kg = c / hBar / m
	kg = 2.842788494468e+42 / m
1 = ~G~ = 6.67408e-11 m^3 / (kg s^2) 
	kg = G m^3 / s^2
	kg = G / c^2 m
	kg = 7.4259154861063e-28 m
	...but kg = c / hBar / m
	...so c / hBar / m = G / c^2 m
		m = sqrt(c^3 / (hBar G))
		m = 6.1872444306747e+34
		s = c m
		s = 1.8548892161188e+43
		kg = 45945954.234073
1 = ~ke~ = ke kg m^3 / (s^2 C^2) = 8.9875517873681764e+9 kg m^3 / (s^2 C^2) 
	C^2 = ke kg m^3 / s^2
	C^2 = ke G / c^4 m^2
	C = sqrt(ke G) / c^2 m using G's def
	C = 8.6173751723517e-18 m
	C^2 = ke / (hBar c) using hBar's def
	C = = 5.3317806542168e+17 using either
1 = ~kB~ = kB m^2 kg / (K s^2) = 1.3806488e-23 m^2 kg / (K s^2)
	K = kB kg m^2 / s^2 = kB / c^2 kg 
	K = kB G / c^4 m using G's def
	K = kB / (c hBar m) using hBar's def
	K = 7.0581208407927e-33

1 C = e ~e~ = 6.2415093414e+18 ~e~
	~e~ = C / e = sqrt(ke G) / (c^2 e) m
	~e~ = 1.380655655707e-36 m
	~e~ = 1.380655655707e-36 m
	~e~ = 0.085424540164524


J = kg m^2 / s^2 = kg / c^2 = G / c^4
1 erg = 1e-7 J = 1e-7 kg m^2 / s^2 = 1e-7 G / c^4
J = C V 

*/
real c = 299792458;
real G = 6.67408e-11;
real ke = 8.9875517873681764e+9;	// = 1 / (4 pi epsilon0)
real hBar = 1.05457180013e-34;
real kB = 1.3806488e-23;
real e = 6.2415093414e+18;

/*

J_a = (rho, j_i)

flat space:

F^uv_,v = 4 pi J^u
F_uv = A_v,u - A_u,v
A^v,u_v - A^u,v_v = 4 pi J^u

curved space:

A^a;u_u - A^u_;u^;a - R^a_u A^u = 4 pi J^a

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

void time(const std::string name, std::function<void()> f) {
	std::cout << name << " ... ";
	std::cout.flush();
	auto start = std::chrono::high_resolution_clock::now();
	f();
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end - start;
	std::cout << "(" << diff.count() << "s)" << std::endl;
}

int main() {
	size_t n = 64;
	Vector<int, gridDim> size;
	for (int i = 0; i < gridDim; ++i) {
		size(i) = n;
	}
	Vector<real, gridDim> xmin(-1,-1,-1);
	Vector<real, gridDim> xmax(1,1,1);
	Vector<real, gridDim> dx = (xmax - xmin) / (Vector<real,gridDim>)size;
	Grid<Vector<real, stDim>, gridDim> JU(size);
	Grid<Vector<real, stDim>, gridDim> AU(size);

	//solve for AU
	RangeObj<gridDim> range = JU.range();
#if 0
	parallel.foreach(range.begin(), range.end(), [&](const Vector<int, gridDim>& index) {
		Vector<real,gridDim> x = ((Vector<real,gridDim>)index + .5) - size/2;
		real linfDist = 0;
		real sign = 1;
		for (int i = 0; i < gridDim; ++i) {
			linfDist = std::max<real>(linfDist, fabs(x(i)));
			if (x(i) < 0) sign = -sign;
		}
		const real r = 15;
		JU(index)(0) = linfDist < r ? sign : 0;
	});
#endif
	
	//lazy rasterization
	real q = 1;	//Coulombs (C)
	q *= sqrt(ke * G) / (c * c);	//...times Columbs to meters (m)
	q *= dx.volume();	// ... per meter cubed (1/m^2)
//q = 1e-22 for n=64 (dx=.03125)
//that might be too small for the solvers...
q = 1;
std::cout << "q: " << q << std::endl;
	
	int divs = 8 * (int)n;
	for (int i = 0; i < divs; ++i) {
		real frac = (real)i / (real)divs;
		real th = 2. * M_PI * frac;
		int x = .5 * size(0) + .25 * n * cos(th);
		int y = .5 * size(1) + .25 * n * sin(th);
		int z = .5 * size(2);
		
		JU(x,y,z)(0) = q;
		
		real I = q;	//amps (C / s) => 1/(m^2 s)
		I *= c;		//(1/m^3)
		
		JU(x,y,z)(1) = -I * sin(th);
		JU(x,y,z)(2) = I * cos(th);
		JU(x,y,z)(3) = 0;
	}

	parallel.foreach(range.begin(), range.end(), [&](const Vector<int, gridDim>& index) {
		AU(index) = -JU(index);
	});

	Solvers::Krylov<real>::Func A = [&](real* JU_, const real* AU_) {
		Grid<Vector<real, stDim>, gridDim> JU(size, (Vector<real,stDim>*)JU_);
		Grid<const Vector<real, stDim>, gridDim> AU(size, (const Vector<real,stDim>*)AU_);

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

		//solve for del A^i = J^i
		parallel.foreach(range.begin(), range.end(), [&](const Vector<int, gridDim>& index) {
#if 0	//set boundary to zero?
			for (int k = 0; k < gridDim; ++k) {
				if (index(k) == 0 || index(k) == (int)n-1) {
					//set boundary to zero:
					JU(index) = Vector<real, stDim>();
					return;
				}
			}
#endif
			
			Vector<real, stDim> sum = AU(index) * (-2. * (
				1. / (dx(0) * dx(0))
				+ 1. / (dx(1) * dx(1))
				+ 1. / (dx(2) * dx(2))
			));
			for (int k = 0; k < gridDim; ++k) {
				Vector<int, gridDim> ip = index;
				ip(k) = std::min<int>(ip(k)+1, n-1);
				Vector<int, gridDim> im = index;
				im(k) = std::max<int>(im(k)-1, 0);
				sum += (AU(ip) + AU(im)) / (dx(k) * dx(k));
			}
			JU(index) = sum;
		});
	};

	int volume = size.volume() * stDim;

//hmm, some problems:
//ConjGrad is segfaulting
//ConjRes and GMRes stop immediately if q=1e-22
//both of them bottom-out at .0002 if q=1
//...and the Lua implementations of these works fine, but slower...

	//Solvers::ConjGrad<real> solver(volume, (real*)AU.v, (const real*)JU.v, A, 1e-5, volume);
	
	//ConjRes took 0.0491484s to solve within 1e-7
	Solvers::ConjRes<real> solver(volume, (real*)AU.v, (const real*)JU.v, A, 1e-5, volume);

	//GMRes took 0.507793s to solve within 1e-7 with a restart of 100
	//Solvers::GMRes<real> solver(volume, (real*)AU.v, (const real*)JU.v, A, 1e-5, volume, 100);
	
	solver.stopCallback = [&]()->bool{
		std::cerr << solver.getResidual() << "\t" << solver.getIter() << std::endl;
		return false;
	};

	time("solve", [&](){
		solver.solve();
	});

	//calculate E and B
	Grid<Vector<real, 3>, gridDim> E(size), B(size);

	//E^i = -grad phi - d/dt A
	//B^i = curl A
	parallel.foreach(range.begin(), range.end(), [&](const Vector<int, gridDim>& index) {
		real dAU[gridDim][stDim];
		for (int i = 0; i < gridDim; ++i) {
			Vector<int, gridDim> ip = index;
			ip(i) = std::min<int>(ip(i)+1, n-1);
			Vector<int, gridDim> im = index;
			im(i) = std::max<int>(im(i)-1, 0);
			for (int u = 0; u < stDim; ++u) {
				dAU[i][u] = (AU(ip)(u) - AU(im)(u)) / (2. * dx(i));	//... TODO - d/dt AU
			}
		}
		for (int i = 0; i < 3; ++i) {
			int i2 = (i+1)%3;
			int i3 = (i+2)%3;
			E(index)(i) = -dAU[i][0];	//... TODO - d/dt AU
			B(index)(i) = dAU[i2][i3+1] - dAU[i3][i2+1];
		}
	});

	//hmm, TODO if stDim < 4 then xyz else txyz
	const char* xs[] = {"t", "x", "y", "z"};
	assert(numberof(xs) >= gridDim);
	assert(numberof(xs) >= stDim);
	
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
	for (int i = 0; i < stDim; ++i) {
		cols.push_back(Col(
			std::string("A^") + std::string(xs[i]),
			[&,i](Vector<int,gridDim> index)->real{ return AU(index)(i); }
		));
	}
	for (int i = 0; i < stDim; ++i) {
		cols.push_back(Col(
			std::string("J^") + std::string(xs[i]),
			[&,i](Vector<int,gridDim> index)->real{ return JU(index)(i); }
		));
	}
	
	cols.push_back(Col("E^x", [&](Vector<int,gridDim> index)->real{ return E(index)(0); }));
	cols.push_back(Col("E^y", [&](Vector<int,gridDim> index)->real{ return E(index)(1); }));
	cols.push_back(Col("E^z", [&](Vector<int,gridDim> index)->real{ return E(index)(2); }));
	cols.push_back(Col("B^x", [&](Vector<int,gridDim> index)->real{ return B(index)(0); }));
	cols.push_back(Col("B^y", [&](Vector<int,gridDim> index)->real{ return B(index)(1); }));
	cols.push_back(Col("B^z", [&](Vector<int,gridDim> index)->real{ return B(index)(2); }));
	
	cols.push_back(Col("|E|", [&](Vector<int,gridDim> index)->real{
		return Vector<real,gridDim>::length(E(index));
	}));
	cols.push_back(Col("|B|", [&](Vector<int,gridDim> index)->real{
		return Vector<real,gridDim>::length(B(index));
	}));
	cols.push_back(Col("div_E", [&](Vector<int,gridDim> index)->real{
		real div = 0.;
		for (int i = 0; i < gridDim; ++i) {
			Vector<int, gridDim> ip = index;
			ip(i) = std::min<int>(ip(i)+1, n-1);
			Vector<int, gridDim> im = index;
			im(i) = std::max<int>(im(i)-1, 0);
			div += (E(ip)(i) - E(im)(i)) / (2. * dx(i));	//... TODO - d/dt AU
		}
		return div;
	}));
	cols.push_back(Col("div_B", [&](Vector<int,gridDim> index)->real{
		real div = 0.;
		for (int i = 0; i < gridDim; ++i) {
			Vector<int, gridDim> ip = index;
			ip(i) = std::min<int>(ip(i)+1, n-1);
			Vector<int, gridDim> im = index;
			im(i) = std::max<int>(im(i)-1, 0);
			div += (B(ip)(i) - B(im)(i)) / (2. * dx(i));	//... TODO - d/dt AU
		}
		return div;
	}));


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
