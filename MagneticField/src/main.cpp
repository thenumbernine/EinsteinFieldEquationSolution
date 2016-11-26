#include <iostream>
#include "Parallel/Parallel.h"
#include "Tensor/Grid.h"
#include "Solvers/ConjGrad.h"
Parallel::Parallel parallel(8);
enum { dim = 3 };
typedef double real;
using namespace Tensor;
int main() {
	int n = 32;
	Vector<int,dim> size(n,n,n);
	Grid<real, dim> rho(size);
	real M = 1e+6;
	parallel.foreach(rho.range().begin(), rho.range().end(), [&](const Vector<int,dim>& index) {
		rho(index) = (std::abs(index(0) - size(0)/2) < 1 &&
			std::abs(index(1) - size(1)/2) < 1 &&
			std::abs(index(2) - size(2)/2) < 1)
			? (index(0) >= size(0)/2 ? M : -M) : 0;
	});
	real dxSq = 1/(real)size(0);
	Grid<real,dim> phi = rho;
	Solvers::ConjGrad<real> solver(
		size.volume(),	//n
		phi.v,			//x
		rho.v,			//b
		[&](real* result_, const real* phi_) {
			//hmm, should I provide a function to do this...
			Grid<real,dim> result; result.v=result_; result.size=size;	
			Grid<const real,dim> phi; phi.v=phi_; phi.size=size;
			parallel.foreach(phi.range().begin(), phi.range().end(), [&](const Vector<int,dim>& index) {
				//build the source
				if (index(0) == 0 || index(0) == size(0)-1 ||
					index(1) == 0 || index(1) == size(1)-1 ||
					index(2) == 0 || index(2) == size(2)-1)
				{
					result(index) = 0;
				} else {
					result(index) = (phi(index + Vector<int,dim>(1,0,0))
								+ phi(index - Vector<int,dim>(1,0,0))
								+ phi(index + Vector<int,dim>(0,1,0))
								+ phi(index - Vector<int,dim>(0,1,0))
								+ phi(index + Vector<int,dim>(0,0,1))
								+ phi(index - Vector<int,dim>(0,0,1))
								- 6. * phi(index))
								/ (dxSq * 6. * M_PI);
				}
			});
			phi.v = NULL;	//to skip the dtor.  very bad hack, I know.
		},		//A
		1e-50,	//epsilon
		0	//maxiter
	);
	solver.stopCallback = [&]()->bool{
		std::cout << solver.getResidual() << "\t" << solver.getIter() << std::endl;
		return false;
	};
	solver.solve();

	std::cout << "#x y z rho phi" << std::endl;
	std::for_each(phi.range().begin(), phi.range().end(), [&](const Vector<int,dim>& index) {
		std::cout << index(0) 
			<< "\t" << index(1)
			<< "\t" << index(2)
			<< "\t" << rho(index) 
			<< "\t" << phi(index)
			<< std::endl;
	});
}
