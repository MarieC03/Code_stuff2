#include <iostream>
#include "constrained_nr.hh"

int main() {

  auto const f = [] ( const double x ) {
		   return x*x*x - x - 2. ;
		 };
  auto const df = [] (const double x) {
		    return 3. * x*x - 1. ;
		  };
  int err ;
  std::cout << zero_constrained_newton_raphson(1., 2., 1e-03, f, df, 100, err) << std::endl ;

}
