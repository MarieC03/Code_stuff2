#define STANDALONE
#include <iostream>
#include "../src/m1_closure.hh"
#include "../src/frankfurt_m1.h"
//#include "../src/legacy/M1closure.C" 
//#include "../src/legacy/FIL_M1_headers.h"


void compute_F(double nx, double ny, double nz, double F, double*Fx, double*Fy, double* Fz) {
  double const n = std::sqrt(nx*nx + ny*ny + nz*nz );

  nx/=n; ny/=n; nz/=n;

  *Fx = F * nx ;
  *Fy = F * ny ;
  *Fz = F * nz ;
}

int main() {
  using namespace frankfurt_m1 ;
  
  metric_c gamma( {1.,0.,0.,1.,0.,1.}, {0.,0.,0.}, 1.) ;
  

  //std::array<double,NUM_EAS> eas{1.e-02,1.e-01,0.,0.1,100} ;
  std::array<double,NUM_EAS> eas{1.58456032e-05,362125.31020406,98934.68212067,0.1,100} ;
  
  double E = 1.36412849e-11;
  double Fx = -2.23562081e-12;
  double Fy = 1.04681314e-11;
  double Fz = 6.14514753e-21;
  double velx,vely,velz;

  //compute_F(1.,0.,0.,0.9*E,&Fx,&Fy,&Fz) ;
  compute_F(1.,1.,0.,0.2,&velx,&vely,&velz) ;
  
  double dE, dSx, dSy, dSz, dYe ;
  double N = 1.;
  double zold = 0. ;
  bool compute_zeta = true;


  closure_t cl( { velx,vely,velz}, &gamma) ;
  //auto z = cl.update_closure(E, {Fx,Fy,Fz}, zold,true) ;

  //std::cout << "z= " << z << std::endl ;
  
  m1_implicit_conventions::m1_implicit_error_bits error_bits ;
  auto z = cl.get_collisional_rootfinder( eas,
				     N, E,
				     Fx,Fy,Fz,
				     dE, dSx, dSy, dSz,
				     dYe,
				     zold,
				     1.e-04 ,
				     true,
				     true,
				     error_bits ) ;

  std::cout << z << std::endl ;
  
}
