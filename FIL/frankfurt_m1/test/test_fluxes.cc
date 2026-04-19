#define STANDALONE
#include <iostream>
#include "../src/m1_closure.hh"
#include "../src/frankfurt_m1.h"
//#include "../src/legacy/M1closure.C" 
//#include "../src/legacy/FIL_M1_headers.h"

constexpr int flux_dirn = 1;

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
  

  std::array<double,NUM_EAS> eas{0.1,100,0.,0.1,100} ;

  double E = 1;
  double Fx = 0.;
  double Fy = 0;
  double Fz = 0;
  double velx,vely,velz;

  compute_F(1.,0.,0.,0.*E,&Fx,&Fy,&Fz) ;
  compute_F(1.,1.,1.,0.0,&velx,&vely,&velz) ;

  double zold = 0.5 ;
  closure_t cl( { velx,vely,velz}, &gamma) ;
  cl.update_closure(E, {Fx,Fy,Fz}, zold,true) ;
  
  double cp,cm ;

  cl.get_wavespeeds<flux_dirn>(cp,cm) ;
  double N{E}, Ef, Fxf, Fyf, Fzf, Nf ;
  cl.get_fluxes<flux_dirn>(N,Nf,Ef,Fxf,Fyf,Fzf) ;


  std::cout << "Wavespeeds: " << cp << ", " << cm << "\n\n";


  std::cout << "Vars:   " << E  << ", " << Fx  << ", " << Fy  << ", " << Fz  << "\n\n" ;
  std::cout << "Fluxes: " << Ef << ", " << Fxf << ", " << Fyf << ", " << Fzf << "\n\n" ;
  

  
}
