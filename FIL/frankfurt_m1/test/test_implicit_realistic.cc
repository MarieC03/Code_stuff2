#define STANDALONE
#include <iostream>
#include "../src/m1_closure.hh"
#include "../src/frankfurt_m1.h"
//#include "../src/legacy/M1closure.C" 
//#include "../src/legacy/FIL_M1_headers.h"
#include <hdf5.h>
#include "test_implicit_rootfinder.h"

void compute_F(double nx, double ny, double nz, double F, double*Fx, double*Fy, double* Fz) {
  double const n = std::sqrt(nx*nx + ny*ny + nz*nz );

  nx/=n; ny/=n; nz/=n;

  *Fx = F * nx ;
  *Fy = F * ny ;
  *Fz = F * nz ;
}

int main() {
  using namespace frankfurt_m1 ;
  std::string fname{ "simdata.h5" };
  metric_c gamma( {readvar("gxx"),
		   readvar("gxy"),
		   readvar("gxz"),
		   readvar("gyy"),
		   readvar("gyz"),
		   readvar("gzz")}, {readvar("betax"),readvar("betay"),readvar("betaz")} ,readvar("alp")) ;
  

  //std::array<double,NUM_EAS> eas{1.e-02,1.e-01,0.,0.1,100} ;
  std::array<double,NUM_EAS> eas{readvar("Qnua"),readvar("kappa_nua_a"),readvar("kappa_nua_s"),0.1,100} ;
  
  double E = readvar("Enua");
  double Fx = readvar("Fnua_x");
  double Fy = readvar("Fnua_y");
  double Fz = readvar("Fnua_z");
  double velx = readvar("vel[0]");
  double vely = readvar("vel[1]");
  double velz = readvar("vel[2]") ;
  
  double dE, dSx, dSy, dSz, dYe ;
  double N = 1.;
  double zold = 0. ;
  bool compute_zeta = true;


  closure_t cl( { velx,vely,velz}, &gamma) ;
  //auto z = cl.update_closure(E, {Fx,Fy,Fz}, zold,true) ;

  //std::cout << "z= " << z << std::endl ;

  std::cout << "Initial state: "
	    << E << '\n'
	    << Fx << '\n'
	    << Fy << '\n'
	    << Fz << '\n' ; 
  
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


static double readvar(std::string const& vname, std::string const& fname){
  hid_t file_id, dspace_id, dset_id, group_id ;
  herr_t h5err ;
  file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT ) ;

  dset_id = H5Dopen( file_id, vname.c_str(), H5F_ACC_RDONLY ) ;
  double value ;
  h5err = H5Dread( dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value ) ;
  return value ;
}
