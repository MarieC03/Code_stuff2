#define STANDALONE
#include <iostream>
#include "../src/m1_closure.hh"

//#include "../src/legacy/M1closure.C" 
//#include "../src/legacy/FIL_M1_headers.h"

constexpr CCTK_REAL NOT_USED = 0.0;
constexpr CCTK_REAL gxx = 1.;
constexpr CCTK_REAL gyy = 1.;
constexpr CCTK_REAL gzz = 1.;
constexpr CCTK_REAL gxy = 0.0;
constexpr CCTK_REAL gxz = 0.0;
constexpr CCTK_REAL gyz = 0.0;
constexpr CCTK_REAL gupxx = 1./gxx;
constexpr CCTK_REAL gupyy = 1./gyy;
constexpr CCTK_REAL gupzz = 1./gzz;
constexpr CCTK_REAL gupxy = 0.0;
constexpr CCTK_REAL gupxz = 0.0;
constexpr CCTK_REAL gupyz = 0.0;
constexpr CCTK_REAL lapse = 1.;
constexpr CCTK_REAL psi4 = 1.0;
constexpr CCTK_REAL betax = 0.0;
constexpr CCTK_REAL betay = 0.0;
constexpr CCTK_REAL betaz = 0.0;

constexpr int SPECIES=0;




int main() {

  metric_c gamma( {1.,0.,0.,1.,0.,1.}, {0.,0.,0.}, 1.) ;
  closure_t cl( { 0.1,0.,0.}, &gamma) ;
  double z = cl.update_closure( 1., {1.,0.,0.}, 0.5, true) ;

  std::cout << z << std::endl ;

  
}
