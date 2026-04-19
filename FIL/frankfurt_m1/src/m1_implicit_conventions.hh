
#ifndef _M1_IMPLICIT_
#define _M1_IMPLICIT_

#include<stddef.h>
#include<bitset>

namespace m1_implicit_conventions {

  static constexpr size_t MAXITER_NR_CLOSURE = 10 ;
  static constexpr size_t MAXITER_NR_SOURCES = 20 ;
  static constexpr double CLOSURE_TOL = 1.e-15 ;
  static constexpr double SOURCES_TOL = 1.e-15 ;
  static constexpr double SOURCES_TOL_BACKUP = 1.e-10 ;
  
  enum implicit_sources_errors {
    ERRNOCONV=0,
    FLUXFIX,
    EFIX,
    NUM_ERRORS } ;

  
  
  typedef std::bitset<implicit_sources_errors::NUM_ERRORS> m1_implicit_error_bits ;
}
#endif 
