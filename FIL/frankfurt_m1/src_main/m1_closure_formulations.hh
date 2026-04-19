#ifndef FRANKFURT_M1_CLOSURE_F_HH
#define FRANKFURT_M1_CLOSURE_F_HH
/**
 * @file m1_closure.hh
 * @author Carlo Musolino
 *
 * M1 closure implementations
 * 
 *
 */

namespace frankfurt_m1 {


  class minerbo_closure_f {
  public:

    static double closure_func(double const zeta) {
      return 1./3. + zeta*zeta*(6.-2.*zeta+6.*zeta*zeta)/15. ;
    };

    static double closure_func(double const zeta, double& dchi_dzeta) {
      dchi_dzeta = 2.*zeta*(6.-2.*zeta+6.*zeta*zeta)/15. + zeta*zeta*(12.*zeta-2.)/15.;
      return 1./3. + zeta*zeta*(6.-2.*zeta+6.*zeta*zeta)/15. ;
    };
  };

  class eddington_closure_f {
  public:
    static double closure_func(double const zeta) {
      return 1/3;
    };
    static double closure_func(double const zeta, double& dchi_dzeta) {
      dchi_dzeta = 0.;
      return 1./3.;
    };
  };

  class levermore_closure_f {
  public:
    static double closure_func(double const zeta) {
      return ( 3. + 4.*zeta*zeta) / ( 5. + 2.*sqrt(4.-3.*zeta*zeta) ) ;
    }

    static double closure_func(double const zeta, double& dchi_dzeta) {
      dchi_dzeta = 2.*zeta / sqrt(4.-3.*zeta*zeta) ;
      return ( 3. + 4.*zeta*zeta) / ( 5. + 2.*sqrt(4.-3.*zeta*zeta) ) ;
    }
  } ;
 

  // The following is an analytic fit for the closure
  // in the radiating sphere test, see Murchikova, E. M.+
  // for details
  class analytic_closure_f {
  public:
    static double closure_func(double const zeta) {
      return (1./3. - 1./3.*zeta + 2./3.*zeta*zeta) * (zeta<=0.5)
	+ (1./3. - 2./3.*zeta + 4./3.*zeta*zeta) * (zeta>0.5) ;
    };

    static double closure_func(double const zeta, double& dchi_dzeta) {
      dchi_dzeta = ( -1./3. + 4./3.*zeta) * (zeta<=0.5) +
	( -2./3. + 8./3.*zeta) * (zeta>0.5) ;
      return (1./3. - 1./3.*zeta + 2./3.*zeta*zeta) * (zeta<=0.5)
	+ (1./3. - 2./3.*zeta + 4./3.*zeta*zeta) * (zeta>0.5) ;
    };
    
  } ;

  
}

#endif 
