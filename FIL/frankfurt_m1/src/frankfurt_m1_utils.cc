#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <cstring>
#include <cstddef>



extern "C" void frankfurt_m1_init_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  size_t gfsize = cctk_ash[0] * cctk_ash[1] * cctk_ash[2] ;
  /*
  std::memset(Enue_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Enua_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Enux_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;

  std::memset(Nnue_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Nnua_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Nnux_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;

  std::memset(Fnue_x_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Fnue_y_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Fnue_z_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;

  std::memset(Fnua_x_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Fnua_y_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Fnua_z_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;

  std::memset(Fnux_x_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Fnux_y_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  std::memset(Fnux_z_rhs, 0, sizeof(CCTK_REAL) * gfsize ) ;
  */

#pragma omp /*parallel for schedule(static)*/ simd
  for(ptrdiff_t ijk=0;ijk<gfsize;++ijk) {

        Nnue_rhs[ijk] = 0. ;
        Nnue_bar_rhs[ijk] = 0. ;
        Nnumu_rhs[ijk] = 0. ;
        Nnumu_bar_rhs[ijk] = 0. ;
        Nnux_rhs[ijk] = 0. ;

        Enue_rhs[ijk] = 0. ;
        Enue_bar_rhs[ijk] = 0. ;
        Enumu_rhs[ijk] = 0. ;
        Enumu_bar_rhs[ijk] = 0. ;
        Enux_rhs[ijk] = 0. ;

        Fnue_x_rhs[ijk] = 0.; Fnue_y_rhs[ijk]=0.; Fnue_z_rhs[ijk]=0.;
        Fnue_bar_x_rhs[ijk] = 0.; Fnue_bar_y_rhs[ijk]=0.; Fnue_bar_z_rhs[ijk]=0.;
        Fnumu_x_rhs[ijk] = 0.; Fnumu_y_rhs[ijk]=0.; Fnumu_z_rhs[ijk]=0.;
        Fnumu_bar_x_rhs[ijk] = 0.; Fnumu_bar_y_rhs[ijk]=0.; Fnumu_bar_z_rhs[ijk]=0.;
        Fnux_x_rhs[ijk] = 0.; Fnux_y_rhs[ijk]=0.; Fnux_z_rhs[ijk]=0.;

   }  // for loop
}


extern "C" void frankfurt_m1_add_rhs(CCTK_ARGUMENTS) 
{
    DECLARE_CCTK_ARGUMENTS; 
    DECLARE_CCTK_PARAMETERS;

    ptrdiff_t gfsize = cctk_ash[0] * cctk_ash[1] * cctk_ash[2] ;

    int const substep = 2 - (*RK_counter) ;

    CCTK_REAL fact = (substep==0) ? 0.5 * CCTK_DELTA_TIME : CCTK_DELTA_TIME ;
    
    #pragma omp simd 
    for(ptrdiff_t i=0; i<gfsize; ++i){
       //if (use_5_spec_m1) {
                 Enue_star[i] = Enue_star_p[i] + fact*Enue_rhs[i] ;
                 Enue_bar_star[i] = Enue_bar_star_p[i] + fact*Enue_bar_rhs[i] ;
                 Enumu_star[i] = Enumu_star_p[i] + fact*Enumu_rhs[i] ;
                 Enumu_bar_star[i] = Enumu_bar_star_p[i] + fact*Enumu_bar_rhs[i] ;
                 Enux_star[i] = Enux_star_p[i] + fact*Enux_rhs[i] ;

                 Nnue_star[i] = Nnue_star_p[i] + fact*Nnue_rhs[i] ;
                 Nnue_bar_star[i] = Nnue_bar_star_p[i] + fact*Nnue_bar_rhs[i] ;
                 Nnumu_star[i] = Nnumu_star_p[i] + fact*Nnumu_rhs[i] ;
                 Nnumu_bar_star[i] = Nnumu_bar_star_p[i] + fact*Nnumu_bar_rhs[i] ;
                 Nnux_star[i] = Nnux_star_p[i] + fact*Nnux_rhs[i] ;

                 Ft_nue_x[i] = Ft_nue_x_p[i] + fact*Fnue_x_rhs[i] ;
                 Ft_nue_y[i] = Ft_nue_y_p[i] + fact*Fnue_y_rhs[i] ;
                 Ft_nue_z[i] = Ft_nue_z_p[i] + fact*Fnue_z_rhs[i] ;

                 Ft_nue_bar_x[i] = Ft_nue_bar_x_p[i] + fact*Fnue_bar_x_rhs[i] ;
                 Ft_nue_bar_y[i] = Ft_nue_bar_y_p[i] + fact*Fnue_bar_y_rhs[i] ;
                 Ft_nue_bar_z[i] = Ft_nue_bar_z_p[i] + fact*Fnue_bar_z_rhs[i] ;

                 Ft_numu_x[i] = Ft_numu_x_p[i] + fact*Fnumu_x_rhs[i] ;
                 Ft_numu_y[i] = Ft_numu_y_p[i] + fact*Fnumu_y_rhs[i] ;
                 Ft_numu_z[i] = Ft_numu_z_p[i] + fact*Fnumu_z_rhs[i] ;

                 Ft_numu_bar_x[i] = Ft_numu_bar_x_p[i] + fact*Fnumu_bar_x_rhs[i] ;
                 Ft_numu_bar_y[i] = Ft_numu_bar_y_p[i] + fact*Fnumu_bar_y_rhs[i] ;
                 Ft_numu_bar_z[i] = Ft_numu_bar_z_p[i] + fact*Fnumu_bar_z_rhs[i] ;

                 Ft_nux_x[i] = Ft_nux_x_p[i] + fact*Fnux_x_rhs[i] ;
                 Ft_nux_y[i] = Ft_nux_y_p[i] + fact*Fnux_y_rhs[i] ;
                 Ft_nux_z[i] = Ft_nux_z_p[i] + fact*Fnux_z_rhs[i] ;

      // } else {
      // 		 Enue_star[i] = Enue_star_p[i] + fact*Enue_rhs[i] ;
      // 		 Enue_bar_star[i] = Enue_bar_star_p[i] + fact*Enue_bar_rhs[i] ;
      // 		 Enux_star[i] = Enux_star_p[i] + fact*Enux_rhs[i] ;

      // 		 Nnue_star[i] = Nnue_star_p[i] + fact*Nnue_rhs[i] ;
      // 		 Nnue_bar_star[i] = Nnue_bar_star_p[i] + fact*Nnue_bar_rhs[i] ;
      // 		 Nnux_star[i] = Nnux_star_p[i] + fact*Nnux_rhs[i] ;

      // 		 Ft_nue_x[i] = Ft_nue_x_p[i] + fact*Fnue_x_rhs[i] ;
      // 		 Ft_nue_y[i] = Ft_nue_y_p[i] + fact*Fnue_y_rhs[i] ;
      // 		 Ft_nue_z[i] = Ft_nue_z_p[i] + fact*Fnue_z_rhs[i] ;

      // 		 Ft_nue_bar_x[i] = Ft_nue_bar_x_p[i] + fact*Fnue_bar_x_rhs[i] ;
      // 		 Ft_nue_bar_y[i] = Ft_nue_bar_y_p[i] + fact*Fnue_bar_y_rhs[i] ;
      // 		 Ft_nue_bar_z[i] = Ft_nue_bar_z_p[i] + fact*Fnue_bar_z_rhs[i] ;

      // 		 Ft_nux_x[i] = Ft_nux_x_p[i] + fact*Fnux_x_rhs[i] ;
      // 		 Ft_nux_y[i] = Ft_nux_y_p[i] + fact*Fnux_y_rhs[i] ;
      // 		 Ft_nux_z[i] = Ft_nux_z_p[i] + fact*Fnux_z_rhs[i] ;

      //} // use_5_spec_m1
    } // for loop

}
