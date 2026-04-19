module mod_m1_eas_test

  use mod_m1_eas
  
  implicit none

  integer, parameter :: test_eas_zero  = 1
  integer, parameter :: test_eas_scatt = 2


  
contains
  subroutine m1_eas_test_activate(whichtest)    
    integer, intent(in) :: whichtest

    if( whichtest .eq. test_eas_zero ) then

       m1_get_eas => m1_get_eas_test_zero

    else if ( whichtest .eq. test_eas_scatt ) then

       m1_get_eas => m1_get_eas_test_scatt

    else
      {#IFNDEF UNIT_TESTS
       call mpistop("Error, invalid test case requested for M1 eas.")
       }
    end if
  end subroutine m1_eas_test_activate

  subroutine m1_get_eas_test_zero(wrad,speciesKSP,ix^D,eas_eq,fluid_Prim) !T_fluid,rho_fluid,ye_fluid)
    use mod_m1_internal
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    {#IFNDEF UNIT_TESTS
     include "amrvacdef.f"
     }
     integer, intent(in) :: ix^D, speciesKSP
     double precision, intent(in)    :: wrad(1:m1_numvars_internal)
     double precision, intent(inout) :: fluid_Prim(1:fluid_vars)
     double precision, intent(inout) :: eas_eq(1:m1_num_eas)

     eas_eq(:) = 0.0d0
   end subroutine m1_get_eas_test_zero

   subroutine m1_get_eas_test_scatt(wrad,speciesKSP,ix^D,eas_eq,fluid_Prim) !T_fluid,rho_fluid,ye_fluid)
    use mod_m1_internal
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    {#IFNDEF UNIT_TESTS
     include "amrvacdef.f"
     }
     integer, intent(in) :: ix^D, speciesKSP
     double precision, intent(in)    :: wrad(1:m1_numvars_internal)
     double precision, intent(inout) :: fluid_Prim(1:fluid_vars)
     double precision, intent(inout) :: eas_eq(1:m1_num_eas)
     
     eas_eq(k_s) = 1.0d5
     eas_eq(Q_ems) = 0.0d0
     eas_eq(Q_ems_n) = 0.0d0
     eas_eq(k_a) = 0.0d0
     eas_eq(k_n) = 0.0d0
    
  end subroutine m1_get_eas_test_scatt

  !----------------------------------------------------
  !subroutine m1_get_eas_test_zero(wprims, ixI^L,ixO^L)
  !  include "amrvacdef.f"
  !  integer, intent(in) :: ixI^L, ixO^L
  !  double precision, intent(in) :: wprims(ixI^S,1:nw)
  !  m1_eas(ixO^S,:,:) = zero
  !end subroutine m1_get_eas_test_zero

  !subroutine m1_get_eas_test_scatt(wprims,ixI^L,ixO^L)
  !  include "amrvacdef.f"
  !  integer, intent(in) :: ixI^L, ixO^L
  !  double precision, intent(in) :: wprims(ixI^S,1:nw)
  !  m1_eas(ixO^S,k_s,:) = 1.0d5
  !  m1_eas(ixO^S,eta,:) = zero
  !  m1_eas(ixO^S,eta_n,:) = zero
  !  m1_eas(ixO^S,k_a,:) = zero
  !  m1_eas(ixO^S,k_n,:) = zero  
  !end subroutine m1_get_eas_test_scatt

end module mod_m1_eas_test
