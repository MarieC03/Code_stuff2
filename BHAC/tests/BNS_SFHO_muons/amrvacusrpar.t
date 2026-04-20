!##############################################################################
! include amrvacusrpar - FMtorus

! This file should contain the number of PROBLEM dependent equation parameters,
! the index names for them with values neqpar+1..neqpar+nspecialpar,
! and a string giving the names for the file header. For example:
!
! INTEGER,PARAMETER:: mass_=neqpar+1, nspecialpar=1
! CHARACTER*4,PARAMETER:: specialparname='mass'
!
! By default there are no special parameters
  INTEGER,PARAMETER:: &
       kappa_=neqpar+1, &
       rin_=neqpar+2, &
       rtmax_=neqpar+3, &
       rout_=neqpar+4, &
       rt1_=neqpar+5, &
       q1_=neqpar+6, &
       q2_=neqpar+7, &
       q3_=neqpar+8, &              
       beta_=neqpar+9, &
       pmin_=neqpar+10, &
       rhomin_=neqpar+11, &
       Nloops_=neqpar+12, &
       fliploops_=neqpar+13, &
       lambdar_=neqpar+14
INTEGER,PARAMETER:: nspecialpar=14
CHARACTER*1,PARAMETER:: specialparname=' '
CHARACTER*20,PARAMETER:: typeuser='FMtorus'
! end include amrvacusrpar - WJtorus
!##############################################################################
