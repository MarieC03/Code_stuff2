module mod_m1_minerbo_closure
  use mod_m1_closure

contains
  
  subroutine m1_minerbo_closure_activate()
     m1_closure_func    => minerbo_closure
     m1_closure_deriv   => minerbo_deriv
   end subroutine m1_minerbo_closure_activate
  
  function minerbo_closure(z)
    double precision :: minerbo_closure
    double precision, intent(in) :: z
    minerbo_closure = 1.0d0/3.0d0 + z**2*( 6.0d0 - 2.0d0*z + 6.0d0*z**2 ) / 15.0d0
  end function minerbo_closure

  function minerbo_deriv(z)
    double precision :: minerbo_deriv
    double precision, intent(in) :: z
    minerbo_deriv = 2.0d0*z*(6.0d0-2.0d0*z+6.0d0*z**2)/15.0d0 + z**2*(12.0d0*z -2.0d0)/15.0d0
  end function minerbo_deriv
  
end module mod_m1_minerbo_closure

