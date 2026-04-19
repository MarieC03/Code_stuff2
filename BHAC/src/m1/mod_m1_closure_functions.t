module mod_m1_closure_functions

contains
  
  !> Minerbo closure ------------------
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

  !> Levermore closure ------------------
  function levermore_closure(z)
    double precision :: levermore_closure
    double precision, intent(in) :: z
    levermore_closure = ( 3.0d0 + 4.d0*z*z) / ( 5.0d0 + 2.0d0*dsqrt(4.0d0-3.0d0*z*z) )
  end function levermore_closure

  function levermore_deriv(z)
    double precision :: levermore_deriv
    double precision, intent(in) :: z
    levermore_deriv =  2.0d0*z / dsqrt(4.0d0-3.0d0*z*z) 
  end function levermore_deriv
  
end module mod_m1_closure_functions

