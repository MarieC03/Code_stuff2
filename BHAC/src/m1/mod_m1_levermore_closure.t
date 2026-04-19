module mod_m1_levermore_closure
  !use mod_m1_closure

contains

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
  
end module mod_m1_levermore_closure

