Subroutine dfunc(x,g)
! version 3.0     July 18, 2019
  implicit none

  real*8                          :: x(6),g(6)

  call CalcGrad(x,g)
  return

End Subroutine dfunc
