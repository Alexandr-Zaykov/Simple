Function func(x)
! version 3.0     July 18, 2019
  implicit none

  real*8                          :: func,SF
  integer,parameter               :: NMAX=6
  real*8                          :: x(NMAX)

  call MatElm(x,SF)
  func=SF

End Function func
