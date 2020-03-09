Subroutine MatElm(x,F)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  real*8                          :: x(6),F

  g = x
  call SFmatel
  if (method==4) call Rate
  F=SF

End Subroutine MatElm
