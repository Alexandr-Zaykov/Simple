Subroutine NegCol(N,c1,A)
! version 3.0     July 18, 2019
  implicit none

  integer,intent(in)              :: N,c1
  real*8,intent(inout)            :: A(N,N)
  integer                         :: i

  do i=1,N
    A(i,c1)=-A(i,c1)
  enddo

End Subroutine NegCol
