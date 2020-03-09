Subroutine SwitchCol(N,c1,c2,A)
! version 3.0     July 18, 2019
  implicit none

  integer,intent(in)              :: N,c1,c2
  real*8,intent(inout)            :: A(N,N)
  real*8                          :: tmp
  integer                         :: i

  do i=1,N
    tmp=A(i,c1)
    A(i,c1)=A(i,c2)
    A(i,c2)=tmp
  enddo

End Subroutine SwitchCol
