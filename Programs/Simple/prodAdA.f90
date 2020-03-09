Subroutine prodAdA(B,d,A)
! version 3.0     July 18, 2019
  implicit none

! Multiplication B=A*d*A', A is matrix, d is vector

  integer,parameter               :: N=2

  real*8,intent(in)               :: A(N,N), d(N)
  real*8,intent(out)              :: B(N,N)
  integer                         :: i,j,k

  do i=1,N
    do j=1,N
      B(i,j)=0.0d0
      do k=1,N
        B(i,j)=B(i,j)+d(k)*A(i,k)*A(j,k)
      enddo
    enddo
  enddo

  return

End Subroutine prodAdA
