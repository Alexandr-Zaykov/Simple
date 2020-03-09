Subroutine TrafoGeom(x,n,g,TMP)
! version 3.0     July 18, 2019
  implicit none

  integer,intent(in)              :: n
  real*8,intent(inout)            :: x(3,n)
  real*8,intent(in)               :: g(6)
  real*8,intent(inout)            :: TMP(3,n)

  real*8                          :: T(3,3)
  real*8                          :: d2r
  real*8                          :: a(6)
  integer                         :: i,j

  d2r = acos(-1.0d0)/180.0d0

  a = g
  do i=1,3
    a(i)=a(i)*d2r
  enddo

! T = Z * Y * X * Rz * Ry * Rx (X' = T * X)

  T(1,1) =  cos(a(2))*cos(a(3))
  T(1,2) = -cos(a(1))*sin(a(3)) + sin(a(1))*sin(a(2))*cos(a(3))
  T(1,3) =  sin(a(1))*sin(a(3)) + cos(a(1))*sin(a(2))*cos(a(3))
  T(2,1) =  cos(a(2))*sin(a(3))
  T(2,2) =  cos(a(1))*cos(a(3)) + sin(a(1))*sin(a(2))*sin(a(3))
  T(2,3) = -sin(a(1))*cos(a(3)) + cos(a(1))*sin(a(2))*sin(a(3))
  T(3,1) = -sin(a(2))
  T(3,2) =  sin(a(1))*cos(a(2))
  T(3,3) =  cos(a(1))*cos(a(2))

! Rotation

  TMP= matmul(T,X)

! Translation
  do i=1,n
    do j=1,3
      X(j,i)=TMP(j,i)+a(j+3)
    enddo
  enddo

  return

End Subroutine TrafoGeom
