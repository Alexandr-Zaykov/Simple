Subroutine CalcGrad(p0,grad)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  real*8                          :: p(6),p0(6),F,grad(6)
  real*8                          :: T(-2:2)
  real*8                          :: dx(6)
  real*8                          :: dxRed=1.0D-2
  real*8                          :: dxSave(6)=(/1.0D-3,1.0D-3,1.0D-3,1.0D-5,1.0D-5,1.0D-5/) ! (3 translations, 3 rotations)
  integer                         :: i,j

  dx = dxSave
  if (last) dx = dx * dxRed

  do i=1,6
    do j=-2,2
      if (j/=0) then
        p = p0
        p(i)=p(i)+j*dx(i)
        call MatElm(p,F)
        T(j)=F
      endif
    enddo
    grad(i)=(-T(2)+8.0d0*T(1)-8.0d0*T(-1)+T(-2))/(12.0d0*dx(i))
  enddo

End Subroutine CalcGrad
