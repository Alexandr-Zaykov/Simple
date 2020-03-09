Subroutine Norm(alpha,a,pi,N)
! version 3.0     July 18, 2019
  implicit none

  real*8                          :: dfact(0:10) = [1.0d0,1.0d0,3.0d0,15.0d0,105.0d0,945.0d0,10395.0d0, &
                                                    135135.0d0,2027025.0d0,34459425.0d0,654729075.0d0]
  real*8,intent(in)               :: alpha
  integer,intent(in)              :: a(3)
  real*8,intent(in)               :: pi

  real*8,intent(out)              :: N

  real*8                          :: D,F

  N=(4.0d0*alpha)**((a(1)+a(2)+a(3))/2.0d0)
  D=sqrt(dfact(a(1))*dfact(a(2))*dfact(a(3)))
  F=2.0d0*alpha/pi
  F=F*F*F
  F=sqrt(F)
  F=sqrt(F)
  N=F*N/D

End Subroutine Norm
