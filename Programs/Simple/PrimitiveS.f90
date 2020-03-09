Subroutine PrimitiveS(Ax,Bx,alpha,beta,a,b,pi,S)
! version 3.0     July 18, 2019
  implicit none

  real*8,external                 :: Ovl
  real*8,intent(in)               :: Ax(3),Bx(3)
  real*8,intent(in)               :: alpha,beta
  integer,intent(in)              :: a(3),b(3)
  real*8,intent(in)               :: pi
  real*8,intent(out)              :: S
  real*8                          :: sx,sy,sz
  real*8                          :: Eab,f

  sx=Ovl(Ax(1),Bx(1),alpha,beta,a(1),b(1))
  sy=Ovl(Ax(2),Bx(2),alpha,beta,a(2),b(2))
  sz=Ovl(Ax(3),Bx(3),alpha,beta,a(3),b(3))

  f=pi/(alpha+beta)
  f=f*f*f
  f=sqrt(f)
  Eab=(Ax(1)-Bx(1))*(Ax(1)-Bx(1))+(Ax(2)-Bx(2))*(Ax(2)-Bx(2))+(Ax(3)-Bx(3))*(Ax(3)-Bx(3))
  Eab=exp(-Eab*alpha*beta/(alpha+beta))
  S=Eab*f*sx*sy*sz

End Subroutine PrimitiveS
