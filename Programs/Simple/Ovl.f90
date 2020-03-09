real*8 Function Ovl(Ax,Bx,alpha,beta,a,b)
! version 3.0     July 18, 2019
  implicit none

  real*8,intent(in)               :: Ax,Bx,alpha,beta
  integer,intent(in)              :: a,b
! s(0:a+b,0:b), a<=2, b<=2 up to fo d-functions
  real*8                          :: s(0:4,0:2)
  real*8                          :: P
  integer                         :: i,j

  P=(alpha*Ax+beta*Bx)/(alpha+beta)
  s(0,0)=1.0d0
  s(1,0)=P-Ax

! Recursion relation
  do i=2,a+b
    s(i,0)=(P-Ax)*s(i-1,0)+(i-1)*s(i-2,0)/(2.0d0*(alpha+beta))
  enddo
! Transfer equation
  do j=1,b
    do i=0,a+b-j
      s(i,j)=s(i+1,j-1)+(Ax-Bx)*s(i,j-1)
    enddo
  enddo

  Ovl=s(a,b)

End Function Ovl
