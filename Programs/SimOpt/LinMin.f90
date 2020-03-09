Subroutine LinMin(p,xi,pnew,n,fret)
! version 3.0     July 18, 2019
  implicit none

  external                        :: df1dim
  real*8,external                 :: dbrent,f1dim
  integer                         :: n
  real*8                          :: fret,p(n),pnew(n),xi(n)
  integer,parameter               :: NMAX=6      ! Maximum anticipated n
  real*8,parameter                :: TOL=1.0d-2  ! TOL passed to brent
  real*8,parameter                :: ZERO=0.0d0,ONE=1.0d0
  real*8                          :: ax,bx,xx,fa,fb,fx,xmin
  real*8                          :: pcom(NMAX),xicom(NMAX)
  integer                         :: ncom
  common /f1com/                     pcom,xicom,ncom
  real*8                          :: xLim(6)=(/1.0,1.0,1.0,0.1,0.1,0.1/)
  real*8                          :: A
  integer                         :: j
  
  ncom=n                                    ! Set up the common block.
  do j=1,n
    pcom(j)=p(j)
    xicom(j)=xi(j)
  enddo

  ax=ZERO                                   ! Initial guess for brackets.
  xx=ONE
  call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
  fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)

  A=ONE                                     ! Limit the displacement to xLim
  do j=1,n
    A=min(A,xLim(j)/abs(xmin*xi(j)))
  enddo

  do j=1,n                                  ! Construct the vector results to return.
    xi(j)=A*xmin*xi(j)
    pnew(j)=p(j)+xi(j)
  enddo

  return

End Subroutine LinMin
