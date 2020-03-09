Subroutine FirstGuess(g)
! version 3.0     July 18, 2019
  implicit none

  real*8,external                 :: F1dim
  real*8                          :: gNorm,xmin
  integer,parameter               :: steps=20
  real*8,parameter                :: ZERO=0.0d0,ONE=1.0d0,extent=1.5D0
  real*8                          :: step,Fmin,x,F
  real*8                          :: g(6),grad(6),p(6)
  integer                         :: j
  integer,parameter               :: NMAX=6
  real*8                          :: pcom(NMAX),xicom(NMAX),xt(NMAX)
  integer                         :: ncom
  common /f1com/                     pcom,xicom,ncom

  p = g
  call CalcGrad(p,grad)
  ncom=6
  do j=1,ncom                        ! Fill the common block
    pcom(j)=g(j)
    xicom(j)=-grad(j)
  enddo
  
  gNorm=ZERO                         ! Calculate norm of gradient
  do j=1,6
    gNorm=gNorm+grad(j)*grad(j)
  enddo
  gNorm=sqrt(gNorm)

  xmin=ZERO                          ! Set step size for linearch search
  step=ONE/(gNorm*steps)
  x=ZERO                             ! Initial step
  Fmin=F1Dim(x)

  do x=ZERO+step,extent/gNorm,step   ! Main loop for minimum finding
    F=F1dim(x)
    if(F<Fmin)then
      xmin=x
      Fmin=F
    endif
  enddo

  do j=1,6                           ! Save geometry of found minimum
    g(j)=g(j)-xmin*grad(j)
  enddo
 
  return

End Subroutine FirstGuess
