Subroutine DFPmin(p,n,gtol,fret,test2,conv1,conv2)
! version 3.0     July 18, 2019
  implicit none

  integer,parameter               :: NMAX=6                  ! maximum anticipated value of n
  integer,parameter               :: ITMAX=25                ! maximum allowed number of iterations
  real*8,parameter                :: EPS=1.0d-6,TOLX=4.*EPS  ! convergence criterion on x values
  real*8,parameter                :: ZERO=0.0d0,ONE=1.0d0
  integer                         :: iter,n
  real*8                          :: fret,gtol,p(n),func
  logical                         :: conv1,conv2
  integer                         :: i,its,j
  real*8                          :: den,fac,fad,fae,fp,sumdg,sumxi,temp,test1,test2
  real*8                          :: dg(NMAX),g(NMAX),hdg(NMAX),hessin(NMAX,NMAX)
  real*8                          :: pnew(NMAX),xi(NMAX),a
  real*8                          :: gNorm

! USES dfunc,func,linmin

! Given a starting point p(1:n) that is a vector of length n, the Broyden-Fletcher-GoldfarbShanno variant of
! Davidon-Fletcher-Powell minimization is performed on a function func, using its gradient as calculated by a routine dfunc.
! The convergence requirement on ZEROing the gradient is input as gtol.
! Returned quantities are p(1:n) (the location of the minimum), iter (the number of iterations that were performed),
! and fret (the minimum value of the function). The routine linmin is called to perform approximate line minimizations.

  fp=func(p)                                ! Calculate starting function value and gradient
  call dfunc(p,g)

  call CalcHess(p,hessin)                   ! Initialize the inverse Hessian
  do i=1,n                                  ! Calculate the initial line direction (-Hinv*g)
    xi(i)=ZERO
    do j=1,n
      xi(i)=xi(i)-Hessin(i,j)*g(j)
    enddo
  enddo

  do its=1,ITMAX                            ! MAIN LOOP over the iterations.
    iter=its

    call LinMin(p,xi,pnew,n,fret)

    fp=fret
    do i=1,n
      xi(i)=pnew(i)-p(i)                    ! Update the line direction,
      p(i)=pnew(i)                          ! and the current point.
    enddo

    test1=ZERO                              ! Test for convergence on delta(x).
    do i=1,n
      temp=abs(xi(i))/max(abs(p(i)),ONE)
      if(temp.gt.test1)test1=temp
    enddo

    do i=1,n                                ! Save the old gradient,
      dg(i)=g(i)
    enddo

    call dfunc(p,g)                         ! and get the new gradient.
    fret=func(p)

    gNorm=ZERO                              ! Test for convergence on zero gradient.
    do i=1,n
      gNorm=gNorm+g(i)*g(i)
    enddo
    gNorm=sqrt(gNorm)
    test2=gNorm

    conv1=test1.lt.TOLX
    conv2=test2.lt.gtol

    if(conv1.and.conv2) RETURN
    
    do i=1,n                                ! Compute difference of gradients,
      dg(i)=g(i)-dg(i)
    enddo
    do i=1,n                                ! and difference times current matrix.
      hdg(i)=ZERO
      do j=1,n
        hdg(i)=hdg(i)+hessin(i,j)*dg(j)
      enddo
    enddo
    fac=ZERO                                ! Calculate dot products for the denominators.
    fae=ZERO
    sumdg=ZERO
    sumxi=ZERO
    do i=1,n
      fac=fac+dg(i)*xi(i)
      fae=fae+dg(i)*hdg(i)
      sumdg=sumdg+dg(i)**2
      sumxi=sumxi+xi(i)**2
    enddo
    if(fac.gt.sqrt(EPS*sumdg*sumxi))then    ! Skip update if fac not sufficiently positive.
      fac=ONE/fac
      fad=ONE/fae
      do i=1,n                              ! The vector that makes BFGS different from DFP:
        dg(i)=fac*xi(i)-fad*hdg(i)
      enddo
      do i=1,n                              ! The BFGS updating formula:
        do j=i,n
          hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j)
          hessin(j,i)=hessin(i,j)
        enddo
      enddo
    endif
    do i=1,n                                ! Now calculate the next direction to go,
      xi(i)=ZERO
      do j=1,n
        xi(i)=xi(i)-hessin(i,j)*g(j)
      enddo
    enddo
  enddo                                     ! and GO BACK for another iteration.

! stop 'too many iterations in dfpmin'

End Subroutine DFPmin
