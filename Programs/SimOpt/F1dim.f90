Function F1dim(x)
! version 3.0     July 18, 2019
  implicit none

  real*8                          :: F1dim,func
  real*8                          :: x
  integer,parameter               :: NMAX=6
  real*8                          :: pcom(NMAX),xicom(NMAX),xt(NMAX)
  integer                         :: ncom
  common /f1com/                     pcom,xicom,ncom
  integer                         :: j


  do j=1,ncom
    xt(j)=pcom(j)+x*xicom(j)
  enddo
  f1dim=func(xt)

  return

END

