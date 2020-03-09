Function dF1dim(x)
! version 3.0     July 18, 2019
  implicit none

  real*8                          :: dF1dim
  real*8                          :: x
  integer,parameter               :: NMAX=6
  real*8                          :: pcom(NMAX),xicom(NMAX)
  integer                         :: ncom
  common /f1com/                     pcom,xicom,ncom
  integer                         :: j
  real*8                          :: df(NMAX),xt(NMAX)

  do j=1,ncom
    xt(j)=pcom(j)+x*xicom(j)
  enddo
  call dfunc(xt,df)
  df1dim=0.0d0
  do j=1,ncom
    df1dim=df1dim+df(j)*xicom(j)
  enddo

  return

END
