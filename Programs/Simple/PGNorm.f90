Subroutine PGNorm
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rPGNorm(nPGA,aPrimA,alphaPrimA,GnormA)
  call rPGNorm(nPGB,aPrimB,alphaPrimB,GnormB)

End Subroutine PGNorm
!------------------------------------------------------------------------------
Subroutine rPGNorm(nPG,aPrim,alphaPrim,Gnorm)
  implicit none

  integer,intent(in)              :: nPG
  integer,intent(in)              :: aPrim(nPG,3)
  real*8,intent(in)               :: alphaPrim(nPG)

  real*8,intent(out)              :: Gnorm(nPG)

  real*8                          :: pi
  integer                         :: aa(3)
  integer                         :: i,j

  pi=4.0d0*atan(1.0d0)

  do i=1,nPG
    do j=1,3
      aa(j)=aPrim(i,j)
    enddo
    call Norm(alphaPrim(i),aa,pi,Gnorm(i))
  enddo

End Subroutine rPGNorm
