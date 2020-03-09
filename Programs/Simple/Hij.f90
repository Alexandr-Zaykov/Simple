Subroutine Hij(MOA,nAOA,MOB,nAOB,SaoAB,HiiA,HiiB,H)
! version 3.0     July 18, 2019
  implicit none

  integer,intent(in)              :: nAOA,nAOB
  real*8,intent(in)               :: MOA(nAOA,2),MOB(nAOB,2)
  real*8,intent(in)               :: SaoAB(nAOA,nAOB)
  real*8,intent(in)               :: HiiA(nAOA),HiiB(nAOB)

  real*8,intent(out)              :: H(2,2)

! From Linderberg formula K for distances > sum of vdW radii =~ 1.03
  real*8                          :: K=0.515 ! 1/2 * 1.03
  integer                         :: i,j,m,n

  H = 0.0d0

  do i=1,2
    do j=1,2
      do m=1,nAOA
        do n=1,nAOB
          H(i,j)=H(i,j)+MOA(m,i)*MOB(n,j)*K*(HiiA(m)+HiiB(n))*SaoAB(m,n)
        enddo
      enddo
    enddo
  enddo

End Subroutine Hij
