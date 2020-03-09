Subroutine MOoverlap(MOA,nAOA,MOB,nAOB,SaoAB,Smo)
! version 3.0     July 18, 2019
  implicit none

  integer,intent(in)              :: nAOA,nAOB
  real*8,intent(in)               :: MOA(nAOA,2),MOB(nAOB,2)
  real*8,intent(in)               :: SaoAB(nAOA,nAOB)

  real*8,intent(out)              :: Smo(2,2)

  integer                         :: i,j,m,n

  Smo = 0.0d0

  do i=1,2
    do j=1,2
      do m=1,nAOA
        do n=1,nAOB
          Smo(i,j)=Smo(i,j)+MOA(m,i)*MOB(n,j)*SaoAB(m,n)
        enddo
      enddo
    enddo
  enddo

End Subroutine MOoverlap
