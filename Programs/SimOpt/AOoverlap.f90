Subroutine AOoverlap(mol)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(2),intent(in)         :: mol

  if (mol == 'aa') call rAOoverlap(nA,xA,nPGA,aPrimA,alphaPrimA,GnormA,gPrimA,nNAOA,nAOA,nAOPrimA,nXA,  &
                                  nA,xA,nPGA,aPrimA,alphaPrimA,GnormA,gPrimA,nNAOA,nAOA,nAOPrimA,nXA,SaoA)
  if (mol == 'bb') call rAOoverlap(nB,xB,nPGB,aPrimB,alphaPrimB,GnormB,gPrimB,nNAOB,nAOB,nAOPrimB,nXB,  &
                                  nB,xB,nPGB,aPrimB,alphaPrimB,GnormB,gPrimB,nNAOB,nAOB,nAOPrimB,nXB,SaoB)
  if (mol == 'ab') call rAOoverlap(nA,xA,nPGA,aPrimA,alphaPrimA,GnormA,gPrimA,nNAOA,nAOA,nAOPrimA,nXA,  &
                                  nB,xB,nPGB,aPrimB,alphaPrimB,GnormB,gPrimB,nNAOB,nAOB,nAOPrimB,nXB,SaoAB)

End Subroutine AOoverlap
!------------------------------------------------------------------------------
Subroutine rAOoverlap(nA,xA,nPGA,aPrimA,alphaPrimA,GnormA,gPrimA,nNAOA,nAOA,nAOPrimA,nXA,  &
                      nB,xB,nPGB,aPrimB,alphaPrimB,GnormB,gPrimB,nNAOB,nAOB,nAOPrimB,nXB,S)
  implicit none

  integer,intent(in)              :: nA,nB
  real*8,intent(in)               :: xA(3,nA),xB(3,nB)
  integer,intent(in)              :: nPGA,nPGB
  integer,intent(in)              :: aPrimA(nPGA,3),aPrimB(nPGB,3)
  real*8,intent(in)               :: alphaPrimA(nPGA),alphaPrimB(nPGB)
  real*8,intent(in)               :: GnormA(nPGA),GnormB(nPGB)
  real*8,intent(in)               :: gPrimA(nPGA),gPrimB(nPGB)
  integer,intent(in)              :: nNAOA,nNAOB
  integer,intent(in)              :: nAOA,nAOB
  integer,intent(in)              :: nAOPrimA(nAOA),nAOPrimB(nAOB)
  integer,intent(in)              :: nXA(nAOA),nXB(nAOB)

  real*8,intent(out)              :: S(nAOA,nAOB)

  real*8                          :: pi
  integer                         :: ii,jj
  integer                         :: i,j,k,l,m
  real*8                          :: tmpxA(3),tmpxB(3)
  integer                         :: aa(3),bb(3)
  real*8                          :: Sprim

  S = 0.0d0

  pi=4.0d0*atan(1.0d0)

! Block A
  ii=0
  do i=1,nAOA
    do k=1,nAOPrimA(i)
      ii=ii+1
      do m=1,3
        tmpxA(m)=xA(m,nXA(i))/0.5291772086D+00
        aa(m)=aPrimA(ii,m)
      enddo
! Block B
      jj=0
      do j=1,nAOB
        do l=1,nAOPrimB(j)
          jj=jj+1
          do m=1,3
            tmpxB(m)=xB(m,nXB(j))/0.5291772086D+00
            bb(m)=aPrimB(jj,m)
          enddo
          call PrimitiveS(tmpxA,tmpxB,alphaPrimA(ii),alphaPrimB(jj),aa,bb,pi,Sprim)
          S(i,j)=S(i,j)+gPrimA(ii)*GnormA(ii)*gPrimB(jj)*GnormB(jj)*Sprim
        enddo
      enddo
! End Block B
    enddo
  enddo
! End Block A

End Subroutine rAOoverlap
