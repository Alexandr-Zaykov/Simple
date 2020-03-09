Subroutine AOPrim
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rAOPrim(nA,nNAOA,nOnAtmA,OnAtmA,DimnContA,maxnCA,TypContA,ngA,nAOA,nAOPrimA,nXA)
  call rAOPrim(nB,nNAOB,nOnAtmB,OnAtmB,DimnContB,maxnCB,TypContB,ngB,nAOB,nAOPrimB,nXB)

End Subroutine AOPrim
!------------------------------------------------------------------------------
Subroutine rAOPrim(n,nNAO,nOnAtm,OnAtm,DimnCont,maxnC,TypCont,ng,nAO,nAOPrim,nX)
  implicit none

  integer,intent(in)              :: n
  integer,intent(in)              :: nNAO
  integer,intent(in)              :: nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(in)              :: DimnCont(nNAO)
  integer,intent(in)              :: maxnC
  character(1),intent(in)         :: TypCont(nNAO,maxnC)
  integer,intent(in)              :: nAO
  integer,intent(in)              :: ng(nNAO,maxnC)

  integer,intent(out)             :: nAOPrim(nAO)
  integer,intent(out)             :: nX(nAO)

  integer                         :: i,j,k,l,m
  integer                         :: iAO
  character(1)                    :: tmp
  character(1)                    :: S
  character(2)                    :: Px(3)
  character(3)                    :: Dx(6)
  integer                         :: nn
  integer                         :: ii

  data S  /'s'/
  data Px /'px','py','pz'/
  data Dx /'dxx','dyy','dzz','dxy','dxz','dyz'/

  ii=0
  iAO=0
  do i=1,n
    do j=1,nNAO
      do k=1,OnAtm(j,0)
        if (i == OnAtm(j,k)) then
          do l=1,DimnCont(j)
            tmp=TypCont(j,l)
            call lowercase(tmp)
            if (tmp=='s') nn=1
            if (tmp=='p') nn=3
            if (tmp=='d') nn=6
            do m=1,nn
             ii=ii+1
             nAOPrim(ii)=ng(j,l)
             nX(ii)=i
            enddo
            iAO=iAO+nn
          enddo
        endif
      enddo
    enddo
  enddo

End Subroutine rAOPrim
