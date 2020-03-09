Subroutine AOlabels
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rAOlabels(nA,nNAOA,nOnAtmA,OnAtmA,DimnContA,maxnCA,TypContA,nAOA,AOsA)
  call rAOlabels(nB,nNAOB,nOnAtmB,OnAtmB,DimnContB,maxnCB,TypContB,nAOB,AOsB)

End Subroutine AOlabels
!------------------------------------------------------------------------------
Subroutine rAOlabels(n,nNAO,nOnAtm,OnAtm,DimnCont,maxnC,TypCont,nAO,AOs)
  implicit none

  integer,intent(in)              :: n
  integer,intent(in)              :: nNAO
  integer,intent(in)              :: nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(in)              :: DimnCont(nNAO)
  integer,intent(in)              :: maxnC
  character(1),intent(in)         :: TypCont(nNAO,maxnC)
  integer,intent(in)              :: nAO

  character(3),intent(out)        :: AOs(nAO)

  integer                         :: i,j,k,l,m
  integer                         :: iAO
  integer                         :: ic,nn
  character(1)                    :: tmp

  do m=1,2
    iAO=0
    do i=1,n
      do j=1,nNAO
        do k=1,OnAtm(j,0)
          if (i == OnAtm(j,k)) then
            do l=1,DimnCont(j)
              tmp=TypCont(j,l)
              call lowercase(tmp)
              if (tmp=='s') then
                iAO=iAO+1
                AOs(iAO)='s  '
              endif
              if (tmp=='p') then
                iAO=iAO+1
                AOs(iAO)='px '
                iAO=iAO+1
                AOs(iAO)='py '
                iAO=iAO+1
                AOs(iAO)='pz '
              endif
              if (tmp=='d') then
                iAO=iAO+1
                AOs(iAO)='dxx'
                iAO=iAO+1
                AOs(iAO)='dyy'
                iAO=iAO+1
                AOs(iAO)='dzz'
                iAO=iAO+1
                AOs(iAO)='dxy'
                iAO=iAO+1
                AOs(iAO)='dxz'
                iAO=iAO+1
                AOs(iAO)='dyz'
              endif
            enddo
          endif
        enddo
      enddo
    enddo
  enddo

End Subroutine rAOlabels
