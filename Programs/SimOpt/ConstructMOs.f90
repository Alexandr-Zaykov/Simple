Subroutine ConstructMOs(mol)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(1),intent(in)         :: mol

  if (mol == 'a') call rConstructMOs(MOsA,nNAOA,nOnAtmA,OnAtmA,DimnContA,maxnCA,TypContA,CContA,CA,nAOA,MOA)
  if (mol == 'b') call rConstructMOs(MOsB,nNAOB,nOnAtmB,OnAtmB,DimnContB,maxnCB,TypContB,CContB,CB,nAOB,MOB)

End Subroutine ConstructMOs
!------------------------------------------------------------------------------
Subroutine rConstructMOs(n,nNAO,nOnAtm,OnAtm,DimnCont,maxnC,TypCont,CCont,C,nAO,MO)

  integer,intent(in)              :: n
  integer,intent(in)              :: nNAO
  integer,intent(in)              :: nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(in)              :: DimnCont(nNAO)
  integer,intent(in)              :: maxnC
  character(1),intent(in)         :: TypCont(nNAO,maxnC)
  real*8,intent(in)               :: CCont(nNAO,maxnC,6)
  real*8,intent(in)               :: C(n,2)
  integer,intent(in)              :: nAO

  real*8,intent(out)              :: MO(nAO,2)

  integer                         :: i,j,k,l,m
  integer                         :: iAO
  integer                         :: ic,nn
  character(1)                    :: tmp
  integer                         :: MOs

  integer                         :: OccAt

  MO = 0.0d0


!Find last atom with nonblank NAO
  OccAt=0
  do j=1,nNAO
     do k=1,OnAtm(j,0)
       if (OnAtm(j,k) .gt. OccAt) then
         OccAt=OnAtm(j,k)
       endif
     enddo
  enddo

  do m=1,2
    MOs=0
    iAO=0
    do i=1,OccAt
      do j=1,nNAO
       do k=1,OnAtm(j,0)
          if (i == OnAtm(j,k)) then
            MOs=MOs+1
            do l=1,DimnCont(j)
              tmp=TypCont(j,l)
              call lowercase(tmp)
              if (tmp=='s') nn=1
              if (tmp=='p') nn=3
              if (tmp=='d') nn=6
              do ic=1,nn
                iAO=iAO+1
                MO(iAO,m)=MO(iAO,m)+C(MOs,m)*CCont(j,l,ic)
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
  enddo


End Subroutine rConstructMOs
