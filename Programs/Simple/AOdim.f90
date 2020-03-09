Subroutine AOdim
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rAOdim(nA,nNAOA,nOnAtmA,OnAtmA,DimnContA,nAOA,maxnCA,TypContA,ngA,nPGA)
  call rAOdim(nB,nNAOB,nOnAtmB,OnAtmB,DimnContB,nAOB,maxnCB,TypContB,ngB,nPGB)

End Subroutine AOdim
!------------------------------------------------------------------------------
Subroutine rAOdim(n,nNAO,nOnAtm,OnAtm,DimnCont,nAO,maxnC,TypCont,ng,nPG)
  implicit none

  integer,intent(in)              :: n
  integer,intent(in)              :: nNAO
  integer,intent(in)              :: nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(in)              :: DimnCont(nNAO)
  integer,intent(in)              :: maxnC
  character(1),intent(in)         :: TypCont(nNAO,maxnC)
  integer,intent(in)              :: ng(nNAO,maxnC)

  integer,intent(out)             :: nAO,nPG

  integer                         :: i,j,k,l,m
  character(1)                    :: tmp

  nAO=0
  nPG=0

  do i=1,n
    do j=1,nNAO
      do k=1,OnATM(j,0)
        if (i == OnATM(j,k)) then
          do l=1,DimnCont(j)
            tmp=TypCont(j,l)
            call lowercase(tmp)
            if (tmp=='s')  nAO=nAO+1
            if (tmp=='p')  nAO=nAO+3
            if (tmp=='d')  nAO=nAO+6
            do m=1,ng(j,l)
              if (tmp=='s')  nPG=nPG+1
              if (tmp=='p')  nPG=nPG+3
              if (tmp=='d')  nPG=nPG+6
            enddo
          enddo
        endif
      enddo
    enddo
  enddo

End Subroutine rAOdim
