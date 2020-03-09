Subroutine TypIndex
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rTypIndex(nA,nNAOA,nOnAtmA,OnAtmA,DimnContA,maxnCA,TypContA,nAOA,a)
  call rTypIndex(nB,nNAOB,nOnAtmB,OnAtmB,DimnContB,maxnCB,TypContB,nAOB,b)

End Subroutine TypIndex
!------------------------------------------------------------------------------
Subroutine rTypIndex(n,nNAO,nOnAtm,OnAtm,DimnCont,maxnC,TypCont,nAO,a)
  implicit none

  integer,intent(in)              :: n
  integer,intent(in)              :: nNAO
  integer,intent(in)              :: nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(in)              :: DimnCont(nNAO)
  integer,intent(in)              :: maxnC
  character(1),intent(in)         :: TypCont(nNAO,maxnC)
  integer,intent(in)              :: nAO

  integer,intent(out)             :: a(nAO,3)

  integer                         :: i,j,k,l
  integer                         :: iAO
  character(1)                    :: tmp

  a = 0

  iAO=0
  do i=1,n
    do j=1,nNAO
      do k=1,OnATM(j,0)
        if (i == OnATM(j,k)) then
          do l=1,DimnCont(j)
            tmp=TypCont(j,l)
            call lowercase(tmp)
            iAO=iAO+1
            if (tmp=='p') then
              a(iAO,1)=1
              iAO=iAO+1
              a(iAO,2)=1
              iAO=iAO+1
              a(iAO,3)=1
            endif
            if (tmp=='d') then
              a(iAO,1)=2
              iAO=iAO+1
              a(iAO,2)=2
              iAO=iAO+1
              a(iAO,3)=2
              iAO=iAO+1
              a(iAO,1)=1; a(iAO,2)=1
              iAO=iAO+1
              a(iAO,1)=1; a(iAO,3)=1
              iAO=iAO+1
              a(iAO,2)=1; a(iAO,3)=1
            endif
          enddo
        endif
      enddo
    enddo
  enddo

End Subroutine rTypIndex
