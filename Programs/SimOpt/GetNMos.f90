Subroutine GetNMOs
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rGetNMOs(nNAOA,nOnAtmA,OnAtmA,MOsA)
  call rGetNMOs(nNAOB,nOnAtmB,OnAtmB,MOsB)

End Subroutine GetNMOs
!------------------------------------------------------------------------------
Subroutine rGetNMOs(nNAO,nOnAtm,OnAtm,MOs)
  implicit none

  integer,intent(in)              :: nNAO,nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(out)             :: MOs

  integer                         :: i,j

  MOs=0
  do i=1,nNAO
    do j=1,OnAtm(i,0)
      MOs=MOs+1
    enddo
  enddo

End Subroutine rGetNMOs
