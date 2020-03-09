Subroutine WriteMOs(mol)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(1),intent(in)         :: mol

  if (mol == 'a') call rWriteMOs(nA,nNAOA,nOnAtmA,OnAtmA)
  if (mol == 'b') call rWriteMOs(nB,nNAOB,nOnAtmB,OnAtmB)

End Subroutine WriteMOs
!------------------------------------------------------------------------------
Subroutine rWriteMOs(n,nNAO,nOnAtm,OnAtm)
  implicit none

  integer,intent(in)              :: n,nNAO,nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)

  integer                         :: MOs

  integer                         :: i,j,k

  write(6,*)
  write(6,'(a)')'Order of MO coefficients:'
  write(6,'(a)')'    C    NAO   ATOM'
  MOs=0
  do i=1,n
    do j=1,nNAO
      do k=1,OnAtm(j,0)
        if (i == OnAtm(j,k)) then
          MOs=MOs+1
          write(6,'(i5,i7,i7)')MOs,j,OnAtm(j,k)
        endif
      enddo
    enddo
  enddo
  write(6,'(a,i5)'),'Number of MO coefficients=',MOs
  write(6,*)

End Subroutine rWriteMOs
