Subroutine ReadNPrim(N,nP)
! version 3.0     July 18, 2019
  implicit none

  integer,intent(in)              :: N

  integer,intent(out)             :: nP

  integer,external                :: ReadLine

  character(120)                  :: line
  character(1)                    :: ctmp
  integer                         :: gtmp
  integer                         :: i,j,k

  nP=0

  do i=1,N
    j=ReadLine(line)
    read(line,*)ctmp,gtmp
    nP=nP+gtmp
    do j=1,gtmp
      k=ReadLine(line)
    enddo
  enddo

  return
End Subroutine ReadNPrim
