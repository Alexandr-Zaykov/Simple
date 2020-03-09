Subroutine WriteMolden(fret)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  real*8                          :: fret
  integer                         :: i,j

 write(9,'(i4)')nA+nB
 write(9,'(a,es15.7)')'F=',fret
 do i=1,nA
   write(9,'(a,3f15.7)')labA(i),(xA(j,i),j=1,3)
 enddo
 do i=1,nB
   write(9,'(a,3f15.7)')labB(i),(xB(j,i),j=1,3)
 enddo

 return

End Subroutine WriteMolden
