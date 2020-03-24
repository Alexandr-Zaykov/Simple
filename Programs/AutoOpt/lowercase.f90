Subroutine lowercase(line)
! version 1.7     May 12, 2017
  implicit none

! Convert upercase characters in line into lowercase

  character(*),intent(inout)      :: line
  integer                         :: jgap
  integer                         :: i

  jgap=ichar('a')-ichar('A')

  do i=1,len_trim(line)
    if (line(i:i) <= 'Z') then
      if (line(i:i) >= 'A') line(i:i)=char(ichar(line(i:i))+jgap)
    endif
  enddo

  return

End Subroutine lowercase
