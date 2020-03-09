function isblank(line)
! version 3.0     July 18, 2019
  implicit none

  integer                         :: isblank
  character(*),intent(in)         :: line
  integer                         :: i,j

  j=len_trim(line)
  isblank=j+1

  do i=1,j
    if (line(i:i)==' ') then
      isblank=i
      return
    endif
  enddo

  return
End function isblank
