Function nonblank(l,line)
! version 3.0     July 18, 2019
  implicit none

  integer                         :: nonblank
  integer,intent(in)              :: l
  character(*),intent(in)         :: line
  integer                         :: j,i

  j=len_trim(line)
  nonblank=0

  do i=l,j
    if (line(i:i) /= ' ') then
      nonblank=i
      return
    endif
  enddo

  return
End Function nonblank
