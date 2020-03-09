function ReadLine(line)
! version 3.0     July 18, 2019
  implicit none

  integer                         :: ReadLine

  integer,external                :: nonblank

  character(*),intent(inout)      :: line

1 read(3,'(a)',end=2) line

  if (nonblank(1,line)==0) goto 1
  if (line(1:1)=='#') goto 1
  ReadLine=0 ! line read
  if (index(line,':')>0) ReadLine=2 ! line with :

  return

2 ReadLine=1 ! End of file

  return
End function ReadLine
