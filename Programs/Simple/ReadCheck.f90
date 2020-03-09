Function ReadCheck()
! version 3.0     July 18, 2019
  implicit none

  logical                         :: ReadCheck
  character(120)                  :: line

  rewind(3)

  ReadCheck=.false.
1 read(3,'(a)',end=2) line
    if (line(1:1)=='#') goto 1
    call lowercase(line)
    if (index(line,'check')>0) then
      ReadCheck=.true.
      goto 2
    endif
  goto 1

2 continue

  rewind(3)

End Function ReadCheck
