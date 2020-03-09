Subroutine ReadReduce
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line
  integer                         :: n

  Reduce=0.0d0

  rewind(3)

1 read(3,'(a)',end=99)line
    if (line(1:1) == '#') goto 1
    if (index(line,':') > 0) then
      call lowercase(line)
      if (index(line,'red') > 0) then
        n=index(line,':')+1
        read(line(n:),*)Reduce
        goto 1
      endif
    endif
  goto 1

99 continue

  rewind(3)

  return

End Subroutine ReadReduce
