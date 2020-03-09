Subroutine GetConst()
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line

  g = 0.0d0

1 read(3,'(a)',end=9)line
    if (line(1:1) == '#') goto 1
    if (index(line,':') > 0) then
      call lowercase(line)
      if (index(line,'const') > 0) then
2       read(3,'(a)',end=9)line
          if (line(1:1) == '#') goto 2
          if (index(line,':') > 0) goto 9
          if (line(1:3) == 'R X') read(line(4:),*)g(1)
          if (line(1:3) == 'R Y') read(line(4:),*)g(2)
          if (line(1:3) == 'R Z') read(line(4:),*)g(3)
          if (line(1:3) == 'T X') read(line(4:),*)g(4)
          if (line(1:3) == 'T Y') read(line(4:),*)g(5)
          if (line(1:3) == 'T Z') read(line(4:),*)g(6)
        goto 2
      endif
    endif
  goto 1

9 continue

  rewind(3)

End Subroutine GetConst
