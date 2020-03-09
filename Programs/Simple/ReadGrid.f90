Subroutine ReadGrid
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line

  grid = .false.
  fidelity = 2

1 read(3,'(a)',end=9)line
    if (line(1:1) == '#') goto 1
    if (index(line,':') > 0) then
      call lowercase(line)
      if (index(line,'grid') > 0) then
        grid=.true.
2       read(3,'(a)',end=9)line
          if (line(1:1) == '#') goto 2
          if (index(line,':') > 0) goto 9
          if (index(line,'fine')  >0) fidelity=1
          if (index(line,'medium')>0) fidelity=2
          if (index(line,'coarse')>0) fidelity=3
        goto 2
      endif
    endif
  goto 1

9 continue

  rewind(3)

End Subroutine ReadGrid
