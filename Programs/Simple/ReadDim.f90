Subroutine ReadDim
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line

  nA=0
  nB=0

1 read(3,'(a)',end=8)line  
    if (line(1:1) == '#') goto 1
    call lowercase(line)
    if (index(line,':') > 0) then
      if (index(line,'geoma') > 0) then
2       read(3,'(a)',end=8)line
          if (line(1:1) == '#') goto 2
          if (index(line,':') <= 0) then
            nA=nA+1
            goto 2
          else
            backspace(3)
            goto 1
          endif
      endif
      if (index(line,'geomb') > 0) then
3       read(3,'(a)',end=8)line
          if (line(1:1) == '#') goto 3
          if (index(line,':') <= 0) then
            nB=nB+1
            goto 3 
          else
            backspace(3)
            goto 1 
          endif
      endif
    endif
  goto 1

8 continue

  rewind(3)

  if (nA == 0) then
    write(6,'(a)')'No coordinates found for A'
    call exit(8)
  endif
  if (nB == 0) then
    write(6,'(a)')'No coordinates found for B'
    call exit(8)
  endif

  return

End Subroutine ReadDim
