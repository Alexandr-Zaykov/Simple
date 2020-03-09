Subroutine GetConst
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line
  logical                         :: notOK(6),OK
  integer                         :: i

  g  = 0.0d0
  notOK = .true.

1 read(3,'(a)',end=9)line
    if (line(1:1) == '#') goto 1
    if (index(line,':') > 0) then
      call lowercase(line)
      if (index(line,'const') > 0) then
2       read(3,'(a)',end=9)line
          if (line(1:1) == '#') goto 2
          if (index(line,':') > 0) goto 9
          if (line(1:3) == gs(1)) then
            read(line(4:),*)g(1)
            notOK(1)=.false.
          endif
          if (line(1:3) == gs(2)) then
            read(line(4:),*)g(2)
            notOK(2)=.false.
          endif
          if (line(1:3) == gs(3)) then
            read(line(4:),*)g(3)
            notOK(3)=.false.
          endif
          if (line(1:3) == gs(4)) then
            read(line(4:),*)g(4)
            notOK(4)=.false.
          endif
          if (line(1:3) == gs(5)) then
            read(line(4:),*)g(5)
            notOK(5)=.false.
          endif
          if (line(1:3) == gs(6)) then
            read(line(4:),*)g(6)
            notOK(6)=.false.
          endif
        goto 2
      endif
    endif
  goto 1

9 continue

  OK=.false.
  do i=1,6
    if (notOK(i)) then
      write(6,'(3a)')'Definition ',gs(i),' missing in Const section'
      OK=.true.
    endif
  enddo
  if (OK)then
    write(6,'(a)')'Program stops'
    call exit(8)
  endif

  rewind(3)

End Subroutine GetConst
