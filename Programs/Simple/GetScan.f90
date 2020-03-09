Subroutine GetScan
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line

  Sini = 0.0d0; Sdelta = 0.0d0; nScan = 1

1 read(3,'(a)',end=9)line
    if (line(1:1) == '#') goto 1
    if (index(line,':') > 0) then
      call lowercase(line)
      if (index(line,'scan') > 0) then
2       read(3,'(a)',end=9)line
          if (line(1:1) == '#') goto 2
          if (index(line,':') > 0) goto 9
          if (line(1:3) == 'R X') call readScan(line(4:),Sini,Sdelta,nScan,1)
          if (line(1:3) == 'R Y') call readScan(line(4:),Sini,Sdelta,nScan,2)
          if (line(1:3) == 'R Z') call readScan(line(4:),Sini,Sdelta,nScan,3)
          if (line(1:3) == 'T X') call readScan(line(4:),Sini,Sdelta,nScan,4)
          if (line(1:3) == 'T Y') call readScan(line(4:),Sini,Sdelta,nScan,5)
          if (line(1:3) == 'T Z') call readScan(line(4:),Sini,Sdelta,nScan,6)
        goto 2
      endif
    endif
  goto 1

9 continue

  rewind(3)

End Subroutine GetScan
