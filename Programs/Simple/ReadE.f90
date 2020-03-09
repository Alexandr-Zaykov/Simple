Subroutine ReadE
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line
  integer                         :: n
  real*8                          :: conv


  rewind(3)

1 read(3,'(a)',end=99)line
    if (line(1:1) == '#') goto 1
    if (index(line,':') > 0) then
      call lowercase(line)
      conv=1.0d0
      if (index(line,'ev'  )>0) conv=1.0d0
      if (index(line,'mev' )>0) conv=1.0d-3
      if (index(line,'cal' )>0) conv=4.33641153d-5
      if (index(line,'kcal')>0) conv=4.33641153d-2
      if (index(line,'eac') > 0) then
        n=index(line,':')+1
        read(line(n:),*)eAC
        eAC=eAC*conv
      endif
      if (index(line,'eca') > 0) then
        n=index(line,':')+1
        read(line(n:),*)eCA
        eCA=eCA*conv
      endif
      if (index(line,'de(ct)') > 0) then
        n=index(line,':')+1
        read(line(n:),*)eCT
        eCT=eCT*conv
      endif
      goto 1
    endif
  goto 1

99 continue

  if ((eAC==0.0d0).and.(eCA==0.0d0)) then
    if (eCT==0.0d0) then
      stop 'Energy of charge-transfer states not spcified'
    else
     eAC=eCT; eCA=eCT
    endif
  else
    if (eAC==0.0d0) stop 'Energy of AC charge-transfer states not spcified'
    if (eCA==0.0d0) stop 'Energy of CA charge-transfer states not spcified'
  endif
  if (eCT/=0.0d0) then
    if ((eCT/=eAC).and.(eCT/=eCA)) then
      stop 'Conflict between dE(CT) and eAC and/or eCA values'
    endif
  endif

  rewind(3)

  return

End Subroutine ReadE
