Subroutine GetMarcus
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line
  real*8                          :: conv
  integer                         :: n

! Default reorganization energy (approximately OK for most of the structures)
  lambda_marcus=0.2d0
  endoerg=0.0d0
  conv=1.0d0

  rewind(3)

1 read(3,'(a)',end=99)line
    if (line(1:1) == '#') goto 1
    if (index(line,':') > 0) then
      if (index(line,'lambda') > 0) then
        n=index(line,':')+1
        read(line(n:),*)lambda_marcus
        !Default unit is eV
        if (index(line,'eV'  )>0) conv=1.0d0
        if (index(line,'meV' )>0) conv=1.0d-3
        if (index(line,'cal' )>0) conv=4.33641153d-5
        if (index(line,'kcal')>0) conv=4.33641153d-2
        lambda_marcus=lambda_marcus*conv
        goto 1
      endif
      conv=1.00d0
      if (index(line,'dE(int)') > 0) then
        n=index(line,':')+1
        read(line(n:),*)endoerg
        !Default unit is eV
        if (index(line,'eV'  )>0) conv=1.0d0
        if (index(line,'meV' )>0) conv=1.0d-3
        if (index(line,'cal' )>0) conv=4.33641153d-5
        if (index(line,'kcal')>0) conv=4.33641153d-2
        endoerg=endoerg*conv
        goto 1
      endif
    endif
  goto 1

99 continue

  rewind(3)
  return
End Subroutine GetMarcus
