Subroutine ReadInp
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line
  integer                         :: i,j,k,l
  integer,parameter               :: nElm=18
  character(2)                    :: atmlab(nElm)
  real*8                          :: rW(nElm)
  character(3)                    :: tmp

! Paulig ...
! Bondi ...
! S.S.Batsanov, Inorganic Materials, Vol. 37, No. 9, 2001, pp. 871-885
  data atmlab /'h ',                              'he',  &
               'li','be','b ','c ','n ','o ','f ','ne',  &
               'na','mg','al','si','p ','s ','cl','ar'/
  data rW     /1.20,                              1.40,  &
               1.90,1.38,1.20,1.70,1.50,1.40,1.35,1.54,  &
               2.32,1.96,1.75,1.68,1.63,1.85,1.80,1.88/

1 read(3,'(a)',end=8)line
    if (line(1:1) == '#') goto 1
    call lowercase(line)
    if (index(line,':') > 0) then
      if (index(line,'title') > 0) then
        read(3,'(a)')Title
      endif
      if (index(line,'geoma') > 0) then
        call GetX(labA,xA,nA)
      endif
      if (index(line,'geomb') > 0) then
        call GetX(labB,xB,nB)
      endif
    endif
  goto 1

8 continue

  rewind(3)

  do i=1,nA
    tmp=labA(i)
    call lowercase(tmp)
    rWA(i)=1.8d0 ! Generic value
    do j=1,nElm
      if (index(tmp,atmlab(j)) > 0) rWA(i)=rW(j)
    enddo
  enddo
  do i=1,nB
    tmp=labB(i)
    call lowercase(tmp)
    rWB(i)=1.8d0 ! Generic value
    do j=1,nElm
      if (index(tmp,atmlab(j)) > 0) rWB(i)=rW(j)
    enddo
  enddo

  return

End Subroutine ReadInp
