Subroutine ReadMethod
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(120)                  :: line
  real*8                          :: x
  integer                         :: i
  logical                         :: OK

  rewind(3)

  method=1      ! default method
  MixLJ=1.0d0   ! Mixing of LJ and SF

  OK=.false.

1 read(3,'(a)',end=2)line
    if(line(1:1)=='#')goto 1
    if(len_trim(line)==0)goto 1
    call lowercase(line)
    if (OK) goto 5
    if((index(line,'method')>0).and.(index(line,':')>0)) OK=.true.
    goto 1
5   continue
    if( index(line,'sf')   >0)method=1 ! Calculate SF matrix elements
    if( index(line,'sf2')  >0)method=2 ! Calculate SF**2 matrix elements
    if( index(line,'lj')   >0)method=3 ! Calculate Lennard-Jones 6-12 potential only
    if( index(line,'sf-lj')>0)method=4 ! Calculate MixLJ*LJ - SF**2  (For method >3, repulsion = hard-sphere)
    if( index(line,'rate') >0)method=5 ! Calculate rate and lifetime in Marcus approximation
    if( index(line,'sf-lj')>0) then
      i=index(line,'sf-lj')+5
      read(line(i:),*,end=3)x
        MixLJ=x
3     continue
    endif
    if( index(line,'rate')>0) then
      i=index(line,'rate')+4
      read(line(i:),*,end=4)x
        MixLJ=x
4     continue
    endif
    if (index(line,':')>0) goto 2
  goto 1

2 continue
  rewind(3)

End Subroutine ReadMethod
