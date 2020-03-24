Subroutine ReadOctant(octant)
! version 1.7     May 12, 2017
  implicit none

  logical,intent(out)             :: octant

  character(120)                  :: line
  real*8                          :: x
  integer                         :: i
  logical                         :: OK

  rewind(3)

  octant=.false. ! By default do NOT sort and look for T Z > 0

  OK=.false.

1 read(3,'(a)',end=2)line
    if(line(1:1)=='#')goto 1
    if(len_trim(line)==0)goto 1
    call lowercase(line)
    if (OK) goto 4
    if((index(line,'method')>0).and.(index(line,':')>0)) OK=.true.
    goto 1
4   continue
      if( index(line,'octant')   >0)octant=.true. ! search for T Z > 0
    if (index(line,':')>0) goto 2
  goto 1

2 continue
  rewind(3)

End Subroutine ReadOctant
