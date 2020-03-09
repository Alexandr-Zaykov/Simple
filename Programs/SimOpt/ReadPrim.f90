Subroutine ReadPrim(maxnC,maxnP,nC,nP,TypCont,CCont,ng,alpha,g)
! version 3.0     July 18, 2019
  implicit none

  integer,intent(in)              :: maxnC,maxnP
  integer,intent(in)              :: nC

  integer,intent(out)             :: nP
  character(1),intent(out)        :: TypCont(maxnC)
  real*8,intent(out)              :: CCont(maxnC,6)
  integer,intent(out)             :: ng(maxnC)
  real*8,intent(out)              :: alpha(maxnP)
  real*8,intent(out)              :: g(maxnP)

  integer,external                :: ReadLine

  character(120)                  :: line
  character(1)                    :: tmp
  integer                         :: i,j,k
  integer                         :: n

  n=0
  nP=0

  do i=1,nC
    j=ReadLine(line)
    read(line,*)tmp
    call lowercase(tmp)
    if (tmp=='s') read(line,*)TypCont(i),ng(i),CCont(i,1)
    if (tmp=='p') read(line,*)TypCont(i),ng(i),(CCont(i,j),j=1,3)
    if (tmp=='d') read(line,*)TypCont(i),ng(i),(CCont(i,j),j=1,6)
    nP=nP+ng(i)
    do j=1,ng(i)
      k=ReadLine(line)
      n=n+1
      read(line,*)alpha(n),g(n)
    enddo
  enddo

  return

End Subroutine ReadPrim
