subroutine ReadBas(mol,n,lab,nNAO,maxnC,maxnP,DimnCont,DimnPrim,TypCont,CCont,ng,alpha,g)
! version 3.0     July 18, 2019
  implicit none

  character(1),intent(in)         :: mol
  integer,intent(in)              :: n,nNAO,maxnC,maxnP
  character(2),intent(in)         :: lab(n)

  integer,intent(out)             :: DimnCont(nNAO),DimnPrim(nNAO)
  character(1),intent(out)        :: TypCont(nNAO,maxnC)
  real*8,intent(out)              :: CCont(nNAO,maxnC,6)
  integer,intent(out)             :: ng(nNAO,maxnC)
  real*8,intent(out)              :: alpha(nNAO,maxnP)
  real*8,intent(out)              :: g(nNAO,maxnP)

  integer,external                :: nonblank
  logical,external                :: isnumber
  integer,external                :: ReadLine

  character(1),allocatable        :: tmpTypCont(:)
  real*8,allocatable              :: tmpCCont(:,:)
  integer,allocatable             :: tmpng(:)
  real*8,allocatable              :: tmpalpha(:),tmpg(:)
  integer,allocatable             :: list(:)

  character(120)                  :: line
  integer                         :: l,i,j,k,m
  character(2)                    :: AtmSymbol
  integer                         :: nC
  integer                         :: nP
  integer                         :: nN
  integer                         :: nlist
  logical                         :: OK
  character(2)                    :: tmplab
  
  allocate (tmpTypCont(maxnC), tmpCCont(maxnC,6), tmpng(maxnC), tmpalpha(maxnP), tmpg(maxnP))
  allocate (list(n))

  nN=0

! Find the NAO section in input

1 i=ReadLine(line)
  if (i==1) goto 9999
  if (i==2) then
    call lowercase(line)
    if (index(line,'nao')>0) goto 2
  endif
  goto 1

! Read NAO

2 continue
  i=ReadLine(line) ! Read the NAO definition line (e.g. C 4, or 1-3 A 1)
  if (i==2) goto 3

  l=1
  i=nonblank(l,line)

  if (.not. isnumber(line(i:i))) then
! Case: Atom label
    OK=.false.
    nP=0
    read(line,*) AtmSymbol,nC
    call lowercase(AtmSymbol)
    do i=1,n
      tmplab=lab(i)
      call lowercase(tmplab)
      if (AtmSymbol==tmplab) then
        OK=.true.
      endif
    enddo
    call ReadPrim(maxnC,maxnP,nC,nP,tmpTypCont,tmpCCont,tmpng,tmpalpha,tmpg)
    if (OK) then
      nN=nN+1
      DimnCont(nN)=nC
      DimnPrim(nN)=nP
      do j=1,nC
        TypCont(nN,j)=tmpTypCont(j)
        call lowercase(tmpTypCont(j))
        if (tmpTypCont(j)=='s') k=1
        if (tmpTypCont(j)=='p') k=3
        if (tmpTypCont(j)=='d') k=6
        do m=1,k
          CCont(nN,j,m)=tmpCCont(j,m)
        enddo
        ng(nN,j)=tmpng(j)
      enddo
      do j=1,nP
        alpha(nN,j)=tmpalpha(j)
        g(nN,j)=tmpg(j)
      enddo
    endif
  else
! Case: Atom numbers
    call ReadList(line,n,nC,nlist,list)
    if (index(line,mol)>0) then
!     Read of primitives
      call ReadPrim(maxnC,maxnP,nC,nP,tmpTypCont,tmpCCont,tmpng,tmpalpha,tmpg)
      nN=nN+1
      DimnCont(nN)=nC
      DimnPrim(nN)=nP
      do j=1,nC
        TypCont(nN,j)=tmpTypCont(j)
        call lowercase(tmpTypCont(j))
        if (tmpTypCont(j)=='s') k=1
        if (tmpTypCont(j)=='p') k=3
        if (tmpTypCont(j)=='d') k=6
        do m=1,k
          CCont(nN,j,m)=tmpCCont(j,m)
        enddo
        ng(nN,j)=tmpng(j)
      enddo
      do j=1,nP
        alpha(nN,j)=tmpalpha(j)
        g(nN,j)=tmpg(j)
      enddo
    else
!     Blank read of primitives
      call ReadPrim(maxnC,maxnP,nC,nP,tmpTypCont,tmpCCont,tmpng,tmpalpha,tmpg)
    endif
  endif
  goto 2

9999 write(6,'(a)')'Definition of NAO not found'
  call exit(8)

3 continue

  deallocate (tmpTypCont, tmpCCont, tmpng, tmpalpha, tmpg)
  deallocate (list)

End subroutine ReadBas
