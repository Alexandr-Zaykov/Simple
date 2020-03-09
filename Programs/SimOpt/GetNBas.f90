subroutine GetNBas(mol,n,lab,nNAO,nOnAtm,maxnC,maxnP)
! version 3.0     July 18, 2019
  implicit none

  character(1),intent(in)         :: mol
  integer,intent(in)              :: n
  character(2),intent(in)         :: lab(n)

  integer,intent(out)             :: nNAO
  integer,intent(out)             :: nOnAtm
  integer,intent(out)             :: maxnC,maxnP

  integer,external                :: nonblank,isblank
  logical,external                :: isnumber
  integer,external                :: ReadLine

  integer,allocatable             :: list(:)

  character(120)                  :: line
  integer                         :: l,i
  character(2)                    :: AtmSymbol
  integer                         :: nC,nP
  integer                         :: nlist
  logical                         :: OK
  integer                         :: x
  
  allocate (list(n))

  nNAO=0
  nOnAtm=0
  maxnC=0
  maxnP=0

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
    x=0
    do i=1,n
      if (AtmSymbol==lab(i)) then
        OK=.true.
        x=x+1
      endif
    enddo
    call ReadNPrim(nC,nP)
    if (OK) then
      nNAO=nNAO+1
      nOnAtm=max(nOnAtm,x)
      maxnC=max(maxnC,nC)
      maxnP=max(maxnP,nP)
    endif
  else
! Case: Atom numbers
    call ReadList(line,n,nC,nlist,list)
    call ReadNPrim(nC,nP)
    call lowercase(line)
    if (index(line,mol)>0) then
      nNAO=nNAO+1
      nOnAtm=max(nOnAtm,nlist)
      maxnC=max(maxnC,nC)
      maxnP=max(maxnP,nP)
    else
      nC=0
      nP=0
    endif
  endif

  goto 2


9999 write(6,'(a)')'Definition of NAO not found'
  call exit(8)

3 continue

  deallocate (list)

End subroutine GetNBas
