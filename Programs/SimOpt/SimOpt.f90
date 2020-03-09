program SimOpt
! version 3.0     July 18, 2019
  use Declaration
  implicit none

! Geometry and input variables
  integer                         :: n
  character(120)                  :: fnameinp,fname,scrname

! Optimization variables
  real*8,external                 :: func
  integer,parameter               :: MaxCycle=100
  real*8,parameter                :: Gtol = 1.0d-8
! real*8,parameter                :: Gtol = 1.0d-7
  real*8,parameter                :: hbar=6.582119514d-16
  integer                         :: nCycle
  real*8                          :: p(6)
  real*8                          :: fret,gNorm
  logical                         :: conv1,conv2
  logical                         :: ConvGeom(MaxCycle)

  integer                         :: i,j

  n = iargc()
  if (n /= 1) then
    write(6,'(a)')'usage: SimOpt input_file_name'
    call exit (8)
  endif
  call getarg(1,fnameinp)
  scrname = trim(fnameinp(1:len_trim(fnameinp)))
  n = scan(scrname,".",BACK=.true.)
  if (n > 0) fname = scrname(1:n)

! Open input file
  open(3,file=fnameinp,status='old',err=9998)

! Open geometry optimization file
  open(9,file=fname(1:len_trim(fname))//'molden',status='unknown')

! Read method
  call ReadMethod

! Geometry variables
  call GetConst

! Reduce LJ (hard-sphere) potential
  call ReadReduce

! Marcus theory
  if (method == 4) call GetMarcus
  
! Read number of atoms
  call ReadDim
  allocate (xA(3,nA), xB(3,nB),xBsave(3,nB))
  allocate (labA(nA), labB(nB))
  allocate (rWA(nA), rWB(nB))

! Read geometry
  call ReadInp
  xBsave = xB

! Read eAC,eCA
  call ReadE

! Convert atomic symbols to atomic numbers
  allocate (listA(nA), listB(nB))
  call List

! Get dimensions for basis variables
  call BasDim

! Allocate basis variables

! Allocate arrays for basis
!  List of centers for each NAO
!  Structure:
!          0    1    2
!        ----------------------
!  NAO=1 | n1 | c1 | c2 | cn1 |
!        | .. | .. | ........ |
!  NAO=m | nm | c1 | c2 | cnm |
!        ----------------------
!      n1, n1, ..., nm = number of elements

  allocate (OnAtmA(nNAOA,0:nOnAtmA),OnAtmB(nNAOB,0:nOnAtmB))

! (nNAO) -> nCont, nPrim
  allocate (DimnContA(nNAOA),DimnContB(nNAOB))
  allocate (DimnPrimA(nNAOA),DimnPrimB(nNAOB))
! (nNao,nCont) -> type, C, ng
  allocate (TypContA(nNAOA,maxnCA),TypContB(nNAOB,maxnCB))
  allocate (CContA(nNAOA,maxnCA,6),CContB(nNAOB,maxnCB,6))
  allocate (ngA(nNAOA,maxnCA),ngB(nNAOB,maxnCB))
! (nNao,nPrim) -> alpha, g
  allocate (alphaA(nNAOA,maxnPA),alphaB(nNAOB,maxnPB))
  allocate (gA(nNAOA,maxnPA),gB(nNAOB,maxnPB))

! Read Basis
  call Basis

! Write header of log
  write(6,'(a)')'========================================'
  write(6,'(a)')'           Program SimOpt'
  write(6,'(a)')'========================================'
  write(6,'(a)')'Title='//Title(1:len_trim(Title))
  if (method==1) write(6,'(a)')'Calculation of T^A**2 matrix elements'
  if (method==2) write(6,'(a)')'Calculation of Lennard-Jones potential'
  if (method==3) write(6,'(a)')'Calculation of mixed (hard-sphere repulsion - T^A**2) surface'
  if (method==4) write(6,'(a)')'Calculation of rate on the (hard-sphere repulsion - rate) surface'

  if (method>=3) write(6,'(a,es15.7)')'MixLJ =',MixLJ
  if (method==4) write(6,'(a, f9.3, a)')'Internal endoergicity is set to', endoerg, ' eV.'
  write(6,*)
  write(6,'(a)')'Geometry A:'
  do i=1,nA
    write(6,'(a3,3f12.7)')labA(i),(xA(j,i),j=1,3)
  enddo
  write(6,'(a)')'Geometry B:'
  do i=1,nA
    write(6,'(a3,3f12.7)')labB(i),(xB(j,i),j=1,3)
  enddo

! Get number of MOs
  call GetNMOs

  write(6,'(a)')'========================================'
  write(6,*)
  write(6,'(a)')'Mol A:'
  write(6,'(a,i4)')'Number of NAOs                 =',nNAOA
  do i=1,nNAOA
    write(6,'(a,i3,a,100i4)')'NAO:',i,' located on:',(OnAtmA(i,j),j=1,OnAtmA(i,0))
  enddo
  write(6,'(a)')'BASIS'
  call WriteBasis('a')

  write(6,*)
  write(6,'(a)')'Mol B:'
  write(6,'(a,i4)')'Number of NAOs                 =',nNAOB
  do i=1,nNAOB
    write(6,'(a,i3,a,100i4)')'NAO:',i,' located on:',(OnAtmB(i,j),j=1,OnAtmB(i,0))
  enddo
  write(6,'(a)')'BASIS'
  call WriteBasis('b')

! Number of AOs
  call AOdim

! alloacte temporary arrays
  allocate (moTMPA(nAOA,2),moTMPB(nAOB,2))
  allocate (SsaveA(nAOA,nAOA),SsaveB(nAOB,nAOB))
  allocate (gTMP(3,nB))

! Convert function types into integer indices
  allocate (a(nAOA,3), b(nAOB,3))
  call TypIndex

! set Hii for Wolfberg-Helmhols formula (ICON default values)
  allocate (HiiA(nAOA), HiiB(nAOB))
  call SetHii

! copy data into Primitive structures
  allocate (aPrimA(nPGA,3), alphaPrimA(nPGA), gPrimA(nPGA), CPrimA(nPGA))
  allocate (aPrimB(nPGB,3), alphaPrimB(nPGB), gPrimB(nPGB), CPrimB(nPGB))
  call Primitive
  allocate (nAOPrimA(nAOA), nAOPrimB(nAOB))
  allocate (nXA(nAOA), nXB(nAOB))
  call AOPrim
  allocate (AOsA(nAOA), AOsB(nAOB))
  call AOlabels

! Normalization of primitive functins
  allocate (GnormA(nPGA),GnormB(nPGB))
  call PGNorm

! Allocate overlap matrices
  allocate (SaoA(nAOA,nAOA), SaoB(nAOB,nAOB), SaoAB(nAOA,nAOB))

! AO overlap:

  call AOoverlap('aa')
  call AOoverlap('bb')

!!! Molecule A !!!
  allocate (CA(MOsA,2), MOA(nAOA,2))
  call ReadMOs('a')
  call ConstructMOs('a')
! Normalization A
  call Normalize('a')
! Orthogonalization A
  call Orthogonalize('a',moTMPA,SsaveA)

!!! Molecule B !!!
  allocate (CB(MOsB,2), MOB(nAOB,2), MOBsave(nAOB,2))
  call ReadMOs('b')
  call ConstructMOs('b')
! Normalize B
  call Normalize('b')
! Orthogonalize B
  call Orthogonalize('b',moTMPB,SsaveB)
! Save orbitals of B
  MOBsave = MOB

! Print geometry parameters
  write(6,'(a)')'Initial structure:'
  do i=1,6
    write(6,'(a,2x,f13.6)')gs(i),g(i)
  enddo
  write(6,*)

  write(6,'(10x,a,48x,a)')'==================== OPTIMIZATION ==================', &
                          ' Convergence: geom grad'
  ConvGeom = .false.
  nCycle=0
  p = g
  fret=func(p)
  write(6,'(i3,a,es15.7,a,6f12.6,a)')nCycle,'  F=',fret,'  geom=',p, &
           '          ----------------------'
  call FirstGuess(p)

  do i=1,MaxCycle
    call DFPmin(p,6,Gtol,fret,gNorm,conv1,conv2)
    call WriteMolden(fret)
    nCycle=nCycle+1
    ConvGeom(i)=conv1
    if (i>3) then
      if (ConvGeom(i-2).and.ConvGeom(i-1).and.ConvGeom(i)) conv2=.true.
    endif
    write(6,'(i3,a,es15.7,a,6f12.6,a,es14.6,l4,l5)')nCycle,'  F=', &
                 fret,'  geom=',p,'  gNorm=',gNorm,conv1,conv2
    if (conv1) last=.true.
    if (conv1.and.conv2) then
      write(6,'(a)')'     ---- Final Cycle ----'
      call DFPmin(p,6,Gtol,fret,gNorm,conv1,conv2)
      call WriteMolden(fret)
      conv2=.true.
      write(6,'(i3,a,es15.7,a,6f12.6,a,es14.6,l4,l5)')nCycle+1,'  F=', &
                   fret,'  geom=',p,'  gNorm=',gNorm,conv1,conv2
      write(6,'(a,i5,a)')'CONVERGED in',nCycle+1,'  cycles'
      goto 1
    else
    endif
  enddo

  write(6,'(a)')'NO CONVERGENCE'

1 continue

  write(6,*)
  write(6,'(a)')'Optimized structure:'
  do i=1,6
    write(6,'(a,2x,f13.6)')gs(i),g(i)
  enddo
  write(6,*)
  if (method==1) write(6,'(4(a,es15.7))')'TA^2=',SF
  if (method==1) write(6,'(4(a,es15.7))')'E-LJ=',SF
  if (Method==3) write(6,'(4(a,es15.7))')'TA^2=',SFa*SFa,'  E-Rep=',EnLJ,'  mix*Rep-TA^2=',SF,' (eV^2)'
  if (Method==4) write(6,'(4(a,es15.7))')'k=',k_rate,'  E-Rep=',EnLJ,'  mix*Rep-k=',SF,'  k(s^-1)=',k_rate/hbar

  goto 9999
9998 write(6,'(a,a)')'Error in the opening of the input file ',fnameinp(1:len_trim(fnameinp))
  call exit(8)

9999 continue

  write(6,*)
  write(6,'(a)')'Program SimOpt finished.'

  call Cpu()

End program SimOpt
