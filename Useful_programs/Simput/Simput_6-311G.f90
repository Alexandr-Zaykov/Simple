!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                               !!!
!!!     PROGRAM: Simput                 Version: 1                !!!
!!!     written by:                     Eric Buchanan             !!!
!!!     6-311G and cutoff minipatch by: Alexandr "Sasa" Zaykov    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM Simput
  ! This program generates input files for the Simple
  ! program from NBO6 output files.
  ! It requires 4 NBO files for monomer A:
  ! FILE.31 : Coordinates
  ! FILE.33 : NAO vectors in AO basis
  ! FILE.46 : Basis set info
  ! FILE.49 : MO vectors in NAO basis
  ! A file (Basis.simput) with the basis is required
  ! in Turbomole format from EMSL (6-311+G for B,C,N,O,F right now)
  ! Only systems with HOMO/LUMO comprised of 2P valence orbitals right now

  IMPLICIT NONE

  character(50)               :: inputForSimpleInp
  character(50)               :: DIR_A,DIR_B,BASIS
  character(50)               :: title
  character(10)               :: method,grid,planarity
  integer                     :: numProc
  character(85)               :: fileline
  integer                     :: IOs,i,j,k,l,m,n
  integer                     :: numA_A,numAO_A, numNAO_A, numMO_A,numVAL_A
  integer,allocatable         :: AO_A(:,:),NAO_A(:,:),Valence_A(:,:)
  integer                     :: numA_B,numAO_B, numNAO_B, numMO_B,numVAL_B
  integer,allocatable         :: AO_B(:,:),NAO_B(:,:),Valence_B(:,:)
  integer,allocatable         :: PH_A(:,:,:),PH_B(:,:,:),buffer(:,:)
  character(3)                :: A1,A2,A3,A4,A5,A6,A7
  integer                     :: I1,I2,I3,I4,I5,I6,I7
  integer                     :: I1a,I2a,I3a,I4a,I5a,I6a,I7a
  real*8                      :: F1,F2,F3,F4,F5,tx,ty,tz,rx,ry,rz,eac,eca
  real*8,allocatable          :: AONAOm_A(:,:),NAOMOm_A(:,:),COORDS_A(:,:)
  real*8,allocatable          :: AONAOm_B(:,:),NAOMOm_B(:,:),COORDS_B(:,:)
  real*8,allocatable          :: BASISm(:,:)
  integer                     :: nHOMO_A,nLUMO_A,nHOMO_B,nLUMO_B
! Cutoff minipatch
  integer                     :: j_trunc
  real*8                      :: cutoff_val
! 2.0 Fidelity Patch
  real*8                      :: E_nine, Reduce, Lambda
  logical                     :: noljwarn
  ! Open input for Simput.
  ! This file should look as follows:

  ! HOMO A:71
  ! LUMO A:72
  ! HOMO B:71
  ! LUMO B:72
  ! DIR A:M1/
  ! DIR B:M1/
  ! Basis Set:<6-311+G>.basis
  ! Title:<A1B1>
  ! Method:SF-LJ 1.0
  ! Grid:<coarse,medium,fine>
  ! T X:
  ! T Y:
  ! T Z:
  ! R X:
  ! R Y:
  ! R Z:
  ! EAC:
  ! ECA:
  ! Pz Orbitals Only:<yes/no>
  ! Processors:
  ! dE(CT):
  ! Cutoff value:
  ! Reduce:
  ! Lambda:
  ! noLJwarn (Optional!)
  ! %end

  i = iargc()
  if (i /= 1) then
    write(6,'(A)')'Usage: Simput input_file_name'
    call exit (8)
  endif
  call getarg(1,inputForSimpleInp)

  open(unit=8,file=inputForSimpleInp,status='old',IOStat=IOs)
  if (IOs /= 0) then
     write(6,'(A,i1)')'Error number', IOs
     write(6,'(A)')'Failure of input opening'
     stop
  end if

  !Init
  cutoff_val=0.0d0
  Reduce=1.0d0
  noLJwarn=.false.
  ! Read the Simput input file:
  do
    read (8,'(A)', IOstat=IOs) fileline
    if (IOs /= 0) then
      write(6,'(A)')'Failure reading simput.inp'
      write(6,'(A)')'Error number',IOs
      stop
    end if
    if (fileline(1:7) == "HOMO A:") then
      read (fileline(8:),'(i4)') nHOMO_A
    elseif (fileline(1:7) == "LUMO A:") then
      read (fileline(8:),'(i4)') nLUMO_A
    elseif (fileline(1:7) == "HOMO B:") then
      read (fileline(8:),'(i4)') nHOMO_B
    elseif (fileline(1:7) == "LUMO B:") then
      read (fileline(8:),'(i4)') nLUMO_B
    elseif (fileline(1:6) == "DIR A:") then
      read (fileline(7:),'(A)') DIR_A
    elseif (fileline(1:6) == "DIR B:") then
      read (fileline(7:),'(A)') DIR_B
    elseif (fileline(1:10) == "Basis Set:") then
      read (fileline(11:),'(A)') BASIS
    elseif (fileline(1:6) == "Title:") then
      read (fileline(7:),'(A)') title
    elseif (fileline(1:7) == "Method:") then
      read (fileline(8:),'(A)') method
    elseif (fileline(1:5) == "Grid:") then
      read (fileline(6:),'(A)') grid
    elseif (fileline(1:4) == "T X:") then
      read (fileline(5:),'(F12.7)') tx
    elseif (fileline(1:4) == "T Y:") then
      read (fileline(5:),'(F12.7)') ty
    elseif (fileline(1:4) == "T Z:") then
      read (fileline(5:),'(F12.7)') tz
    elseif (fileline(1:4) == "R X:") then
      read (fileline(5:),'(F12.7)') rx
    elseif (fileline(1:4) == "R Y:") then
      read (fileline(5:),'(F12.7)') ry
    elseif (fileline(1:4) == "R Z:") then
      read (fileline(5:),'(F12.7)') rz
    elseif (fileline(1:4) == "EAC:") then
      read (fileline(5:),'(F12.7)') eac
    elseif (fileline(1:4) == "ECA:") then
      read (fileline(5:),'(F12.7)') eca
    elseif (fileline(1:7) == "Planar:") then
      read (fileline(8:),'(A)') planarity
    elseif (fileline(1:11) == "Processors:") then
      read (fileline(12:16),'(i4)') numProc
    elseif (fileline(1:7) == "dE(CT):") then
      read (fileline(8:),'(F12.7)') E_nine
! Read cutoff, default = 0.0d0
    elseif (fileline(1:13) == "Cutoff value:") then
      read (fileline(14:),'(F12.7)') cutoff_val
    elseif (fileline(1:7) == "Reduce:") then
      read (fileline(8:),'(F12.7)') Reduce
    elseif (fileline(1:7) == "Lambda:") then
      read (fileline(8:),'(F12.7)') Lambda
    elseif (fileline(1:8) == "noLJwarn") then
      noljwarn=.true.
! End of cutoff minipatch read
    elseif (fileline(1:4) == "%end") then
      exit
    endif
  enddo
  close(8)


  ! Open FILE.31 file for A
  open (unit=13,file=trim(DIR_A)//'FILE.31',status='old',IOStat=IOs)
  if (IOs /= 0) then
    write(6,'(A)')'There is no FILE.31/Coordinate info file for monomer A'
    stop
  end if
  ! Read Atoms & Coordinates
  do
    read(13,'(A)', IOstat=IOs) fileline
    if (fileline(1:22) == ' Basis set information') then
      read(13,'(/,4X,I3)') numA_A
      exit
    endif
  enddo
  read(13,'(A)')
  ! COORDS_A(i,1)=ith Atom type
  ! COORDS_A(i,2)=ith Atom X Coord
  ! COORDS_A(i,3)=ith Atom Y Coord
  ! COORDS_A(i,4)=ith Atom Z Coord
  allocate(COORDS_A(numA_A,4))

  do i=1,numA_A
    read(13,'(3X,I2,2X,F12.9,2X,F12.9,2X,F12.9)')I1,F1,F2,F3
    COORDS_A(i,1)=I1
    COORDS_A(i,2)=F1
    COORDS_A(i,3)=F2
    COORDS_A(i,4)=F3
  enddo
  close(13)
  ! Open FILE.31 file for B
  open (unit=18,file=trim(DIR_B)//'FILE.31',status='old',IOStat=IOs)
  if (IOs /= 0) then
    write(6,'(A)')'There is no FILE.31/Coordinate info file for monomer B'
    stop
  end if

  ! Read Atoms & Coordinates
  do
    read(18,'(A)', IOstat=IOs) fileline
    if (fileline(1:22) == ' Basis set information') then
      read(18,'(/,4X,I3)') numA_B
      exit
    endif
  enddo
  read(18,'(A)')
  ! COORDS_B(i,1)=ith Atom type
  ! COORDS_B(i,2)=ith Atom X Coord
  ! COORDS_B(i,3)=ith Atom Y Coord
  ! COORDS_B(i,4)=ith Atom Z Coord
  allocate(COORDS_B(numA_B,4))

  do i=1,numA_B
    read(18,'(3X,I2,2X,F12.9,2X,F12.9,2X,F12.9)')I1,F1,F2,F3
    COORDS_B(i,1)=I1
    COORDS_B(i,2)=F1
    COORDS_B(i,3)=F2
    COORDS_B(i,4)=F3
  enddo
  close(18)

  ! Open simput.basis file
  open (unit=14,file=BASIS,status='old',IOStat=IOs)
  if (IOs /= 0) then
    write(6,'(A)')'There is no simput.basis/Basis set file'
    stop
  end if

  ! Read Basis set (6-311+G only now)
  ! BASISm(i,1)=ith Atom type
  ! BASISm(i,2-13)=ith Atoms 6 sets of 
  ! coefficients in pairs by shell (4 P shells)
  ! B=5,C=1,N=2,O=3,F=4 only now
  ! (Most of this should be redone and is crap)
  allocate(BASISm(5,11))
  ! Find B Atom
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(1:1) == 'b') then
      exit
    endif
  enddo
  ! Find P orbitals of B
  read(14,'(A)')
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(5:8) == '3  p') then
      exit
    endif
  enddo
  BASISm(5,1)=5
  do i=1,3
    read(14,'(4X,F11.7,13X,F11.7)')F1,F2
    BASISm(5,(i-1)*2+2)=F1
    BASISm(5,(i-1)*2+3)=F2
  enddo
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(5,8)=F1
  BASISm(5,9)=F2
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(5,10)=F1
  BASISm(5,11)=F2
!  read(14,'(A)')
! Diffuse functions
!  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
!  BASISm(5,12)=F1
!  BASISm(5,13)=F2
  ! Find C Atom
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(1:1) == 'c') then
      exit
    endif
  enddo
  ! Find P orbitals of C
  read(14,'(A)')
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(5:8) == '3  p') then
      exit
    endif
  enddo
  BASISm(1,1)=6
  do i=1,3
    read(14,'(4X,F11.7,13X,F11.7)')F1,F2
    BASISm(1,(i-1)*2+2)=F1
    BASISm(1,(i-1)*2+3)=F2
  enddo
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(1,8)=F1
  BASISm(1,9)=F2
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(1,10)=F1
  BASISm(1,11)=F2
!  read(14,'(A)')
! Diffuse functions
!  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
!  BASISm(1,12)=F1
!  BASISm(1,13)=F2
  ! Find N Atom
  rewind(14)
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(1:1) == 'n') then
      exit
    endif
  enddo
  ! Find P orbitals of N
  read(14,'(A)')
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(5:8) == '3  p') then
      exit
    endif
  enddo
  BASISm(2,1)=7
  do i=1,3
    read(14,'(4X,F11.7,13X,F11.7)')F1,F2
    BASISm(2,(i-1)*2+2)=F1
    BASISm(2,(i-1)*2+3)=F2
  enddo
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(2,8)=F1
  BASISm(2,9)=F2
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(2,10)=F1
  BASISm(2,11)=F2
!  read(14,'(A)')
! Diffuse functions
!  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
!  BASISm(2,12)=F1
!  BASISm(2,13)=F2
  ! Find O Atom
  rewind(14)
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(1:1) == 'o') then
      exit
    endif
  enddo
  ! Find P orbitals of O
  read(14,'(A)')
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(5:8) == '3  p') then
      exit
    endif
  enddo
  BASISm(3,1)=8
  do i=1,3
    read(14,'(4X,F11.7,13X,F11.7)')F1,F2
    BASISm(3,(i-1)*2+2)=F1
    BASISm(3,(i-1)*2+3)=F2
  enddo
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(3,8)=F1
  BASISm(3,9)=F2
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(3,10)=F1
  BASISm(3,11)=F2
!  read(14,'(A)')
!  Diffuse functions
!  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
!  BASISm(3,12)=F1
!  BASISm(3,13)=F2
  ! Find F Atom
  rewind(14)
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(1:1) == 'f') then
      exit
    endif
  enddo
  ! Find P orbitals of F
  read(14,'(A)')
  do
    read(14,'(A)', IOstat=IOs) fileline
    if (fileline(5:8) == '3  p') then
      exit
    endif
  enddo
  BASISm(4,1)=9
  do i=1,3
    read(14,'(4X,F11.7,13X,F11.7)')F1,F2
    BASISm(4,(i-1)*2+2)=F1
    BASISm(4,(i-1)*2+3)=F2
  enddo
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(4,8)=F1
  BASISm(4,9)=F2
  read(14,'(A)')
  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
  BASISm(4,10)=F1
  BASISm(4,11)=F2
!  read(14,'(A)')
!  Diffuse functions
!  read(14,'(4X,F11.7,13X,F11.7)')F1,F2
!  BASISm(4,12)=F1
!  BASISm(4,13)=F2


  close(14)
  ! Open FILE.46 file for A
  open (unit=9,file=trim(DIR_A)//'FILE.46',status='old',IOStat=IOs)
  if (IOs /= 0) then
    write(6,'(A)')'There is no FILE.46/Basis info file for monomer A'
    stop
  end if

  ! Read number of AOs, NAOs, & MOs in monomer A
  do
    read(9,'(A)', IOstat=IOs) fileline
    if (fileline(3:4) == 'AO') then
      read(fileline(12:),'(i4)') numAO_A
      exit
    endif
  enddo
  rewind(9)
  do
    read(9,'(A)', IOstat=IOs) fileline
    if (fileline(3:5) == 'NAO') then
      read(fileline(12:),'(i4)') numNAO_A
      exit
    endif
  enddo

  numMO_A=numNAO_A
  rewind(9)
  allocate(AO_A(numAO_A,2),NAO_A(numNAO_A,3))
  AO_A=0
  NAO_A=0

  ! Read and store AO basis for A
  ! s=1,pz=2,px=3,py=4
  ! AO_A(i,1)=ith AO Atom
  ! AO_A(i,2)=ith AO Type
  if (mod(numAO_A,7)==0) then
    k=floor(numAO_A/7.0)
  else
    k=floor(numAO_A/7.0)+1
  endif
  i=0
  read(9,'(A)') fileline
  do j=1,k
    read(9,'(3X,I2,2X,A2,4X,I2,2X,A2,4X,I2,2X,A2,4X,I2,2X,A2,                           &
    & 4X,I2,2X,A2,4X,I2,2X,A2,4X,I2,2X,A2)')I1,A1,I2,A2,I3,A3,I4,A4,I5,A5,I6,A6,I7,A7
    if (A1=='s)') then
      i=i+1
      AO_A(i,1)=I1
      AO_A(i,2)=1
    elseif (A1=='pz') then
      i=i+1
      AO_A(i,1)=I1
      AO_A(i,2)=2
    elseif (A1=='px') then
      i=i+1
      AO_A(i,1)=I1
      AO_A(i,2)=3
    elseif (A1=='py') then
      i=i+1
      AO_A(i,1)=I1
      AO_A(i,2)=4
    elseif (A1/='s '.and.A1/='pz'.and.A1/='px'.and.A1/='py') then
      exit
    endif
    if (A2=='s)') then
      i=i+1
      AO_A(i,1)=I2
      AO_A(i,2)=1
    elseif (A2=='pz') then
      i=i+1
      AO_A(i,1)=I2
      AO_A(i,2)=2
    elseif (A2=='px') then
      i=i+1
      AO_A(i,1)=I2
      AO_A(i,2)=3
    elseif (A2=='py') then
      i=i+1
      AO_A(i,1)=I2
      AO_A(i,2)=4
    elseif (A2/='s '.and.A2/='pz'.and.A2/='px'.and.A2/='py') then
      exit
    endif
    if (A3=='s)') then
      i=i+1
      AO_A(i,1)=I3
      AO_A(i,2)=1
    elseif (A3=='pz') then
      i=i+1
      AO_A(i,1)=I3
      AO_A(i,2)=2
    elseif (A3=='px') then
      i=i+1
      AO_A(i,1)=I3
      AO_A(i,2)=3
    elseif (A3=='py') then
      i=i+1
      AO_A(i,1)=I3
      AO_A(i,2)=4
    elseif (A3/='s '.and.A3/='pz'.and.A3/='px'.and.A3/='py') then
      exit
    endif
    if (A4=='s)') then
      i=i+1
      AO_A(i,1)=I4
      AO_A(i,2)=1
    elseif (A4=='pz') then
      i=i+1
      AO_A(i,1)=I4
      AO_A(i,2)=2
    elseif (A4=='px') then
      i=i+1
      AO_A(i,1)=I4
      AO_A(i,2)=3
    elseif (A4=='py') then
      i=i+1
      AO_A(i,1)=I4
      AO_A(i,2)=4
    elseif (A4/='s '.and.A4/='pz'.and.A4/='px'.and.A4/='py') then
      exit
    endif
    if (A5=='s)') then
      i=i+1
      AO_A(i,1)=I5
      AO_A(i,2)=1
    elseif (A5=='pz') then
      i=i+1
      AO_A(i,1)=I5
      AO_A(i,2)=2
    elseif (A5=='px') then
      i=i+1
      AO_A(i,1)=I5
      AO_A(i,2)=3
    elseif (A5=='py') then
      i=i+1
      AO_A(i,1)=I5
      AO_A(i,2)=4
    elseif (A5/='s '.and.A5/='pz'.and.A5/='px'.and.A5/='py') then
      exit
    endif
    if (A6=='s)') then
      i=i+1
      AO_A(i,1)=I6
      AO_A(i,2)=1
    elseif (A6=='pz') then
      i=i+1
      AO_A(i,1)=I6
      AO_A(i,2)=2
    elseif (A6=='px') then
      i=i+1
      AO_A(i,1)=I6
      AO_A(i,2)=3
    elseif (A6=='py') then
      i=i+1
      AO_A(i,1)=I6
      AO_A(i,2)=4
    elseif (A6/='s '.and.A6/='pz'.and.A6/='px'.and.A6/='py') then
      exit
    endif
    if (A7=='s)') then
      i=i+1
      AO_A(i,1)=I7
      AO_A(i,2)=1
    elseif (A7=='pz') then
      i=i+1
      AO_A(i,1)=I7
      AO_A(i,2)=2
    elseif (A7=='px') then
      i=i+1
      AO_A(i,1)=I7
      AO_A(i,2)=3
    elseif (A7=='py') then
      i=i+1
      AO_A(i,1)=I7
      AO_A(i,2)=4
    elseif (A7/='s '.and.A7/='pz'.and.A7/='px'.and.A7/='py') then
      exit
    endif
  enddo

  ! Read and store NAO basis for A
  ! s=1,pz=2,px=3,py=4
  ! NAO_A(i,1)=ith NAO Atom
  ! NAO_A(i,2)=ith NAO Type
  ! NAO_A(i,3)=ith NAO Shell
  if (mod(numNAO_A,7)==0) then
    k=floor(numNAO_A/7.0)
  else
    k=floor(numNAO_A/7.0)+1
  endif
  i=0
  ! Positioning
  ! Planarity if (NAO_A(i,3)==2.and.NAO_A(i,2)==2) then
  read(9,'(A)') fileline
  do j=1,k
    read(9,'(3X,I2,2X,I1,A2,3X,I2,2X,I1,A2,3X,I2,2X,I1,A2,3X,I2,2X,I1,A2,               &
         &3X,I2,2X,I1,A2,3X,I2,2X,I1,A2,3X,I2,2X,I1,A2)')I1,I1a,A1,I2,I2a,A2,I3,        &
         &I3a,A3,I4,I4a,A4,I5,I5a,A5,I6,I6a,A6,I7,I7a,A7
 ! Read:  C 1( 1s ) C 1( 2s ) C 1( 3s ) C 1( 4s ) C 1( 2px) C 1( 3px) C 1( 4px)
    if (A1=='s ') then
      i=i+1
      NAO_A(i,1)=I1
      NAO_A(i,2)=1
      NAO_A(i,3)=I1a
    elseif (A1=='pz') then
      i=i+1
      NAO_A(i,1)=I1
      NAO_A(i,2)=2
      NAO_A(i,3)=I1a
    elseif (A1=='px') then
      i=i+1
      NAO_A(i,1)=I1
      NAO_A(i,2)=3
      NAO_A(i,3)=I1a
    elseif (A1=='py') then
      i=i+1
      NAO_A(i,1)=I1
      NAO_A(i,2)=4
      NAO_A(i,3)=I1a
    elseif (A1/='s '.and.A1/='pz'.and.A1/='px'.and.A1/='py') then
      exit
    endif
    if (A2=='s ') then
      i=i+1
      NAO_A(i,1)=I2
      NAO_A(i,2)=1
      NAO_A(i,3)=I2a
    elseif (A2=='pz') then
      i=i+1
      NAO_A(i,1)=I2
      NAO_A(i,2)=2
      NAO_A(i,3)=I2a
    elseif (A2=='px') then
      i=i+1
      NAO_A(i,1)=I2
      NAO_A(i,2)=3
      NAO_A(i,3)=I2a
    elseif (A2=='py') then
      i=i+1
      NAO_A(i,1)=I2
      NAO_A(i,2)=4
      NAO_A(i,3)=I2a
    elseif (A2/='s '.and.A2/='pz'.and.A2/='px'.and.A2/='py') then
      exit
    endif
    if (A3=='s ') then
      i=i+1
      NAO_A(i,1)=I3
      NAO_A(i,2)=1
      NAO_A(i,3)=I3a
    elseif (A3=='pz') then
      i=i+1
      NAO_A(i,1)=I3
      NAO_A(i,2)=2
      NAO_A(i,3)=I3a
    elseif (A3=='px') then
      i=i+1
      NAO_A(i,1)=I3
      NAO_A(i,2)=3
      NAO_A(i,3)=I3a
    elseif (A3=='py') then
      i=i+1
      NAO_A(i,1)=I3
      NAO_A(i,2)=4
      NAO_A(i,3)=I3a
    elseif (A3/='s '.and.A3/='pz'.and.A3/='px'.and.A3/='py') then
      exit
    endif
    if (A4=='s ') then
      i=i+1
      NAO_A(i,1)=I4
      NAO_A(i,2)=1
      NAO_A(i,3)=I4a
    elseif (A4=='pz') then
      i=i+1
      NAO_A(i,1)=I4
      NAO_A(i,2)=2
      NAO_A(i,3)=I4a
    elseif (A4=='px') then
      i=i+1
      NAO_A(i,1)=I4
      NAO_A(i,2)=3
      NAO_A(i,3)=I4a
    elseif (A4=='py') then
      i=i+1
      NAO_A(i,1)=I4
      NAO_A(i,2)=4
      NAO_A(i,3)=I4a
    elseif (A4/='s '.and.A4/='pz'.and.A4/='px'.and.A4/='py') then
      exit
    endif
    if (A5=='s ') then
      i=i+1
      NAO_A(i,1)=I5
      NAO_A(i,2)=1
      NAO_A(i,3)=I5a
    elseif (A5=='pz') then
      i=i+1
      NAO_A(i,1)=I5
      NAO_A(i,2)=2
      NAO_A(i,3)=I5a
    elseif (A5=='px') then
      i=i+1
      NAO_A(i,1)=I5
      NAO_A(i,2)=3
      NAO_A(i,3)=I5a
    elseif (A5=='py') then
      i=i+1
      NAO_A(i,1)=I5
      NAO_A(i,2)=4
      NAO_A(i,3)=I5a
    elseif (A5/='s '.and.A5/='pz'.and.A5/='px'.and.A5/='py') then
      exit
    endif
    if (A6=='s ') then
      i=i+1
      NAO_A(i,1)=I6
      NAO_A(i,2)=1
      NAO_A(i,3)=I6a
    elseif (A6=='pz') then
      i=i+1
      NAO_A(i,1)=I6
      NAO_A(i,2)=2
      NAO_A(i,3)=I6a
    elseif (A6=='px') then
      i=i+1
      NAO_A(i,1)=I6
      NAO_A(i,2)=3
      NAO_A(i,3)=I6a
    elseif (A6=='py') then
      i=i+1
      NAO_A(i,1)=I6
      NAO_A(i,2)=4
      NAO_A(i,3)=I6a
    elseif (A6/='s '.and.A6/='pz'.and.A6/='px'.and.A6/='py') then
      exit
    endif
    if (A7=='s ') then
      i=i+1
      NAO_A(i,1)=I7
      NAO_A(i,2)=1
      NAO_A(i,3)=I7a
    elseif (A7=='pz') then
      i=i+1
      NAO_A(i,1)=I7
      NAO_A(i,2)=2
      NAO_A(i,3)=I7a
    elseif (A7=='px') then
      i=i+1
      NAO_A(i,1)=I7
      NAO_A(i,2)=3
      NAO_A(i,3)=I7a
    elseif (A7=='py') then
      i=i+1
      NAO_A(i,1)=I7
      NAO_A(i,2)=4
      NAO_A(i,3)=I7a
    elseif (A7/='s '.and.A7/='pz'.and.A7/='px'.and.A7/='py') then
      exit
    endif
  enddo

  close(9)

  ! Get AO > NAO coefficients from FILE.33
  ! Open FILE.33 file for A
  open(unit=10,file=trim(DIR_A)//'FILE.33',status='old',IOStat=IOs)
  if (IOs /= 0) then
    write(6,'(A)')'There is no FILE.33/AO to NAO transformation file for monomer A'
    stop
  end if

  ! Positioning
  read(10,'(/,/,A)') fileline

  ! Read and store AO > NAO coefficients for A
  ! AONAOm_A(i,j)=Coeff of ith AO in jth NAO
  allocate(AONAOm_A(numAO_A,numNAO_A))
  k=floor(numAO_A/5.0)
  l=mod(numAO_A,5)
  do j=1,numNAO_A
    do i=1,k
      read(10,'(X,5F15.9)') F1,F2,F3,F4,F5
      AONAOm_A((i-1)*5+1,j)=F1
      AONAOm_A((i-1)*5+2,j)=F2
      AONAOm_A((i-1)*5+3,j)=F3
      AONAOm_A((i-1)*5+4,j)=F4
      AONAOm_A((i-1)*5+5,j)=F5
    enddo
    if (l/=0) then
      read(10,'(X,5F15.9)') F1,F2,F3,F4,F5
      AONAOm_A((k)*5+1,j)=F1
      AONAOm_A((k)*5+2,j)=F2
      AONAOm_A((k)*5+3,j)=F3
      AONAOm_A((k)*5+4,j)=F4
      AONAOm_A((k)*5+5,j)=F5
    endif
  enddo
  close(10)

  ! Open FILE.46 file for B
  open (unit=15,file=trim(DIR_B)//'FILE.46',status='old',IOStat=IOs)
  if (IOs /= 0) then
    write(6,'(A)')'There is no FILE.46/Basis info file for monomer B'
    stop
  end if

  ! Read number of AOs, NAOs, & MOs in monomer B
  do
    read(15,'(A)', IOstat=IOs) fileline
    if (fileline(3:4) == 'AO') then
      read(fileline(12:),'(i4)') numAO_B
      exit
    endif
  enddo
  rewind(15)
  do
    read(15,'(A)', IOstat=IOs) fileline
    if (fileline(3:5) == 'NAO') then
      read(fileline(12:),'(i4)') numNAO_B
      exit
    endif
  enddo

  numMO_B=numNAO_B
  rewind(15)
  allocate(AO_B(numAO_B,2),NAO_B(numNAO_B,3))
  AO_B=0
  NAO_B=0

  ! Read and store AO basis for B
  ! s=1,pz=2,px=3,py=4
  ! AO_B(i,1)=ith AO Atom
  ! AO_B(i,2)=ith AO Type
  if (mod(numAO_B,7)==0) then
    k=floor(numAO_B/7.0)
  else
    k=floor(numAO_B/7.0)+1
  endif
  i=0
  read(15,'(A)') fileline
  do j=1,k
    read(15,'(3X,I2,2X,A2,4X,I2,2X,A2,4X,I2,2X,A2,4X,I2,2X,A2,                          &
    &4X,I2,2X,A2,4X,I2,2X,A2,4X,I2,2X,A2)')I1,A1,I2,A2,I3,A3,I4,A4,I5,A5,I6,A6,I7,A7
    if (A1=='s)') then
      i=i+1
      AO_B(i,1)=I1
      AO_B(i,2)=1
    elseif (A1=='pz') then
      i=i+1
      AO_B(i,1)=I1
      AO_B(i,2)=2
    elseif (A1=='px') then
      i=i+1
      AO_B(i,1)=I1
      AO_B(i,2)=3
    elseif (A1=='py') then
      i=i+1
      AO_B(i,1)=I1
      AO_B(i,2)=4
    elseif (A1/='s '.and.A1/='pz'.and.A1/='px'.and.A1/='py') then
      exit
    endif
    if (A2=='s)') then
      i=i+1
      AO_B(i,1)=I2
      AO_B(i,2)=1
    elseif (A2=='pz') then
      i=i+1
      AO_B(i,1)=I2
      AO_B(i,2)=2
    elseif (A2=='px') then
      i=i+1
      AO_B(i,1)=I2
      AO_B(i,2)=3
    elseif (A2=='py') then
      i=i+1
      AO_B(i,1)=I2
      AO_B(i,2)=4
    elseif (A2/='s '.and.A2/='pz'.and.A2/='px'.and.A2/='py') then
      exit
    endif
    if (A3=='s)') then
      i=i+1
      AO_B(i,1)=I3
      AO_B(i,2)=1
    elseif (A3=='pz') then
      i=i+1
      AO_B(i,1)=I3
      AO_B(i,2)=2
    elseif (A3=='px') then
      i=i+1
      AO_B(i,1)=I3
      AO_B(i,2)=3
    elseif (A3=='py') then
      i=i+1
      AO_B(i,1)=I3
      AO_B(i,2)=4
    elseif (A3/='s '.and.A3/='pz'.and.A3/='px'.and.A3/='py') then
      exit
    endif
    if (A4=='s)') then
      i=i+1
      AO_B(i,1)=I4
      AO_B(i,2)=1
    elseif (A4=='pz') then
      i=i+1
      AO_B(i,1)=I4
      AO_B(i,2)=2
    elseif (A4=='px') then
      i=i+1
      AO_B(i,1)=I4
      AO_B(i,2)=3
    elseif (A4=='py') then
      i=i+1
      AO_B(i,1)=I4
      AO_B(i,2)=4
    elseif (A4/='s '.and.A4/='pz'.and.A4/='px'.and.A4/='py') then
      exit
    endif
    if (A5=='s)') then
      i=i+1
      AO_B(i,1)=I5
      AO_B(i,2)=1
    elseif (A5=='pz') then
      i=i+1
      AO_B(i,1)=I5
      AO_B(i,2)=2
    elseif (A5=='px') then
      i=i+1
      AO_B(i,1)=I5
      AO_B(i,2)=3
    elseif (A5=='py') then
      i=i+1
      AO_B(i,1)=I5
      AO_B(i,2)=4
    elseif (A5/='s '.and.A5/='pz'.and.A5/='px'.and.A5/='py') then
      exit
    endif
    if (A6=='s)') then
      i=i+1
      AO_B(i,1)=I6
      AO_B(i,2)=1
    elseif (A6=='pz') then
      i=i+1
      AO_B(i,1)=I6
      AO_B(i,2)=2
    elseif (A6=='px') then
      i=i+1
      AO_B(i,1)=I6
      AO_B(i,2)=3
    elseif (A6=='py') then
      i=i+1
      AO_B(i,1)=I6
      AO_B(i,2)=4
    elseif (A6/='s '.and.A6/='pz'.and.A6/='px'.and.A6/='py') then
      exit
    endif
    if (A7=='s)') then
      i=i+1
      AO_B(i,1)=I7
      AO_B(i,2)=1
    elseif (A7=='pz') then
      i=i+1
      AO_B(i,1)=I7
      AO_B(i,2)=2
    elseif (A7=='px') then
      i=i+1
      AO_B(i,1)=I7
      AO_B(i,2)=3
    elseif (A7=='py') then
      i=i+1
      AO_B(i,1)=I7
      AO_B(i,2)=4
    elseif (A7/='s '.and.A7/='pz'.and.A7/='px'.and.A7/='py') then
      exit
    endif
  enddo

  ! Read and store NAO basis for B
  ! s=1,pz=2,px=3,py=4
  ! NAO_B(i,1)=ith NAO Atom
  ! NAO_B(i,2)=ith NAO Type
  ! NAO_B(i,3)=ith NAO Shell
  if (mod(numNAO_B,7)==0) then
    k=floor(numNAO_B/7.0)
  else
    k=floor(numNAO_B/7.0)+1
  endif
  i=0
  ! Positioning
  ! Planarity, if (NAO_A(i,3)==2.and.NAO_A(i,2)==2)
  read(15,'(A)') fileline
  do j=1,k
    read(15,'(3X,I2,2X,I1,A2,3X,I2,2X,I1,A2,3X,I2,2X,I1,A2,3X,I2,2X,I1,A2,              &
         &3X,I2,2X,I1,A2,3X,I2,2X,I1,A2,3X,I2,2X,I1,A2)')I1,I1a,A1,I2,I2a,A2,I3,        &
         &I3a,A3,I4,I4a,A4,I5,I5a,A5,I6,I6a,A6,I7,I7a,A7
    if (A1=='s ') then
      i=i+1
      NAO_B(i,1)=I1
      NAO_B(i,2)=1
      NAO_B(i,3)=I1a
    elseif (A1=='pz') then
      i=i+1
      NAO_B(i,1)=I1
      NAO_B(i,2)=2
      NAO_B(i,3)=I1a
    elseif (A1=='px') then
      i=i+1
      NAO_B(i,1)=I1
      NAO_B(i,2)=3
      NAO_B(i,3)=I1a
    elseif (A1=='py') then
      i=i+1
      NAO_B(i,1)=I1
      NAO_B(i,2)=4
      NAO_B(i,3)=I1a
    elseif (A1/='s '.and.A1/='pz'.and.A1/='px'.and.A1/='py') then
      exit
    endif
    if (A2=='s ') then
      i=i+1
      NAO_B(i,1)=I2
      NAO_B(i,2)=1
      NAO_B(i,3)=I2a
    elseif (A2=='pz') then
      i=i+1
      NAO_B(i,1)=I2
      NAO_B(i,2)=2
      NAO_B(i,3)=I2a
    elseif (A2=='px') then
      i=i+1
      NAO_B(i,1)=I2
      NAO_B(i,2)=3
      NAO_B(i,3)=I2a
    elseif (A2=='py') then
      i=i+1
      NAO_B(i,1)=I2
      NAO_B(i,2)=4
      NAO_B(i,3)=I2a
    elseif (A2/='s '.and.A2/='pz'.and.A2/='px'.and.A2/='py') then
      exit
    endif
    if (A3=='s ') then
      i=i+1
      NAO_B(i,1)=I3
      NAO_B(i,2)=1
      NAO_B(i,3)=I3a
    elseif (A3=='pz') then
      i=i+1
      NAO_B(i,1)=I3
      NAO_B(i,2)=2
      NAO_B(i,3)=I3a
    elseif (A3=='px') then
      i=i+1
      NAO_B(i,1)=I3
      NAO_B(i,2)=3
      NAO_B(i,3)=I3a
    elseif (A3=='py') then
      i=i+1
      NAO_B(i,1)=I3
      NAO_B(i,2)=4
      NAO_B(i,3)=I3a
    elseif (A3/='s '.and.A3/='pz'.and.A3/='px'.and.A3/='py') then
      exit
    endif
    if (A4=='s ') then
      i=i+1
      NAO_B(i,1)=I4
      NAO_B(i,2)=1
      NAO_B(i,3)=I4a
    elseif (A4=='pz') then
      i=i+1
      NAO_B(i,1)=I4
      NAO_B(i,2)=2
      NAO_B(i,3)=I4a
    elseif (A4=='px') then
      i=i+1
      NAO_B(i,1)=I4
      NAO_B(i,2)=3
      NAO_B(i,3)=I4a
    elseif (A4=='py') then
      i=i+1
      NAO_B(i,1)=I4
      NAO_B(i,2)=4
      NAO_B(i,3)=I4a
    elseif (A4/='s '.and.A4/='pz'.and.A4/='px'.and.A4/='py') then
      exit
    endif
    if (A5=='s ') then
      i=i+1
      NAO_B(i,1)=I5
      NAO_B(i,2)=1
      NAO_B(i,3)=I5a
    elseif (A5=='pz') then
      i=i+1
      NAO_B(i,1)=I5
      NAO_B(i,2)=2
      NAO_B(i,3)=I5a
    elseif (A5=='px') then
      i=i+1
      NAO_B(i,1)=I5
      NAO_B(i,2)=3
      NAO_B(i,3)=I5a
    elseif (A5=='py') then
      i=i+1
      NAO_B(i,1)=I5
      NAO_B(i,2)=4
      NAO_B(i,3)=I5a
    elseif (A5/='s '.and.A5/='pz'.and.A5/='px'.and.A5/='py') then
      exit
    endif
    if (A6=='s ') then
      i=i+1
      NAO_B(i,1)=I6
      NAO_B(i,2)=1
      NAO_B(i,3)=I6a
    elseif (A6=='pz') then
      i=i+1
      NAO_B(i,1)=I6
      NAO_B(i,2)=2
      NAO_B(i,3)=I6a
    elseif (A6=='px') then
      i=i+1
      NAO_B(i,1)=I6
      NAO_B(i,2)=3
      NAO_B(i,3)=I6a
    elseif (A6=='py') then
      i=i+1
      NAO_B(i,1)=I6
      NAO_B(i,2)=4
      NAO_B(i,3)=I6a
    elseif (A6/='s '.and.A6/='pz'.and.A6/='px'.and.A6/='py') then
      exit
    endif
    if (A7=='s ') then
      i=i+1
      NAO_B(i,1)=I7
      NAO_B(i,2)=1
      NAO_B(i,3)=I7a
    elseif (A7=='pz') then
      i=i+1
      NAO_B(i,1)=I7
      NAO_B(i,2)=2
      NAO_B(i,3)=I7a
    elseif (A7=='px') then
      i=i+1
      NAO_B(i,1)=I7
      NAO_B(i,2)=3
      NAO_B(i,3)=I7a
    elseif (A7=='py') then
      i=i+1
      NAO_B(i,1)=I7
      NAO_B(i,2)=4
      NAO_B(i,3)=I7a
    elseif (A7/='s '.and.A7/='pz'.and.A7/='px'.and.A7/='py') then
      exit
    endif
  enddo

  close(15)

  ! Get AO > NAO coefficients from FILE.33
  ! Open FILE.33 file for B
  open(unit=16,file=trim(DIR_B)//'FILE.33',status='old',IOStat=IOs)
  if (IOs /= 0) then
    write(6,'(A)')'There is no FILE.33/AO to NAO transformation file for monomer B'
    stop
  end if

  ! Positioning
  read(16,'(/,/,A)') fileline

  ! Read and store AO > NAO coefficients for B
  ! AONAOm_B(i,j)=Coeff of ith AO in jth NAO
  allocate(AONAOm_B(numAO_B,numNAO_B))
  k=floor(numAO_B/5.0)
  l=mod(numAO_B,5)
  do j=1,numNAO_B
    do i=1,k
      read(16,'(X,5F15.9)') F1,F2,F3,F4,F5
      AONAOm_B((i-1)*5+1,j)=F1
      AONAOm_B((i-1)*5+2,j)=F2
      AONAOm_B((i-1)*5+3,j)=F3
      AONAOm_B((i-1)*5+4,j)=F4
      AONAOm_B((i-1)*5+5,j)=F5
    enddo
    if (l/=0) then
      read(16,'(X,5F15.9)') F1,F2,F3,F4,F5
      AONAOm_B((k)*5+1,j)=F1
      AONAOm_B((k)*5+2,j)=F2
      AONAOm_B((k)*5+3,j)=F3
      AONAOm_B((k)*5+4,j)=F4
      AONAOm_B((k)*5+5,j)=F5
    endif
  enddo
  close(16)

  ! Get NAO > MO coefficients from FILE.49
  ! Open FILE.49 file for A
  open(unit=11,file=trim(DIR_A)//'FILE.49',status='old',IOStat=IOs)
  if (IOs /= 0) then
    write(6,'(A)')'There is no FILE.49/NAO to MO transformation file for monomer A'
    stop
  end if

  ! Positioning
  read(11,'(/,/,A)') fileline

  ! Read and store NAO > MO coefficients for A
  ! NAOMOm_A(i,j)=Coeff of ith NAO in jth MO
  allocate(NAOMOm_A(numNAO_A,numMO_A))
  k=floor(numAO_A/5.0)
  l=mod(numAO_A,5)
  do j=1,numMO_A
    do i=1,k
      read(11,'(X,5F15.9)') F1,F2,F3,F4,F5
      NAOMOm_A((i-1)*5+1,j)=F1
      NAOMOm_A((i-1)*5+2,j)=F2
      NAOMOm_A((i-1)*5+3,j)=F3
      NAOMOm_A((i-1)*5+4,j)=F4
      NAOMOm_A((i-1)*5+5,j)=F5
    enddo
      if (l/=0) then
        read(11,'(X,5F15.9)') F1,F2,F3,F4,F5
        NAOMOm_A((k)*5+1,j)=F1
        NAOMOm_A((k)*5+2,j)=F2
        NAOMOm_A((k)*5+3,j)=F3
        NAOMOm_A((k)*5+4,j)=F4
        NAOMOm_A((k)*5+5,j)=F5
      endif
  enddo
  close(11)

  ! Get NAO > MO coefficients from FILE.49
  ! Open FILE.49 file for B
  open(unit=17,file=trim(DIR_B)//'FILE.49',status='old',IOStat=IOs)
  if (IOs /= 0) then
    write(6,'(A)')'There is no FILE.49/NAO to MO transformation file for monomer B'
    stop
  end if

  ! Positioning
  read(17,'(/,/,A)') fileline

  ! Read and store NAO > MO coefficients for B
  ! NAOMOm_B(i,j)=Coeff of ith NAO in jth MO
  allocate(NAOMOm_B(numNAO_B,numMO_B))
  k=floor(numAO_B/5.0)
  l=mod(numAO_B,5)
  do j=1,numMO_B
    do i=1,k
      read(17,'(X,5F15.9)') F1,F2,F3,F4,F5
      NAOMOm_B((i-1)*5+1,j)=F1
      NAOMOm_B((i-1)*5+2,j)=F2
      NAOMOm_B((i-1)*5+3,j)=F3
      NAOMOm_B((i-1)*5+4,j)=F4
      NAOMOm_B((i-1)*5+5,j)=F5
    enddo
      if (l/=0) then
        read(17,'(X,5F15.9)') F1,F2,F3,F4,F5
        NAOMOm_B((k)*5+1,j)=F1
        NAOMOm_B((k)*5+2,j)=F2
        NAOMOm_B((k)*5+3,j)=F3
        NAOMOm_B((k)*5+4,j)=F4
        NAOMOm_B((k)*5+5,j)=F5
      endif
  enddo
  close(17)

  ! All P orbitals if planarity=no 
  ! Just Pz orbitals if planarity=yes
  ! Get Valence P NAOs of A
  allocate(buffer(numNAO_A,4))
  k=0
  ! NAO xyz stored in Valence
!    16 A 3
!    P 3  0.3925  0.0012 -0.0089
!    20.9642000  0.0402487
!    4.8033100  0.2375940
!    16 A 3
!    P 3  0.0002  0.3880 -0.0035
!    20.9642000  0.0402487
!       ...
!    16 A 3
!       ...
  if (TRIM(planarity)=='no') then
    do i=1,numNAO_A
      if (NAO_A(i,3)==2.and.NAO_A(i,2)/=1) then
        k=k+1
        buffer(k,:)=NAO_A(i,:)
        buffer(k,4)=i
      endif
    enddo
    numVAL_A=k
    allocate(Valence_A(numVAL_A,4))
    Valence_A=buffer
  elseif (TRIM(planarity)=='yes') then
    do i=1,numNAO_A
      if (NAO_A(i,3)==2.and.NAO_A(i,2)==2) then
        k=k+1
        buffer(k,:)=NAO_A(i,:)
        buffer(k,4)=i
      endif
    enddo
    numVAL_A=k
    allocate(Valence_A(numVAL_A,4))
    Valence_A=buffer
!   do i=1,numVAL_A
!     do j=1,2
!       do k=1,4
!         AONAOm_A(PH_A(Valence_A(i,4),j,k),Valence_A(i,4))=0.0d0
!       enddo
!     enddo
!   enddo
  endif
  deallocate(buffer)
  ! Get Valence P NAOs of B
  allocate(buffer(numNAO_B,4))
  k=0
  if (TRIM(planarity)=='no') then
    do i=1,numNAO_B
      if (NAO_B(i,3)==2.and.NAO_B(i,2)/=1) then
        k=k+1
        buffer(k,:)=NAO_B(i,:)
        buffer(k,4)=i
      endif
    enddo
    numVAL_B=k
    allocate(Valence_B(numVAL_B,4))
    Valence_B=buffer
  elseif (TRIM(planarity)=='yes') then
    do i=1,numNAO_B
      if (NAO_B(i,3)==2.and.NAO_B(i,2)==2) then
        k=k+1
        buffer(k,:)=NAO_B(i,:)
        buffer(k,4)=i
      endif
    enddo
    numVAL_B=k
    allocate(Valence_B(numVAL_B,4))
    Valence_B=buffer
!   do i=1,numVAL_B
!     do j=1,2
!       do k=1,4
!         AONAOm_B(PH_B(Valence_B(i,4),j,k),Valence_B(i,4))=0.0d0
!       enddo
!     enddo
!   enddo
  endif
  deallocate(buffer)
  ! Print output
  ! A
  write(6,'(/,A,I3,A,/)'),'There are ',numVAL_A,' valence P type NAOs for A:'
  write(6,'(4A8,/)'),'NAO #','Atom #:','Type:','Shell:'
  do i=1,numVAL_A
    if (Valence_A(i,2)==2) then
      write(6,'(2I8,A8,I8)'),Valence_A(i,4),Valence_A(i,1),'pz',Valence_A(i,3)
    elseif (Valence_A(i,2)==3) then
      write(6,'(2I8,A8,I8)'),Valence_A(i,4),Valence_A(i,1),'px',Valence_A(i,3)
    elseif (Valence_A(i,2)==4) then
      write(6,'(2I8,A8,I8)'),Valence_A(i,4),Valence_A(i,1),'py',Valence_A(i,3)
    endif 
  enddo
  ! B
  write(6,'(/,A,I3,A,/)'),'There are ',numVAL_B,' valence P type NAOs for B:'
  write(6,'(4A8,/)'),'NAO #','Atom #:','Type:','Shell:'
  do i=1,numVAL_B
    if (Valence_B(i,2)==2) then
      write(6,'(2I8,A8,I8)'),Valence_B(i,4),Valence_B(i,1),'pz',Valence_B(i,3)
    elseif (Valence_B(i,2)==3) then
      write(6,'(2I8,A8,I8)'),Valence_B(i,4),Valence_B(i,1),'px',Valence_B(i,3)
    elseif (Valence_B(i,2)==4) then
      write(6,'(2I8,A8,I8)'),Valence_B(i,4),Valence_B(i,1),'py',Valence_B(i,3)
    endif
  enddo

  ! Stupid PH tensor to hold list of px, py, and pz AOs in NAOs for A
  allocate(PH_A(numNAO_A,3,4))
  PH_A=0.0d0
  write(6,'(/,A,I3,A,/)'),'The non-zero coefficients of the AOs in the valence P type NAOs of A:'
  write(6,'(A15,4A8,/)'),'Coefficient:','AO #','Atom #:','Type:'
  do i=1,numVAL_A
    l=1
    m=1
    n=1
    write(6,*),'NAO',Valence_A(i,4)
    do j=1,numAO_A
      if(AO_A(j,2)==3.and.abs(AONAOm_A(j,Valence_A(i,4)))>=0.0d0.and.AO_A(j,1)==Valence_A(i,1)) then
        write(6,'(F15.9,2I8,A8,I8)'),AONAOm_A(j,Valence_A(i,4)),j,AO_A(j,1),'px'
        PH_A(Valence_A(i,4),1,l)=j
        l=l+1
      elseif(AO_A(j,2)==4.and.abs(AONAOm_A(j,Valence_A(i,4)))>=0.0d0.and.AO_A(j,1)==Valence_A(i,1)) then
        write(6,'(F15.9,2I8,A8,I8)'),AONAOm_A(j,Valence_A(i,4)),j,AO_A(j,1),'py'
        PH_A(Valence_A(i,4),2,m)=j
        m=m+1
      elseif(AO_A(j,2)==2.and.abs(AONAOm_A(j,Valence_A(i,4)))>=0.0d0.and.AO_A(j,1)==Valence_A(i,1)) then
        write(6,'(F15.9,2I8,A8,I8)'),AONAOm_A(j,Valence_A(i,4)),j,AO_A(j,1),'pz'
        PH_A(Valence_A(i,4),3,n)=j
        n=n+1
      endif
    enddo
  enddo
  ! Stupid PH tensor to hold list of px, py, and pz AOs in NAOs for B
  allocate(PH_B(numNAO_B,3,4))
  PH_B=0.0d0
  write(6,'(/,A,I3,A,/)'),'The non-zero coefficients of the AOs in the valence P type NAOs of B:'
  write(6,'(A15,4A8,/)'),'Coefficient:','AO #','Atom #:','Type:'
  do i=1,numVAL_B
    l=1
    m=1
    n=1
    write(6,*),'NAO',Valence_B(i,4)
    do j=1,numAO_B
      if(AO_B(j,2)==3.and.abs(AONAOm_B(j,Valence_B(i,4)))>=0.0d0.and.AO_B(j,1)==Valence_B(i,1)) then
        write(6,'(F15.9,2I8,A8,I8)'),AONAOm_B(j,Valence_B(i,4)),j,AO_B(j,1),'px'
        PH_B(Valence_B(i,4),1,l)=j
        l=l+1
      elseif(AO_B(j,2)==4.and.abs(AONAOm_B(j,Valence_B(i,4)))>=0.0d0.and.AO_B(j,1)==Valence_B(i,1)) then
        write(6,'(F15.9,2I8,A8,I8)'),AONAOm_B(j,Valence_B(i,4)),j,AO_B(j,1),'py'
        PH_B(Valence_B(i,4),2,m)=j
        m=m+1
      elseif(AO_B(j,2)==2.and.abs(AONAOm_B(j,Valence_B(i,4)))>=0.0d0.and.AO_B(j,1)==Valence_B(i,1)) then
        write(6,'(F15.9,2I8,A8,I8)'),AONAOm_B(j,Valence_B(i,4)),j,AO_B(j,1),'pz'
        PH_B(Valence_B(i,4),3,n)=j
        n=n+1
      endif
    enddo
  enddo
  !  If Planar then only Pz contributions to NAOs
  if (TRIM(planarity)=='yes') then
    do i=1,numVAL_A
      do j=1,2
        do k=1,4
          AONAOm_A(PH_A(Valence_A(i,4),j,k),Valence_A(i,4))=0.0d0
        enddo
      enddo
    enddo
    do i=1,numVAL_B
      do j=1,2
        do k=1,4
          AONAOm_B(PH_B(Valence_B(i,4),j,k),Valence_B(i,4))=0.0d0
        enddo
      enddo
    enddo
  endif

  ! HOMO A
  write(6,'(/,A,I3,/)'),'HOMO for A is MO #',nHOMO_A
  write(6,'(A,I3,A,/)'),'The non-zero coefficients of the',numVAL_A,' NAOs in the HOMO of A:'
  write(6,'(A15,4A8,/)'),'Coefficient:','NAO #','Atom #:','Type:','Shell:'
  do j=1,numVAL_A
    if (Valence_A(j,2)==2.and.NAOMOm_A(Valence_A(j,4),nHOMO_A)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_A(Valence_A(j,4),nHOMO_A),Valence_A(j,4),Valence_A(j,1),'pz',Valence_A(j,3)
    elseif (Valence_A(j,2)==3.and.NAOMOm_A(Valence_A(j,4),nHOMO_A)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_A(Valence_A(j,4),nHOMO_A),Valence_A(j,4),Valence_A(j,1),'px',Valence_A(j,3)
    elseif (Valence_A(j,2)==4.and.NAOMOm_A(Valence_A(j,4),nHOMO_A)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_A(Valence_A(j,4),nHOMO_A),Valence_A(j,4),Valence_A(j,1),'py',Valence_A(j,3)
    endif
  enddo
  ! LUMO A
  write(6,'(/,A,I3,/)'),'LUMO for A is MO #',nLUMO_A
  write(6,'(A,I3,A,/)'),'The non-zero coefficients of the',numVAL_A,' NAOs in the LUMO of A:'
  write(6,'(A15,4A8,/)'),'Coefficient:','NAO #','Atom #:','Type:','Shell:'
  do j=1,numVAL_A
    if (Valence_A(j,2)==2.and.NAOMOm_A(Valence_A(j,4),nLUMO_A)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_A(Valence_A(j,4),nLUMO_A),Valence_A(j,4),Valence_A(j,1),'pz',Valence_A(j,3)
    elseif (Valence_A(j,2)==3.and.NAOMOm_A(Valence_A(j,4),nLUMO_A)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_A(Valence_A(j,4),nLUMO_A),Valence_A(j,4),Valence_A(j,1),'px',Valence_A(j,3)
    elseif (Valence_A(j,2)==4.and.NAOMOm_A(Valence_A(j,4),nLUMO_A)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_A(Valence_A(j,4),nLUMO_A),Valence_A(j,4),Valence_A(j,1),'py',Valence_A(j,3)
    endif
  enddo
  ! HOMO B
  write(6,'(/,A,I3,/)'),'HOMO for B is MO #',nHOMO_B
  write(6,'(A,I3,A,/)'),'The non-zero coefficients of the',numVAL_B,' NAOs in the HOMO of B:'
  write(6,'(A15,4A8,/)'),'Coefficient:','NAO #','Atom #:','Type:','Shell:'
  do j=1,numVAL_B
    if (Valence_B(j,2)==2.and.NAOMOm_B(Valence_B(j,4),nHOMO_B)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_B(Valence_B(j,4),nHOMO_B),Valence_B(j,4),Valence_B(j,1),'pz',Valence_B(j,3)
    elseif (Valence_B(j,2)==3.and.NAOMOm_B(Valence_B(j,4),nHOMO_B)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_B(Valence_B(j,4),nHOMO_B),Valence_B(j,4),Valence_B(j,1),'px',Valence_B(j,3)
    elseif (Valence_B(j,2)==4.and.NAOMOm_B(Valence_B(j,4),nHOMO_B)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_B(Valence_B(j,4),nHOMO_B),Valence_B(j,4),Valence_B(j,1),'py',Valence_B(j,3)
    endif
  enddo
  ! LUMO B
  write(6,'(/,A,I3,/)'),'LUMO for B is MO #',nLUMO_B
  write(6,'(A,I3,A,/)'),'The non-zero coefficients of the',numVAL_B,' NAOs in the LUMO of B:'
  write(6,'(A15,4A8,/)'),'Coefficient:','NAO #','Atom #:','Type:','Shell:'
  do j=1,numVAL_B
    if (Valence_B(j,2)==2.and.NAOMOm_B(Valence_B(j,4),nLUMO_B)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_B(Valence_B(j,4),nLUMO_B),Valence_B(j,4),Valence_B(j,1),'pz',Valence_B(j,3)
    elseif (Valence_B(j,2)==3.and.NAOMOm_B(Valence_B(j,4),nLUMO_B)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_B(Valence_B(j,4),nLUMO_B),Valence_B(j,4),Valence_B(j,1),'px',Valence_B(j,3)
    elseif (Valence_B(j,2)==4.and.NAOMOm_B(Valence_B(j,4),nLUMO_B)/=0.0d0) then
      write(6,'(F15.9,2I8,A8,I8)'),NAOMOm_B(Valence_B(j,4),nLUMO_B),Valence_B(j,4),Valence_B(j,1),'py',Valence_B(j,3)
    endif
  enddo

  ! Create Simput.out file for Simple input
  open(unit=12,file='Simple.inp',status='new')

  ! Write Simput.out header
  write(12,'(A7,X,A)') 'Title :',title
  write(12,'(A12,X,i4)')'Processors :',numProc
  write(12,'(A8)') 'Method :'
  write(12,'(A)') method
  if (noLJwarn) write(12,'(A)') 'noLJwarn'
  ! Write GeomA
  write(12,'(A7)') 'GeomA :'
  do i=1,numA_A
    if (COORDS_A(i,1) == 5) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'B',COORDS_A(i,2),COORDS_A(i,3),COORDS_A(i,4)
    elseif (COORDS_A(i,1) == 6) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'C',COORDS_A(i,2),COORDS_A(i,3),COORDS_A(i,4)
    elseif (COORDS_A(i,1) == 7) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'N',COORDS_A(i,2),COORDS_A(i,3),COORDS_A(i,4)
    elseif (COORDS_A(i,1) == 8) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'O',COORDS_A(i,2),COORDS_A(i,3),COORDS_A(i,4)
    elseif (COORDS_A(i,1) == 1) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'H',COORDS_A(i,2),COORDS_A(i,3),COORDS_A(i,4)
    elseif (COORDS_A(i,1) == 9) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'F',COORDS_A(i,2),COORDS_A(i,3),COORDS_A(i,4)
    endif
  enddo
  ! Write GeomB
  write(12,'(A7)') 'GeomB :'
  do i=1,numA_B
    if (COORDS_B(i,1) == 5) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'C',COORDS_B(i,2),COORDS_B(i,3),COORDS_B(i,4)
    elseif (COORDS_B(i,1) == 6) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'C',COORDS_B(i,2),COORDS_B(i,3),COORDS_B(i,4)
    elseif (COORDS_B(i,1) == 7) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'N',COORDS_B(i,2),COORDS_B(i,3),COORDS_B(i,4)
    elseif (COORDS_B(i,1) == 8) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'O',COORDS_B(i,2),COORDS_B(i,3),COORDS_B(i,4)
    elseif (COORDS_B(i,1) == 1) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'H',COORDS_B(i,2),COORDS_B(i,3),COORDS_B(i,4)
    elseif (COORDS_B(i,1) == 9) then
      write(12,'(1X,A1,6X,F12.9,3X,F12.9,3X,F12.9)')'F',COORDS_B(i,2),COORDS_B(i,3),COORDS_B(i,4)
    endif
  enddo
  ! Write Constants
  write(12,'(A7)') 'Const :'
  write(12,'(A6,F12.7)') 'T X   ',tx
  write(12,'(A6,F12.7)') 'T Y   ',ty
  write(12,'(A6,F12.7)') 'T Z   ',tz
  write(12,'(A6,F12.7)') 'R X   ',rx
  write(12,'(A6,F12.7)') 'R Y   ',ry
  write(12,'(A6,F12.7)') 'R Z   ',rz
  ! Scan and Grid left blank for manual input
  write(12,'(A6)') 'Scan :'
  ! Write NAOs
  write(12,'(A5)') 'NAO :'
! ! Thresholds for AO>NAO coefficients
! do i=1,numVAL_A
!   do j=1,3
!     do k=1,3
!       if (ABS(AONAOm_A(PH_A(Valence_A(i,4),j,k),Valence_A(i,4)))<0.1d0) then
!         AONAOm_A(PH_A(Valence_A(i,4),j,k),Valence_A(i,4))=0.0d0
!       endif
!     enddo
!     if (ABS(AONAOm_A(PH_A(Valence_A(i,4),j,4),Valence_A(i,4)))<0.01d0) then
!       AONAOm_A(PH_A(Valence_A(i,4),j,4),Valence_A(i,4))=0.0d0
!     endif
!   enddo
! enddo
! do i=1,numVAL_B
!   do j=1,3
!     do k=1,3
!       if (ABS(AONAOm_B(PH_B(Valence_B(i,4),j,k),Valence_B(i,4)))<0.1d0) then
!         AONAOm_B(PH_B(Valence_B(i,4),j,k),Valence_B(i,4))=0.0d0
!       endif
!     enddo
!     if (ABS(AONAOm_B(PH_B(Valence_B(i,4),j,4),Valence_B(i,4)))<0.01d0) then
!       AONAOm_B(PH_B(Valence_B(i,4),j,4),Valence_B(i,4))=0.0d0
!     endif
!   enddo
! enddo


  ! NAOs for A
  do i=1,numVAL_A
! Cutoff minipatch
    if ((ABS(NAOMOm_A(Valence_A(i,4),nHOMO_A))>=cutoff_val).or.(ABS(NAOMOm_A(Valence_A(i,4),nLUMO_A))>=cutoff_val)) then
      write(12,'(I4,X,A1,X,I1)')Valence_A(i,1),'A',3
      write(12,'(A3,X,F7.4,X,F7.4,X,F7.4)')'P 3',AONAOm_A(PH_A(Valence_A(i,4),1,1),Valence_A(i,4)),AONAOm_A(PH_A    &
          &(Valence_A(i,4),2,1),Valence_A(i,4)),AONAOm_A(PH_A(Valence_A(i,4),3,1),Valence_A(i,4))
      write(12,'(2F11.7)')BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,2),BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,3)
      write(12,'(2F11.7)')BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,4),BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,5)
      write(12,'(2F11.7)')BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,6),BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,7)
      write(12,'(A3,X,F7.4,X,F7.4,X,F7.4)')'P 1',AONAOm_A(PH_A(Valence_A(i,4),1,2),Valence_A(i,4)),AONAOm_A(PH_A    &
          &(Valence_A(i,4),2,2),Valence_A(i,4)),AONAOm_A(PH_A(Valence_A(i,4),3,2),Valence_A(i,4))
      write(12,'(2F11.7)')BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,8),BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,9)
      write(12,'(A3,X,F7.4,X,F7.4,X,F7.4)')'P 1',AONAOm_A(PH_A(Valence_A(i,4),1,3),Valence_A(i,4)),AONAOm_A(PH_A    &
          &(Valence_A(i,4),2,3),Valence_A(i,4)),AONAOm_A(PH_A(Valence_A(i,4),3,3),Valence_A(i,4))
      write(12,'(2F11.7)')BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,10),BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,11)
   ! Diffuse funct.:
   ! write(12,'(A3,X,F7.4,X,F7.4,X,F7.4)')'P 1',AONAOm_A(PH_A(Valence_A(i,4),1,4),Valence_A(i,4)),AONAOm_A(PH_A    &
   !       &(Valence_A(i,4),2,4),Valence_A(i,4)),AONAOm_A(PH_A(Valence_A(i,4),3,4),Valence_A(i,4))
   ! write(12,'(2F11.7)')BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,12),BASISm(floor(COORDS_A(Valence_A(i,1),1))-5,13)
    endif
  enddo
  ! NAOs for B
  do i=1,numVAL_B
! Cutoff minipatch
    if ((ABS(NAOMOm_B(Valence_B(i,4),nHOMO_B))>=cutoff_val).or.(ABS(NAOMOm_B(Valence_B(i,4),nLUMO_B))>=cutoff_val)) then
      write(12,'(I4,X,A1,X,I1)')Valence_B(i,1),'B',3
      write(12,'(A3,X,F7.4,X,F7.4,X,F7.4)')'P 3',AONAOm_B(PH_B(Valence_B(i,4),1,1),Valence_B(i,4)),AONAOm_B(PH_B    &
          &(Valence_B(i,4),2,1),Valence_B(i,4)),AONAOm_B(PH_B(Valence_B(i,4),3,1),Valence_B(i,4))
      write(12,'(2F11.7)')BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,2),BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,3)
      write(12,'(2F11.7)')BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,4),BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,5)
      write(12,'(2F11.7)')BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,6),BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,7)
      write(12,'(A3,X,F7.4,X,F7.4,X,F7.4)')'P 1',AONAOm_B(PH_B(Valence_B(i,4),1,2),Valence_B(i,4)),AONAOm_B(PH_B    &
          &(Valence_B(i,4),2,2),Valence_B(i,4)),AONAOm_B(PH_B(Valence_B(i,4),3,2),Valence_B(i,4))
      write(12,'(2F11.7)')BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,8),BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,9)
      write(12,'(A3,X,F7.4,X,F7.4,X,F7.4)')'P 1',AONAOm_B(PH_B(Valence_B(i,4),1,3),Valence_B(i,4)),AONAOm_B(PH_B    &
          &(Valence_B(i,4),2,3),Valence_B(i,4)),AONAOm_B(PH_B(Valence_B(i,4),3,3),Valence_B(i,4))
      write(12,'(2F11.7)')BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,10),BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,11)
   ! Diffuse function
   ! write(12,'(A3,X,F7.4,X,F7.4,X,F7.4)')'P 1',AONAOm_B(PH_B(Valence_B(i,4),1,4),Valence_B(i,4)),AONAOm_B(PH_B    &
   !       &(Valence_B(i,4),2,4),Valence_B(i,4)),AONAOm_B(PH_B(Valence_B(i,4),3,4),Valence_B(i,4))
   ! write(12,'(2F11.7)')BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,12),BASISm(floor(COORDS_B(Valence_B(i,1),1))-5,13)
    endif
  enddo
  ! Write MOs for A
  write(12,'(A5)') 'MOs :'
  ! HOMO A
  write(12,'(A2)') 'hA'
  j_trunc=1
  do j=1,numVAL_A
    if (Valence_A(j,2)==2.and.(ABS(NAOMOm_A(Valence_A(j,4),nHOMO_A))>=cutoff_val).or.(ABS(NAOMOm_A(Valence_A(j,4),nLUMO_A))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_A(Valence_A(j,4),nHOMO_A)
    elseif (Valence_A(j,2)==3.and.(ABS(NAOMOm_A(Valence_A(j,4),nHOMO_A))>=cutoff_val).or.(ABS(NAOMOm_A(Valence_A(j,4),nLUMO_A))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_A(Valence_A(j,4),nHOMO_A)
    elseif (Valence_A(j,2)==4.and.(ABS(NAOMOm_A(Valence_A(j,4),nHOMO_A))>=cutoff_val).or.(ABS(NAOMOm_A(Valence_A(j,4),nLUMO_A))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_A(Valence_A(j,4),nHOMO_A)
! Cutoff minipatch
    else
      j_trunc=j_trunc-1
    endif
    j_trunc=j_trunc+1
  enddo
  ! LUMO A
  j_trunc=1
  write(12,'(A2)') 'lA'
  do j=1,numVAL_A
    if (Valence_A(j,2)==2.and.(ABS(NAOMOm_A(Valence_A(j,4),nHOMO_A))>=cutoff_val).or.(ABS(NAOMOm_A(Valence_A(j,4),nLUMO_A))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_A(Valence_A(j,4),nLUMO_A)
    elseif (Valence_A(j,2)==3.and.(ABS(NAOMOm_A(Valence_A(j,4),nHOMO_A))>=cutoff_val).or.(ABS(NAOMOm_A(Valence_A(j,4),nLUMO_A))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_A(Valence_A(j,4),nLUMO_A)
    elseif (Valence_A(j,2)==4.and.(ABS(NAOMOm_A(Valence_A(j,4),nHOMO_A))>=cutoff_val).or.(ABS(NAOMOm_A(Valence_A(j,4),nLUMO_A))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_A(Valence_A(j,4),nLUMO_A)
! Cutoff minipatch
    else 
      j_trunc=j_trunc-1
    endif
    j_trunc=j_trunc+1
  enddo
  ! HOMO B
  j_trunc=1
  write(12,'(A2)') 'hB'
  do j=1,numVAL_B
    if (Valence_B(j,2)==2.and.(ABS(NAOMOm_B(Valence_B(j,4),nHOMO_B))>=cutoff_val).or.(ABS(NAOMOm_B(Valence_B(j,4),nLUMO_B))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_B(Valence_B(j,4),nHOMO_B)
    elseif (Valence_B(j,2)==3.and.(ABS(NAOMOm_B(Valence_B(j,4),nHOMO_B))>=cutoff_val).or.(ABS(NAOMOm_B(Valence_B(j,4),nLUMO_B))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_B(Valence_B(j,4),nHOMO_B)
    elseif (Valence_B(j,2)==4.and.(ABS(NAOMOm_B(Valence_B(j,4),nHOMO_B))>=cutoff_val).or.(ABS(NAOMOm_B(Valence_B(j,4),nLUMO_B))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_B(Valence_B(j,4),nHOMO_B)
! Cutoff minipatch
    else 
      j_trunc=j_trunc-1
    endif
    j_trunc=j_trunc+1
  enddo
  ! LUMO B
  j_trunc=1
  write(12,'(A2)') 'lB'
  do j=1,numVAL_B
    if (Valence_B(j,2)==2.and.(ABS(NAOMOm_B(Valence_B(j,4),nHOMO_B))>=cutoff_val).or.(ABS(NAOMOm_B(Valence_B(j,4),nLUMO_B))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_B(Valence_B(j,4),nLUMO_B)
    elseif (Valence_B(j,2)==3.and.(ABS(NAOMOm_B(Valence_B(j,4),nHOMO_B))>=cutoff_val).or.(ABS(NAOMOm_B(Valence_B(j,4),nLUMO_B))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_B(Valence_B(j,4),nLUMO_B)
    elseif (Valence_B(j,2)==4.and.(ABS(NAOMOm_B(Valence_B(j,4),nHOMO_B))>=cutoff_val).or.(ABS(NAOMOm_B(Valence_B(j,4),nLUMO_B))>=cutoff_val)) then
      write(12,'(I3,3X,F9.6)')j_trunc,NAOMOm_B(Valence_B(j,4),nLUMO_B)
! Cutoff minipatch
    else 
      j_trunc=j_trunc-1
    endif
    j_trunc=j_trunc+1
  enddo
  ! Write dE(CA),dE(AC)
  write(12,'(A6,F12.7)') 'EAC : ',eac
  write(12,'(A6,F12.7)') 'ECA : ',eca
  write(12,'(A,F12.7)')  'dE(CT) :', E_nine
  write(12,'(A,F12.7)')  'Reduce :', Reduce
  write(12,'(A,F12.7)')  'Lambda :', Lambda
  ! Cloes simput.out
  close(12)

END PROGRAM Simput

