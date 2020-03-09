program Simple
! version 3.0     July 18, 2019
  use Declaration
  implicit none

! Function declaration
  logical                         :: IsMin,ReadCheck

  real*8,allocatable              :: xBsave(:,:),xTMP(:,:)
  character(120)                  :: fnameinp
  character(120)                  :: fname
  character(120)                  :: header
  character(120)                  :: scrname
  character(256)                  :: command
  integer                         :: i,j,k,n

  n = iargc()
  if (n /= 1) then
    write(6,'(a)')'usage: Simple input_file_name'
    call exit (8)
  endif
  call getarg(1,fnameinp)
  scrname = trim(fnameinp(1:len_trim(fnameinp)))
  n = scan(scrname,".",BACK=.true.)
  if (n > 0) fname = scrname(1:n)

! Open input file
  open(3,file=fnameinp,status='old',err=9998)

! Read method
  call ReadMethod

! Read Reduce for reduction of repulsion
  call ReadReduce

! Geometry variables
  call GetConst
  call GetScan

! Read number of atoms
  call ReadDim
  allocate (xA(3,nA), xB(3,nB),xBsave(3,nB),xTMP(3,nB))
  allocate (labA(nA), labB(nB))
  allocate (rWA(nA), rWB(nB))

! Read geometry
  call ReadInp
  xBsave = xB

! Read grid
  call ReadGrid
  if (grid) then
    call GetGrid
  endif
  allocate( R(nScan(6),nScan(5),nScan(4),nScan(3),nScan(2),nScan(1)) )

! Read Monomer state energies
  call ReadE

! Read Marcus parameters
  call GetMarcus

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
  write(6,'(a)')'           Program Simple'
  write(6,'(a)')'========================================'
  write(6,'(a)')'Title='//Title(1:len_trim(Title))
  if (method==1) write(6,'(a)')'Calcuation of TRP matrix elements'
  if (method==2) write(6,'(a)')'Calcuation of TRP**2 matrix elements'
  if (method==3) write(6,'(a)')'Calcuation of Lennard-Jones potential'
  if (method==4) write(6,'(a)')'Calcuation of mixed hard-sphere - TRP**2 surface'
  if (method==5) write(6,'(a)')'Calcuation of mixed hard-sphere - k_rate surface'
  if (method>=4) write(6,'(a,es15.7)')'MixLJ =',MixLJ

  write(6,*)
  write(6,'(a)')'Geometry A:'
  do i=1,nA
    write(6,'(a3,3f12.7)')labA(i),(xA(j,i),j=1,3)
  enddo
  write(6,'(a)')'Geometry B:'
  do i=1,nA
    write(6,'(a3,3f12.7)')labB(i),(xB(j,i),j=1,3)
  enddo

! Do the calculations for all selected geometries
  write(6,'(a)')'========================================'
  nScans=nScan(1)*nScan(2)*nScan(3)*nScan(4)*nScan(5)*nScan(6)
  do i=1,6
    write(6,'(i3,2x,a,i7)')i,gs(i),nScan(i)
  enddo
  write(6,'(a,i16)')'Number of calculations:',nScans
  write(6,'(a)')'========================================'
  write(6,*)
  write(6,'(a)')'Mol A:'
  write(6,'(a,i4)')'Number of NAOs                 =',nNAOA
  do i=1,nNAOA
    write(6,'(a,i3,a,100i4)')'NAO:',i,' located on:',(OnAtmA(i,j),j=1,OnAtmA(i,0))
  enddo

! Read if the run is only checking the input
  check=ReadCheck

! Get number of MOs
  call GetNMOs('a')
  if (check) call WriteMOs('a')

! Write basis set
  write(6,'(a)')'BASIS'
  call WriteBasis('a')

  write(6,*)
  write(6,'(a)')'Mol B:'
  write(6,'(a,i4)')'Number of NAOs                 =',nNAOB
  do i=1,nNAOB
    write(6,'(a,i3,a,100i4)')'NAO:',i,' located on:',(OnAtmB(i,j),j=1,OnAtmB(i,0))
  enddo

! Get number of MOs
  call GetNMOs('b')
  if (check) call WriteMOs('b')

! Write basis set
  write(6,'(a)')'BASIS'
  call WriteBasis('b')

! Stop here if it is only test run
  if (check) call exit(0)

! Number of AOs
  call AOdim

! alloacte temporary arrays
  allocate (moTMPA(nAOA,2),moTMPB(nAOB,2))
  allocate (SsaveA(nAOA,nAOA),SsaveB(nAOB,nAOB))
  allocate (gTMP(3,nB))

! Convert function types into integer indices
  allocate (a(nAOA,3), b(nAOB,3))
  call TypIndex

! set Hii using Wolfberg-Helmhols formula (ICON default values)
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
! Normalize A
  call Normalize('a')
! Orthogonalize A
  call Orthogonalize('a')
! Orthogonalization A

!!! Molecule B !!!
  allocate (CB(MOsB,2), MOB(nAOB,2), MOBsave(nAOB,2))
  call ReadMOs('b')
  call ConstructMOs('b')
! Normalize B
  call Normalize('b')
! Orthogonalize B
  call Orthogonalize('b')
! Save orbitals of B
  MOBsave = MOB

! Print geom parameters (if only one calculation is requested)
  if (nScans == 1) then
    do i=1,6
      write(6,'(a,2x,f12.6)')gs(i),g(i)
    enddo
    write(6,*)
    do i=1,3
      do j=1,nB
        xTMP(i,j)=xB(i,j)
      enddo
    enddo
    call TrafoGeom(xTMP,nB,g,gTMP)
    write(6,*)
    write(6,'(a)')'Transformed geometry of B'
    do i=1,nB
      write(6,'(i4,2x,a,3f12.6)')i,labB(i),(xTMP(j,i),j=1,3)
    enddo
    write(6,*)
  endif

!!! Scan

  if (grid) then
    write(6,'(a)')'=========================================================='
    if (fidelity==1) write(6,'(a)')'AUTOMATIC GRID SELECTED : fine'
    if (fidelity==2) write(6,'(a)')'AUTOMATIC GRID SELECTED : medium'
    if (fidelity==3) write(6,'(a)')'AUTOMATIC GRID SELECTED : coarse'
    write(6,'(a,f10.4,a,f10.4,a,f8.3,a,i4)')'R X: from', Sini(1),' to',Sini(1)+(nScan(1)-1)*Sdelta(1),' step=',Sdelta(1),' points:',nScan(1)
    write(6,'(a,f10.4,a,f10.4,a,f8.3,a,i4)')'R Y: from', Sini(2),' to',Sini(2)+(nScan(2)-1)*Sdelta(2),' step=',Sdelta(2),' points:',nScan(2)
    write(6,'(a,f10.4,a,f10.4,a,f8.3,a,i4)')'R Z: from', Sini(3),' to',Sini(3)+(nScan(3)-1)*Sdelta(3),' step=',Sdelta(3),' points:',nScan(3)
    write(6,'(a,f10.4,a,f10.4,a,f8.3,a,i4)')'T X: from', Sini(4),' to',Sini(4)+(nScan(4)-1)*Sdelta(4),' step=',Sdelta(4),' points:',nScan(4)
    write(6,'(a,f10.4,a,f10.4,a,f8.3,a,i4)')'T Y: from', Sini(5),' to',Sini(5)+(nScan(5)-1)*Sdelta(5),' step=',Sdelta(5),' points:',nScan(5)
    write(6,'(a,f10.4,a,f10.4,a,f8.3,a,i4)')'T Z: from', Sini(6),' to',Sini(6)+(nScan(6)-1)*Sdelta(6),' step=',Sdelta(6),' points:',nScan(6)
    points=nScan(1)*nScan(2)*nScan(3)*nScan(4)*nScan(5)*nScan(6)
    Kpoints=points/1000
    if (Kpoints<1000) then
      write(6,'(a,i6,a,i12,a)')'Number of points',Kpoints,'K  (',points,')'
    else
      Kpoints=Kpoints/1000
      if (Kpoints<1000) then
        write(6,'(a,i6,a,i12,a)')'Number of points',Kpoints,'M  (',points,')'
      else
        Kpoints=Kpoints/1000
        if (Kpoints<1000) then
          write(6,'(a,i6,a,i12,a)')'Number of points',Kpoints,'G  (',points,')'
        else
          Kpoints=Kpoints/1000
          write(6,'(a,i6,a,i12,a)')'Number of points',Kpoints,'T  (',points,')'
        endif
      endif
    endif
    write(6,'(a)')'=========================================================='
    write(6,*)
  endif

  do v1=1,nScan(6)                                         ! T Z
    if (nScan(6) > 1) then
      g(6)=Sini(6)+(v1-1)*Sdelta(6)
      write(6,'(a,i4)')'T Z step:',v1
    endif
    do v2=1,nScan(5)                                       ! T Y
      if (nScan(5) > 1) then
        g(5)=Sini(5)+(v2-1)*Sdelta(5)
        write(6,'(a,i4)')'  T Y step:',v2
      endif
      do v3=1,nScan(4)                                     ! T X
        if (nScan(4) > 1) then
          g(4)=Sini(4)+(v3-1)*Sdelta(4)
!         write(6,'(a,i4)')'    T X step:',v3
        endif
        do v4=1,nScan(3)                                   ! R Z
          if (nScan(3) > 1) g(3)=Sini(3)+(v4-1)*Sdelta(3)
          do v5=1,nScan(2)                                 ! R Y
            if (nScan(2) > 1) g(2)=Sini(2)+(v5-1)*Sdelta(2)
            do v6=1,nScan(1)                               ! R X
              if (nScan(1) > 1) g(1)=Sini(1)+(v6-1)*Sdelta(1)
              call SFmatel(xBsave,moTMPB,SsaveB)
              if (method==1) R(v1,v2,v3,v4,v5,v6)=SFa                ! in eV
              if (method==2) R(v1,v2,v3,v4,v5,v6)=SFa*SFa            ! in eV**2
              if (method==3) R(v1,v2,v3,v4,v5,v6)=EnLJ               ! in eV
              if (method==4) R(v1,v2,v3,v4,v5,v6)=MixLJ*EnLJ-SFa*SFa ! in eV**2
              if (method==5) then
                call Rate(.false.)
                R(v1,v2,v3,v4,v5,v6)=mixLJ*EnLJ-k_rate
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  if (nScans == 1) then
! results details for SP calculations
    call Rate(.true.)
    call PrintRate
    call MoldenAB(fname(1:len_trim(fname))//'molden',xTMP) 
    call Json(fname(1:len_trim(fname))//'json',g,R(1,1,1,1,1,1),SFa,Biexc_binding,Overall_endoerg,k_rate)
  endif
!END OF SP CALC

!FIND MINIMA
  if ((method>=3) .and. (nScans>1)) then
    open(8,status='scratch')
    inquire(8,name=scrname)
    if ((nScan(1)<3).or.(nScan(2)<3).or.(nScan(3)<3).or.(nScan(4)<3).or.(nScan(5)<3).or.(nScan(6)<3)) then
      write(6,'(a)')'Limited number of points in one or more directions; no serach for minima'
      write(6,'(a)')'Data are written into file Surf'
      open(4,file='Surf',status='unknown',form='binary')
      do i=1,6
        write(4)Sini(i),Sdelta(i),nScan(i)
      enddo
      write(4) R
write(6,'(f20.7)')R*10**7
      close(4)
    else
      write(8,'(a)')'              F            T Z          T Y          T X          R Z          R Y          R X'
      if (grid) then
        do v1=2,nScan(6)-1                                       ! T Z
          do v2=2,nScan(5)-1                                     ! T Y
            do v3=2,nScan(4)-1                                   ! T X
              do v4=1,nScan(3)
                do v5=1,nScan(2)
                  do v6=1,nScan(1)
                    do q1=-1,1
                      q1v=v1+q1
                      do q2=-1,1
                        q2v=v2+q2
                        do q3=-1,1
                          q3v=v3+q3
                          do q4=-1,1
                            q4v=v4+q4
                            if (q4v<1) q4v=nScan(3)
                            if (q4v>nScan(3)) q4v=1
                            do q5=-1,1
                              q5v=v5+q5
                              if (q5v<1) q5v=nScan(2)
                              if (q5v>nScan(2)) q5v=1
                              do q6=-1,1
                                q6v=v6+q6
                                if (q6v<1) q6v=nScan(1)
                                if (q6v>nScan(1)) q6v=1
                                Q(q1,q2,q3,q4,q5,q6)=R(q1v,q2v,q3v,q4v,q5v,q6v)
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                    if (IsMin(Q)) then
                      g(1)=Sini(1)+(v6-1)*Sdelta(1)
                      g(2)=Sini(2)+(v5-1)*Sdelta(2)
                      g(3)=Sini(3)+(v4-1)*Sdelta(3)
                      g(4)=Sini(4)+(v3-1)*Sdelta(4)
                      g(5)=Sini(5)+(v2-1)*Sdelta(5)
                      g(6)=Sini(6)+(v1-1)*Sdelta(6)
                      write(8,'(5x,e16.7,6f13.7)')R(v1,v2,v3,v4,v5,v6),g(6),g(5),g(4),g(3),g(2),g(1)
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      else
        do v1=2,nScan(6)-1                                       ! T Z
          do v2=2,nScan(5)-1                                     ! T Y
            do v3=2,nScan(4)-1                                   ! T X
              do v4=2,nScan(3)-1                                 ! R Z
                do v5=2,nScan(2)-1                               ! R Y
                  do v6=2,nScan(1)-1                             ! R X
                    if ((R(v1-1,v2,v3,v4,v5,v6)>R(v1,v2,v3,v4,v5,v6)).and.(R(v1+1,v2,v3,v4,v5,v6)>R(v1,v2,v3,v4,v5,v6))) then
                      if ((R(v1,v2-1,v3,v4,v5,v6)>R(v1,v2,v3,v4,v5,v6)).and.(R(v1,v2+1,v3,v4,v5,v6)>R(v1,v2,v3,v4,v5,v6))) then
                        if ((R(v1,v2,v3-1,v4,v5,v6)>R(v1,v2,v3,v4,v5,v6)).and.(R(v1,v2,v3+1,v4,v5,v6)>R(v1,v2,v3,v4,v5,v6))) then
                          if ((R(v1,v2,v3,v4-1,v5,v6)>R(v1,v2,v3,v4,v5,v6)).and.(R(v1,v2,v3,v4+1,v5,v6)>R(v1,v2,v3,v4,v5,v6))) then
                            if ((R(v1,v2,v3,v4,v5-1,v6)>R(v1,v2,v3,v4,v5,v6)).and.(R(v1,v2,v3,v4,v5+1,v6)>R(v1,v2,v3,v4,v5,v6))) then
                              if ((R(v1,v2,v3,v4,v5,v6-1)>R(v1,v2,v3,v4,v5,v6)).and.(R(v1,v2,v3,v4,v5,v6+1)>R(v1,v2,v3,v4,v5,v6))) then
                                g(1)=Sini(1)+(v6-1)*Sdelta(1)
                                g(2)=Sini(2)+(v5-1)*Sdelta(2)
                                g(3)=Sini(3)+(v4-1)*Sdelta(3)
                                g(4)=Sini(4)+(v3-1)*Sdelta(4)
                                g(5)=Sini(5)+(v2-1)*Sdelta(5)
                                g(6)=Sini(6)+(v1-1)*Sdelta(6)
                                write(8,'(5x,es16.7,6f13.7)')R(v1,v2,v3,v4,v5,v6),g(6),g(5),g(4),g(3),g(2),g(1)
                              endif
                            endif
                          endif
                        endif
                      endif
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
      command='sort -u -g -k 1 '//scrname(1:len_trim(scrname))//' > '//fname(1:len_trim(fname))//'Minima3'
      call system(command)
      close(8)

      open(3,file=fname(1:len_trim(fname))//'Minima3',status='old')
      read(3,'(a)')header
      n=1
1     read(3,*,end=2)
        n=n+1
      goto 1
2     continue
      n=n-1
      rewind(3)
      allocate(F(n),Tz(n),Ty(n),Tx(n),Rz(n),Ry(n),Rx(n))
      read(3,*)
      do i=1,n
        read(3,'(5x,es16.7,6f13.7)')F(i),Tz(i),Ty(i),Tx(i),Rz(i),Ry(i),Rx(i)
      enddo
      close(3)
      open(4,file=fname(1:len_trim(fname))//'Minima3',status='unknown')
      write(4,'(a)')header(1:len_trim(header))
      do i=1,n
        write(4,'(i5,es16.7,6f13.7)')i,F(i),Tz(i),Ty(i),Tx(i),Rz(i),Ry(i),Rx(i)
      enddo
      close(4)
      write(6,'(i5,a)')n,' structures written to file '//fname(1:len_trim(fname))//'Minima3'
      deallocate(F,Tz,Ty,Tx,Rz,Ry,Rx)

    endif
  endif

  goto 9999
9998 write(6,'(a,a)')'Error opening the input file ',fnameinp(1:len_trim(fnameinp))
  call exit(8)

9999 continue

  write(6,*)
  write(6,'(a)')'Program Simple finished'

  call Cpu()

End program Simple
