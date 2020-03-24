Program AutoOpt
! version 2.x     May 12, 2017
  implicit none

  logical,external                :: OK
  character(120),allocatable      :: CalcData(:),OptData(:)
  character(120)                  :: scrname,scrinp,scrlog
  character(256)                  :: command
  character(120)                  :: Simple_input,Structure_input,line
  integer                         :: nCalc
  integer                         :: i,j,n
  real(8)                         :: g(6)
  character(3)                    :: lab(6)

! Printout
  character(5)                    :: method
  real(8)                         :: S
  integer                         :: nCycle
  real(8),allocatable             :: F(:),Tz(:),Ty(:),Tx(:),Rz(:),Ry(:),Rx(:)
  real(8)                         :: sec,omp_get_wtime
  logical                         :: okConst

! Processors
  integer                         :: nProc,Proc,LUN,LUN1,LUN2
  integer                         :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
! Sorting
  logical                         :: octant

  data lab /'T Z','T Y','T X','R Z','R Y','R X'/

  sec = omp_get_wtime()

! Read command line options
  i=iargc()
  if (i < 2) then
    write(6,'(a)')'Usage: AutoOpt Simple_input_file       structures_input_file        [number of processors]'
    write(6,'(a)')'   structures_input_file is output from Simple, SimplePar, or Reduce programs'
    call exit(8)
  endif
  call getarg(1,Simple_input)
  call getarg(2,Structure_input)
  nProc = 1
  if (i > 2) then
    call getarg(3,line)
    read(line,*)nProc
  endif

! Set number of processors
  write(6,'(a,i4)')'Number of processors =',nProc
  call OMP_SET_NUM_THREADS(nProc)

! Set number of optimization calculations and copy data into CalcData
  nCalc=0
  open(8,file=Structure_input,status='old')
1 read(8,*,end=2)
    nCalc=nCalc+1
  goto 1
2 continue
  nCalc=nCalc-1
  rewind(8)
  read(8,*)
  allocate (CalcData(nCalc),OptData(nCalc))
  allocate (F(nCalc),Tz(nCalc),Ty(nCalc),Tx(nCalc),Rz(nCalc),Ry(nCalc),Rx(nCalc))
  do i=1,nCalc
    read(8,'(a)')CalcData(i)
  enddo
  close(8)
! EXPERIMENTAL - for octant
  open(3,file=Simple_input,status='old')
  call ReadOctant(octant)
  close(3)

!!! Run optimization !!!

!$OMP PARALLEL &
  !$OMP SHARED      (CalcData,OptData,Simple_input,nProc,lab)  &
  !$OMP PRIVATE     (Proc,LUN,LUN1,LUN2,scrname,command,line,scrinp,scrlog,g,j,S,nCycle,okConst)

  !$OMP DO PRIVATE(i) SCHEDULE(DYNAMIC,1)

  do i=1,nCalc

    read(CalcData(i),'(21x,6f13.7)')g

!   Open a scratch files and copy the Simple_input into it
    Proc=OMP_GET_THREAD_NUM()
    write(6,'(a,i6,a,i3)')'Optimizing point',i,' on processor',Proc
    LUN=Proc+10
    open(LUN,status='scratch')
    inquire(LUN,name=scrname)
    LUN1=Proc+nProc+10
    open(LUN1,status='scratch')
    inquire(LUN1,name=scrinp)
    LUN2=Proc+2*nProc+10
    open(LUN2,status='scratch')
    inquire(LUN2,name=scrlog)
    command='cp '//Simple_input(1:len_trim(Simple_input))//' '//scrname(1:len_trim(scrname))
    call System(command)
    rewind(LUN)

!   Modify the input inserting CaclData into it
    okConst=.false.
3   read(LUN,'(a)',end=4)line
      write(LUN1,'(a)')line(1:len_trim(line))
      if (index(line,'Const :')>0)then
        okConst=.true.
        do j=1,6
          write(LUN1,'(a3,f13.7)')lab(j),g(j)
        enddo
5       read(LUN,'(a)',end=4)line
          if (index(line,':')>0) then
            backspace(LUN)
            goto 3
          endif
        goto 5
      endif
    goto 3
4   continue
    if (.not. okConst) then
      write(LUN1,'(a)')'Const :'
      do j=1,6
        write(LUN1,'(a3,f13.7)')lab(j),g(j)
      enddo
    endif
    close(LUN)
    command='SimOpt -nm '//scrinp(1:len_trim(scrinp))//' > '//scrlog(1:len_trim(scrlog))
!   Run optimization
    call system(command)
    close(LUN1)
    rewind(LUN2)

!   Read results
6   read(LUN2,'(a)',end=8)line
      if (index(line,'CONVERGED')>0) read(line,'(13x,i5)')nCycle
      if (index(line,'Too many cycles')>0) nCycle=-1
      if (index(line,'Optimized structure:')>0) then
        do j=6,1,-1
          read(LUN2,'(5x,f15.7)',err=8)g(j)
        enddo
        read(LUN2,*,err=8)
        read(LUN2,'(a)',err=8) line        ! method, TA, EREP, k
        read(line,'(a)') method
        if (method.eq.'TA^2=') then
          read(line,'(58x,e15.7)',err=8) S
        else
          read(line,'(52x,e15.7)',err=8) S
        endif
        goto 7
      endif
    goto 6

7   continue
    write(OptData(i),'(i5,e16.7,6f13.7)')i,S,g
    goto 9
8   continue
    write(OptData(i),'(i5,a)')i,' No convergence'
9   continue
    close(LUN2)
    if (nCycle>0) then
      write(6,'(a,i6,a,i5,a,i3)')'......... Finished point',i,' in ',nCycle,' cycles on processor',Proc
    else
      write(6,'(a,i6,a,i3)')'......... Finished point',i,'; no convergence; on processor',Proc
    endif

  enddo

  !$OMP END DO

!$OMP END PARALLEL

  open(8,file='OptResults',status='unknown')
  write(8,'(a)')'              F            T Z          T Y          T X          R Z          R Y          R X'
  do i=1,nCalc
    write(8,'(a)')OptData(i)(1:len_trim(OptData(i)))
  enddo
  close(8)

  if (octant) then
    call system('sort -u -k2,2g -k3,3gr OptResults > OptResultsSorted')
  else
    call system('sort -u -g -k 2 OptResults > OptResultsSorted')
  endif

  open(8,file='OptResultsSorted',status='old')
  open(9,file='OptResultsReduced',status='unknown')
  read(8,'(a)')line
  write(9,'(a)')line(1:len_trim(line))

  n=1
10 read(8,'(i5,e16.7,6f13.7)',end=11)j,F(n),Tz(n),Ty(n),Tx(n),Rz(n),Ry(n),Rx(n)
    n=n+1
  goto 10
11 continue
  n=n-1

  j=1
  do i=1,n
    if (OK(F,Tz,Ty,Tx,i,nCalc)) then
      write(9,'(i5,e16.7,6f13.7)')j,F(i),Tz(i),Ty(i),Tx(i),Rz(i),Ry(i),Rx(i)
      j=j+1
    endif
  enddo
  close(8)
  close(9)

  sec = omp_get_wtime()-sec
  call wtime(sec)

    
End Program AutoOpt
