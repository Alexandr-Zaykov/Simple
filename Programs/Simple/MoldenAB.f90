Subroutine MoldenAB(fname,xTMP)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(*)                    :: fname
  real*8,intent(in)               :: xTMP(3,nB)

  integer                         :: i,j,k,l,m,ig,im,ix
  integer,parameter               :: ntab=18
  character(2)                    :: tab(ntab)
  integer,allocatable             :: itabA(:),itabB(:)
  character(2),allocatable        :: labtmpA(:),labtmpB(:)
  character(1)                    :: tmp
  real*8                          :: zero

  data tab /'h',                               'he',  &
            'li','be','b ','c ','n ','o ','f ','ne',  &
            'na','mg','al','si','p ','s ','cl','ar'/

  zero=0.0d0

  allocate (labtmpA(nA),labtmpB(nB),itabA(nA),itabB(nB))
  labtmpA = labA
  labtmpB = labB
  itabA = 0
  itabB = 0
  
  do i=1,nA
    call lowercase(labtmpA(i))
    do j=1,ntab
     if (labtmpA(i) == tab(j)) itabA(i)=j
    enddo
  enddo

  do i=1,nB
    call lowercase(labtmpB(i))
    do j=1,ntab
     if (labtmpB(i) == tab(j)) itabB(i)=j
    enddo
  enddo

  deallocate (labtmpA,labtmpB)


  open(8,file=fname,status='unknown')
  write(8,'(a)')'[Molden Format]'
  write(8,'(a)')'[Atoms] Angs'
  do i=1,nA
    write(8,'(a3,2i3,3f15.7)')labA(i),i,itabA(i),(xA(j,i),j=1,3)
  enddo
  do i=1,nB
    write(8,'(a3,2i3,3f15.7)')labB(i),i+nA,itabB(i),(xTMP(j,i),j=1,3)
  enddo

  deallocate (itabA,itabB)

  write(8,'(a)')'[GTO]'
  do i=1,nA
    write(8,'(i4,a)')i,'  0'
    do j=1,nNAOA
      do k=1,OnAtmA(j,0)
        if (i == OnAtmA(j,k)) then
          ig=1
          do l=1,DimnContA(j)
            tmp=TypContA(j,l)
            call lowercase(tmp)
            if (tmp == 's') write(8,'(a,i4,a)')'S',ngA(j,l),'  1.0'
            if (tmp == 'p') write(8,'(a,i4,a)')'P',ngA(j,l),'  1.0'
            if (tmp == 'd') write(8,'(a,i4,a)')'D',ngA(j,l),'  1.0'
            do m=ig,ig+ngA(j,l)-1
              write(8,'(2f20.10)')alphaA(j,m),gA(j,m)
            enddo
            ig=ig+ngA(j,l)
          enddo
        endif
      enddo
    enddo
    write(8,*)
  enddo
  do i=1,nB
    write(8,'(i4,a)')i+nA,'  0'
    do j=1,nNAOB
      do k=1,OnAtmB(j,0)
        if (i == OnAtmB(j,k)) then
          ig=1
          do l=1,DimnContB(j)
            tmp=TypContB(j,l)
            call lowercase(tmp)
            if (tmp == 's') write(8,'(a,i4,a)')'S',ngB(j,l),'  1.0'
            if (tmp == 'p') write(8,'(a,i4,a)')'P',ngB(j,l),'  1.0'
            if (tmp == 'd') write(8,'(a,i4,a)')'D',ngB(j,l),'  1.0'
            do m=ig,ig+ngB(j,l)-1
              write(8,'(2f20.10)')alphaB(j,m),gB(j,m)
            enddo
            ig=ig+ngB(j,l)
          enddo
        endif
      enddo
    enddo
    write(8,*)
  enddo

  write(8,'(a)')'[MO]'
  do  i=1,2
    write(8,'(a)')'Sym=A'
    write(8,'(a)')'Ene= 0.0'
    write(8,'(a)')'Spin=Alpha'
    if (i==1) write(8,'(a)')'Occup=2.0'
    if (i==2) write(8,'(a)')'Occup=0.0'
    do k=1,nAOA
      write(8,'(i3,f15.7)')k,MOA(k,i)
    enddo
    do k=1,nAOB
      write(8,'(i3,f15.7)')k+nAOA,zero
    enddo
  enddo
  do  i=1,2
    write(8,'(a)')'Sym=A'
    write(8,'(a)')'Ene= 0.0'
    write(8,'(a)')'Spin=Alpha'
    if (i==1) write(8,'(a)')'Occup=2.0'
    if (i==2) write(8,'(a)')'Occup=0.0'
    do k=1,nAOA
      write(8,'(i5,f15.7)')k,zero
    enddo
    do k=1,nAOB
      write(8,'(i5,f15.7)')k+nAOA,MOB(k,i)
    enddo
  enddo
! HOMO-A HOMO-B
  write(8,'(a)')'Sym=A'
  write(8,'(a)')'Ene= 2.0'
  write(8,'(a)')'Spin=Alpha'
  write(8,'(a)')'Occup=4.0'
  do k=1,nAOA
    write(8,'(i5,f15.7)')k,MOA(k,1)
  enddo
  do k=1,nAOB
    write(8,'(i5,f15.7)')k+nAOA,MOB(k,1)
  enddo
! HOMO-A LUMO-B
  write(8,'(a)')'Sym=A'
  write(8,'(a)')'Ene= 2.0'
  write(8,'(a)')'Spin=Alpha'
  write(8,'(a)')'Occup=2.0'
  do k=1,nAOA
    write(8,'(i5,f15.7)')k,MOA(k,1)
  enddo
  do k=1,nAOB
    write(8,'(i5,f15.7)')k+nAOA,MOB(k,2)
  enddo
! LUMO-A HOMO-B
  write(8,'(a)')'Sym=A'
  write(8,'(a)')'Ene= 2.0'
  write(8,'(a)')'Spin=Alpha'
  write(8,'(a)')'Occup=2.0'
  do k=1,nAOA
    write(8,'(i5,f15.7)')k,MOA(k,2)
  enddo
  do k=1,nAOB
    write(8,'(i5,f15.7)')k+nAOA,MOB(k,1)
  enddo
! LUMO-A LUMO-B
  write(8,'(a)')'Sym=A'
  write(8,'(a)')'Ene= 2.0'
  write(8,'(a)')'Spin=Alpha'
  write(8,'(a)')'Occup=0.0'
  do k=1,nAOA
    write(8,'(i5,f15.7)')k,MOA(k,2)
  enddo
  do k=1,nAOB
    write(8,'(i5,f15.7)')k+nAOA,MOB(k,2)
  enddo

  close(8)

End Subroutine MoldenAB
