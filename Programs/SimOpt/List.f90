Subroutine List
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  integer,parameter               :: nElm=18
  character(2)                    :: atmlab(nElm)
  character(2)                    :: tmp

  integer                         :: i,j

  data atmlab /'h',                               'he',  &
               'li','be','b ','c ','n ','o ','f ','ne',  &
               'na','mg','al','si','p ','s ','cl','ar'/
  listA = 0
  listB = 0

  do i=1,nA
    tmp=labA(i)
    call lowercase(tmp)
    do j=1,nElm
      if (tmp==atmlab(j)) then
        listA(i)=j
      endif
    enddo
  enddo
  do i=1,nB
    tmp=labB(i)
    call lowercase(tmp)
    do j=1,nElm
      if (tmp==atmlab(j)) then
        listB(i)=j
      endif
    enddo
  enddo

  do i=1,nA
    if (listA(i) == 0) goto 1
  enddo

  do i=1,nB
    if (listB(i) == 0) goto 2
  enddo

  return

1 continue
  write(6,'(3a)')'Atom type ',labA(i),' not yet supported, add to the List subroutine'
  call exit(8)
2 continue
  write(6,'(3a)')'Atom type ',labB(i),' not yet supported, add to the List subroutine'
  call exit(8)

End Subroutine List
