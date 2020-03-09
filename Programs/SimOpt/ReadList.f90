Subroutine ReadList(line,n,nBas,nlist,list)
! version 3.0     July 18, 2019
  implicit none

  character(*),intent(in)         :: line
  integer,intent(in)              :: n

  integer,intent(out)             :: list(n)
  integer,intent(out)             :: nlist,nBas

  integer,external                :: nonblank,isblank

  character(20)                   :: s
  integer                         :: i,j,k,l
  integer                         :: lstart,lend

  nlist=0

  l=1
1 continue
  i=nonblank(l,line)
  j=isblank(line(i:))+i-1
  call lowercase(line(i:j))
  if (index(line(i:j),'a')>0) goto 2
  if (index(line(i:j),'b')>0) goto 2
  k=index(line(i:j),'-')
  if (k>0) then
    k=k+i-1
    s=line(i:k-1)
    read(s,*)lstart
    s=line(k+1:j)
    read(s,*)lend
    do i=lstart,lend
      nlist=nlist+1
      list(nlist)=i
    enddo
    l=j+1
    goto 1
  else
    nlist=nlist+1
    read(line(i:j),*) list(nlist)
    l=j+1
    goto 1
  endif

2 continue
  l=j+1
  i=nonblank(l,line)
  j=isblank(line(i:))+i-1
  read(line(i:j),*)nBas

  return
End Subroutine ReadList
