Subroutine readScan(line,Sini,Sdelta,nScan,n)
! version 3.0     July 18, 2019
  implicit none

  character(*),intent(in)         :: line
  integer,intent(in)              :: n
  real*8,intent(inout)            :: Sini(6),Sdelta(6)
  integer,intent(inout)           :: nScan(6)
  integer                         :: l,k

  do l=1,len_trim(line)
    if (line(l:l) /= ' ') goto 1
  enddo
1 do k=l,len_trim(line)
    if (line(k:k) == ' ') goto 2
  enddo
2 read(line(l:k-1),*)Sini(n)

  do l=k,len_trim(line)
    if (line(l:l) /= ' ') goto 3
  enddo
3 do k=l,len_trim(line)
    if (line(k:k) == ' ') goto 4
  enddo
4 read(line(l:k-1),*)Sdelta(n)

  do l=k,len_trim(line)
    if (line(l:l) /= ' ') goto 5
  enddo
5 do k=l,len_trim(line)+1
    if (line(k:k) == ' ') goto 6
  enddo
6 read(line(l:k-1),*)nScan(n)

  return

End Subroutine readScan
