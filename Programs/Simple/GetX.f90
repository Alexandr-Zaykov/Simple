Subroutine GetX(lab,X,n)
! version 3.0     July 18, 2019
  implicit none

  integer,intent(in)              :: n
  character(*),intent(out)        :: lab(n)
  real*8,intent(out)              :: X(3,n)
  character(120)                  :: line
  integer                         :: i,l,k

  do i=1,n
9   read(3,'(a)')line
    if (line(1:1) == '#') goto 9
    do l=1,len_trim(line)
      if (line(l:l) /= ' ') goto 1
    enddo
1   do k=l,len_trim(line)
      if (line(k:k) == ' ') then
        lab(i)=line(l:k-1)
        goto 2
      endif
    enddo
2   do l=k,len_trim(line)
      if (line(l:l) /= ' ') goto 3
    enddo
3   do k=l,len_trim(line)
      if (line(k:k) == ' ') then
        read(line(l:k-1),*)X(1,i)
        goto 4
      endif
    enddo
4   do l=k,len_trim(line)
      if (line(l:l) /= ' ') goto 5
    enddo
5   do k=l,len_trim(line)
      if (line(k:k) == ' ') then
        read(line(l:k-1),*)X(2,i)
        goto 6
      endif
    enddo
6   do l=k,len_trim(line)
      if (line(l:l) /= ' ') goto 7
    enddo
7   do k=l,len_trim(line)+1
      if (line(k:k) == ' ') then
        read(line(l:k-1),*)X(3,i)
        goto 8
      endif
    enddo
8   continue
  enddo

  return

End Subroutine GetX
