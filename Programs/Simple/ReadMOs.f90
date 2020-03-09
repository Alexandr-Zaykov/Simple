Subroutine ReadMOs(mol)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(1),intent(in)         :: mol

  if (mol == 'a') call rReadMOs(CA,MOsA,'A')
  if (mol == 'b') call rReadMOs(CB,MOsB,'B')

End Subroutine ReadMOs
!------------------------------------------------------------------------------
Subroutine rReadMOs(C,nMOs,sub)
  implicit none

  integer,intent(in)              :: nMOs
  character(1),intent(in)         :: sub

  real*8,intent(out)              :: C(nMOs,2)

  character(120)                  :: line
  integer                         :: m,k,l,n

  C = 0.0d0

  rewind(3)

1 read(3,'(a)',end=99)line
    if (line(1:1) == '#') goto 1
    if (index(line,':') > 0) then
      call lowercase(line)
      if (index(line,'mos') > 0) then
2       read(3,'(a)',end=90)line
          if (index(line,':') > 0) goto 90
          if (line(1:1) == '#') goto 2
          if (sub == 'A') then
            if (index(line,'hA') > 0) then
3             read(3,'(a)',end=90)line
                if (line(1:1) == '#') goto 3
                if (index(line,':') > 0) goto 90
                m=index(line,'hA')+index(line,'lA')+index(line,'hB')+index(line,'lB')
                if (m > 0) then
                  backspace(3)
                  goto 2
                endif
                do l=1,len_trim(line)
                  if (line(l:l) /= ' ') goto 4
                enddo
4               do k=l,len_trim(line)
                  if (line(k:k) == ' ') goto 5
                enddo
5               read(line(l:k),*)n
                do l=k,len_trim(line)
                  if (line(l:l) /= ' ') goto 6
                enddo
6               do k=l,len_trim(line)
                  if (line(k:k) == ' ') goto 7
                enddo
7               read(line(l:k),*)C(n,1)
              goto 3
            endif
            if (index(line,'lA') > 0) then
8             read(3,'(a)',end=90)line
                if (line(1:1) == '#') goto 8
                if (index(line,':') > 0) goto 90
                m=index(line,'hA')+index(line,'lA')+index(line,'hB')+index(line,'lB')
                if (m > 0) then
                  backspace(3)
                  goto 2
                endif
                do l=1,len_trim(line)
                  if (line(l:l) /= ' ') goto 9
                enddo
9               do k=l,len_trim(line)
                  if (line(k:k) == ' ') goto 10
                enddo
10              read(line(l:k),*)n
                do l=k,len_trim(line)
                  if (line(l:l) /= ' ') goto 11
                enddo
11              do k=l,len_trim(line)
                  if (line(k:k) == ' ') goto 12
                enddo
12              read(line(l:k),*)C(n,2)
              goto 8
            endif
          endif
          if (sub == 'B') then
            if (index(line,'hB') > 0) then
13            read(3,'(a)',end=90)line
                if (line(1:1) == '#') goto 13
                if (index(line,':') > 0) goto 90
                m=index(line,'hA')+index(line,'lA')+index(line,'hB')+index(line,'lB')
                if (m > 0) then
                  backspace(3)
                  goto 2
                endif
                do l=1,len_trim(line)
                  if (line(l:l) /= ' ') goto 14
                enddo
14              do k=l,len_trim(line)
                  if (line(k:k) == ' ') goto 15
                enddo
15              read(line(l:k),*)n
                do l=k,len_trim(line)
                  if (line(l:l) /= ' ') goto 16
                enddo
16              do k=l,len_trim(line)
                  if (line(k:k) == ' ') goto 17
                enddo
17              read(line(l:k),*)C(n,1)
              goto 13
            endif
            if (index(line,'lB') > 0) then
18            read(3,'(a)',end=90)line
                if (line(1:1) == '#') goto 18
                if (index(line,':') > 0) goto 90
                m=index(line,'hA')+index(line,'lA')+index(line,'hB')+index(line,'lB')
                if (m > 0) then
                  backspace(3)
                  goto 2
                endif
                do l=1,len_trim(line)
                  if (line(l:l) /= ' ') goto 19
                enddo
19              do k=l,len_trim(line)
                  if (line(k:k) == ' ') goto 20
                enddo
20              read(line(l:k),*)n
                do l=k,len_trim(line)
                  if (line(l:l) /= ' ') goto 21
                enddo
21              do k=l,len_trim(line)
                  if (line(k:k) == ' ') goto 22
                enddo
22              read(line(l:k),*)C(n,2)
              goto 18
            endif
          endif
        goto 2
      endif
    endif
  goto 1

90 continue

  return

99 write(6,'(a)')'ReadMOs: Error, MOs not foud'
  call exit(8)

End Subroutine rReadMOs
