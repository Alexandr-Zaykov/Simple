Function OK(F,Tz,Ty,Tx,i,n)
! version 1.7     May 12, 2017
  implicit none

  logical                         :: OK

  integer,intent(in)              :: i,n
  real*8,intent(in)               :: F(n),Tz(n),Ty(n),Tx(n)

  integer                         :: j
  logical                         :: OK1,OK2,OK3,OK4,OK5,OK6

  OK=.true.

  do j=1,i-1
    OK1=(abs(F(i)-F(j))<abs(F(i)*1.0d-6))
    OK2=(abs(abs(Tz(i))-abs(Tz(j)))<1.0d-5)
    OK3=(abs(abs(Ty(i))-abs(Ty(j)))<1.0d-5)
    OK4=(abs(abs(Tx(i))-abs(Tx(j)))<1.0d-5)
    if (OK1.and.OK2.and.OK3.and.OK4) then
      OK=.false.
      return
    endif
  enddo

End Function OK
