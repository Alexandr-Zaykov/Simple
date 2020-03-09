Subroutine RotEig(A,Eig)
! version 3.0     July 18, 2019
  implicit none

  real*8,intent(inout)            :: A(4,4),Eig(4)
  real*8                          :: tmp

! Convert the eigenvector matrix into form:
!      | a -c      |
!      | b  d      |
!      |      u -x |
!      |      v  y |

  if((A(1,1)*A(2,1))<0.0d0) then
    call SwitchCol(4,1,2,A)
    tmp=Eig(2)
    Eig(2)=Eig(1)
    Eig(1)=tmp
  endif
  if(A(1,1)<0.0d0) then
    call NegCol(4,1,A)
  endif
  if(A(1,2)>0.0d0) then
    call NegCol(4,2,A)
  endif
  if((A(3,3)*A(3,4))<0.0d0) then
    call SwitchCol(4,3,4,A)
    tmp=Eig(4)
    Eig(4)=Eig(3)
    Eig(3)=tmp
  endif
  if(A(3,3)<0.0d0) then
    call NegCol(4,3,A)
  endif
  if(A(3,4)>0.0d0) then
    call NegCol(4,4,A)
  endif

End Subroutine RotEig
