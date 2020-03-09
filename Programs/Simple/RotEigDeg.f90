Subroutine RotEigDeg(A)
! version 3.0     July 18, 2019
  implicit none

  real*8,intent(inout)            :: A(4,4)
  real*8                          :: alpha,beta
  real*8                          :: R(4,4),C(4,4)

! Convert the eigenvector matrix into form:
!      | a -b      |
!      | b  a      |
!      |      c -d |
!      |      d  c |

  if((A(1,1)*A(2,1))<0.0d0) then
    call SwitchCol(4,1,2,A)
  endif
  if(A(1,1)<0.0d0) then
    call NegCol(4,1,A)
  endif
  if(A(1,2)>0.0d0) then
    call NegCol(4,2,A)
  endif
  if((A(3,3)*A(3,4))<0.0d0) then
    call SwitchCol(4,3,4,A)
  endif
  if(A(3,3)<0.0d0) then
    call NegCol(4,3,A)
  endif
  if(A(3,4)>0.0d0) then
    call NegCol(4,4,A)
  endif

! Rotate the degenered eigenvecrors into form:
!      | a -b      |     | x -x      |
!      | b  a      | >>> | x  x      |
!      |      c -d |     |      y -y |
!      |      d  c |     |      y  y |

  alpha=atan((A(1,1)-A(2,1))/(A(1,1)+A(2,1)))
  beta =atan((A(3,3)-A(4,3))/(A(3,3)+A(4,3)))
  R=0.0d0
  R(1,1)=cos(alpha); R(1,2)=-sin(alpha)
  R(2,1)=sin(alpha); R(2,2)=cos(alpha)
                                 R(3,3)=cos(beta); R(3,4)=-sin(beta)
                                 R(4,3)=sin(beta); R(4,4)=cos(beta)

  ! Rotated eigenvalues => A' = A * RotationMatrix

  C=matmul(A,R)
  A=C

End Subroutine RotEigDeg
