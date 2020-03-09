Subroutine RotOrb(MOs,AOs,nAO,g)
! version 3.0     July 18, 2019
  implicit none

  integer,intent(in)              :: nAO
  real*8,intent(inout)            :: MOs(nAO,2)
  character(3),intent(in)         :: AOs(nAO)
  real*8,intent(in)               :: g(6)

  real*8                          :: d2r
  real*8                          :: a(6)
  integer                         :: b
  integer                         :: n,m,i,j,k,l
  real*8                          :: X(3),y(3)
  real*8                          :: T(3,3)
  real*8                          :: R(9,9)
  real*8                          :: d(9),newd(9)

  d2r = acos(-1.0d0)/180.0d0

  a = g
  do i=1,3
    a(i)=a(i)*d2r
  enddo

! T = Z * Y * X * Rz * Ry * Rx (X' = T * X)

  T(1,1) =  cos(a(2))*cos(a(3))
  T(1,2) = -cos(a(1))*sin(a(3)) + sin(a(1))*sin(a(2))*cos(a(3))
  T(1,3) =  sin(a(1))*sin(a(3)) + cos(a(1))*sin(a(2))*cos(a(3))
  T(2,1) =  cos(a(2))*sin(a(3))
  T(2,2) =  cos(a(1))*cos(a(3)) + sin(a(1))*sin(a(2))*sin(a(3))
  T(2,3) = -sin(a(1))*cos(a(3)) + cos(a(1))*sin(a(2))*sin(a(3))
  T(3,1) = -sin(a(2))
  T(3,2) =  sin(a(1))*cos(a(2))
  T(3,3) =  cos(a(1))*cos(a(2))

! Rotational tensor (rotation of d-orbitals) (tensor product T x T)
! XX=R*(xx,xy,xz,yx,yy,yz,zx,zy,zz); yx=xy,etc.

  n=1
  do k=1,3
    do l=1,3
      m=1
      do i=1,3
        do j=1,3
          R(n,m)=T(k,i)*T(l,j)
          m=m+1
        enddo
      enddo
      n=n+1
    enddo
  enddo

  do b=1,2
    n=1
1   continue
! Rotate p-orbitals
      if (AOs(n) == 'px ') then
        do i=1,3
          X(i)=MOs(n+i-1,b)
        enddo
        y=matmul(T,X)
        do i=1,3
          MOs(n+i-1,b)=y(i)
        enddo
        n=n+3
        if (n > nAO) goto 2
        goto 1
      endif
! Rotate d-orbitals
      if (AOs(n) == 'dxx') then

        d(1)=MOs(n,b)
        d(2)=MOs(n+1,b)
        d(3)=MOs(n+2,b)
        d(4)=d(2)
        d(5)=MOs(n+3,b)
        d(6)=MOs(n+4,b)
        d(7)=d(3)
        d(8)=d(6)
        d(9)=MOs(n+5,b)

        newd=matmul(R,d)

        MOs(n  ,b)=newd(1)
        MOs(n+1,b)=newd(2)
        MOs(n+2,b)=newd(3)
        MOs(n+3,b)=newd(5)
        MOs(n+4,b)=newd(6)
        MOs(n+5,b)=newd(9)

        n=n+6
        if (n > nAO) goto 2
        goto 1
      endif
      n=n+1
      if (n > nAO) goto 2
    goto 1
2   continue
  enddo

  return

End Subroutine RotOrb
