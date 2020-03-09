Subroutine Shalf(S,Sh)
! version 3.0     July 18, 2019
  implicit none

!     Contruct S^-1/2

  integer,parameter                :: N=2
  integer,parameter                :: lwork=1+6*N+2*N*N
  integer,parameter                :: liwork=3+5*N

  real*8,intent(in)                :: S(2,2)
  real*8,intent(out)               :: Sh(2,2)

  real*8                           :: Eig(N)
  real*8                           :: Work(lwork)
  integer                          :: iWork(liwork)
  integer                          :: k
  integer                          :: INFO

!     Diagonalize the S matrix (it will be destroied in diagonalization)

!     Lapack+Blas diagonalization
!     Diagonalizes symmetrical matrix S, uses upper triangle
!     On leave the eigenvectors are in the original matrix S,
!     eigenvalues in vector Eig

  call DSYEVD('V','U',N,S,N,Eig,Work,lwork,iWork,liwork,INFO)

!     S(i,k) corresponds to the eigenvalue Eig(k)

  if (INFO.ne.0) then
    write(6,*)'ERROR in diagonalization, INFO=',INFO
    call exit(8)
  endif

!     Construct S^-1/2

  do k=1,N
    Eig(k)=sqrt(Eig(k))
    Eig(k)=1.0d0/Eig(k)
  enddo
  call prodAdA(Sh,Eig,S)

  return
End Subroutine Shalf
