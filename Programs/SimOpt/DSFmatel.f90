Subroutine DSFmatel
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  real*8                          :: Hs(4,4),Ht(3,3),Singular(3),U(5,5),VT(5,5),Ov(3,3)
  real*8                          :: Vect(5,3),VectO(5,5)
  integer                         :: i,j,k,l
  integer                         :: INFO
  real*8,parameter                :: TINY=1.0d-13

  integer,parameter               :: N=5
  integer,parameter               :: lwork=1+6*N+2*N*N
  integer,parameter               :: liwork=3+5*N
  real*8                          :: Eig(N)
  real*8                          :: Work(lwork)
  integer                         :: iWork(liwork)

  integer,parameter               :: N4=4
  integer,parameter               :: lwork4=1+6*N4+2*N4*N4
  integer,parameter               :: liwork4=3+5*N4
  real*8                          :: Eig4(N4)
  real*8                          :: Work4(lwork4)
  integer                         :: iWork4(liwork4)

  integer,parameter               :: N3=3
  integer,parameter               :: lwork3=1+6*N3+2*N3*N3
  integer,parameter               :: liwork3=3+5*N3
  real*8                          :: Eig3(N3)
  real*8                          :: Work3(lwork3)
  integer                         :: iWork3(liwork3)

  ! Construct ZDO Hamiltonian
  H1=0.0d0
  H1(3,3)=eAC     !AC
  H1(4,4)=eCA     !CA

  H1(2,1)=2.0d0*Vab
  H1(3,1)=H(2,2)
  H1(4,1)=-H(1,1)

  H1(1,2)=H1(2,1)
  H1(3,2)=-H(1,1)
  H1(4,2)=H(2,2)

  H1(1,3)=H1(3,1)
  H1(2,3)=H1(3,2)
  H1(5,3)=Sqrt(1.5d0)*H(2,1)

  H1(1,4)=H1(4,1)
  H1(2,4)=H1(4,2)
  H1(5,4)=Sqrt(1.5d0)*H(1,2)

  H1(3,5)=H1(5,3)
  H1(4,5)=H1(5,4)

! Left Diagonalization
  do i=1,4
    do j=1,4
      Hs(i,j)=H1(i,j)
    enddo
  enddo
  call DSYEVD('V','U',N4,Hs,N4,Eig4,Work4,lwork4,iWork4,liwork4,INFO)

  if (abs(Eig4(1)-Eig4(2))<TINY) then
! Degenered states
    call RotEigDeg(Hs)
  else
! NotDegenered states
    call RotEig(Hs,Eig4)
  endif

! Right Diagonalization
  do i=1,3
    do j=1,3
      Ht(i,j)=H1(i+2,j+2)
    enddo
  enddo
  call DSYEVD('V','U',N3,Ht,N3,Eig3,Work3,lwork3,iWork3,liwork3,INFO)

  ! Overlap Matrix
  Vect=0.0d0
  Vect(1:4,1)=Hs(1:4,1)
  Vect(1:4,2)=Hs(1:4,2)
  Vect(3:5,3)=Ht(1:3,1)
  Ov=MatMul(Transpose(Vect),Vect)

  ! Singular Value Decomposition
  CALL DGESDD('Singular vectors',5,3,Vect,5,Singular,U,5,VT,5,Work,lwork,iWork,INFO)
  VectO=MatMul(U,VT)

  ! Orthogonalized Singlet Fission Matrix
  exSF=0.0d0
  do i=1,3
    do j=1,3
      do k=1,5
        do l=1,5
          exSF(i,j)=exSF(i,j)+VectO(k,i)*VectO(l,j)*H1(k,l)
        enddo
      enddo
    enddo
  enddo

End Subroutine DSFmatel
