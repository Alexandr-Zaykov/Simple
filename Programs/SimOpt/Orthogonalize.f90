Subroutine Orthogonalize(mol,TMP,Ssave)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(1),intent(in)         :: mol
  real*8,intent(inout)            :: TMP(nAOB,2),Ssave(nAOB,nAOB)

  if (mol == 'a') call rOrthogonalize(MOA,SaoA,nAOA,TMP,Ssave)
  if (mol == 'b') call rOrthogonalize(MOB,SaoB,nAOB,TMP,Ssave)

End Subroutine Orthogonalize
!------------------------------------------------------------------------------
Subroutine rOrthogonalize(MOs,S,nao,TMP,Ssave)
  implicit none

  integer,parameter               :: nmo=2
  integer,intent(in)              :: nao
  real*8,intent(inout)            :: MOs(nao,nmo)
  real*8,intent(inout)            :: S(nao,nao)
  real*8,intent(inout)            :: TMP(nao,nmo)

  real*8                          :: Smo(nmo,nmo)
  real*8                          :: Sh(nmo,nmo)
  real*8,intent(inout)            :: Ssave(nao,nao)

! Smo (in MO basis)
  call TrafoS(Smo,S,nao,MOs)

  Ssave = S
  call Shalf(Smo,Sh)
  TMP = matmul(MOs,Sh)

  S = Ssave
  call TrafoS(Smo,S,nao,TMP)

  MOs = TMP

  return

End Subroutine rOrthogonalize
