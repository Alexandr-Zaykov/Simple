Subroutine Orthogonalize(mol)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(1),intent(in)         :: mol

  if (mol == 'a') call rOrthogonalize(MOA,SaoA,nAOA,moTMPA,SsaveA)
  if (mol == 'b') call rOrthogonalize(MOB,SaoB,nAOB,moTMPB,SsaveB)

End Subroutine Orthogonalize
!------------------------------------------------------------------------------
Subroutine rOrthogonalize(MOs,S,nao,TMP,Ssave)
  implicit none

  integer,intent(in)              :: nao
  real*8,intent(inout)            :: MOs(nao,2)
  real*8,intent(inout)            :: S(nao,nao)
  real*8,intent(inout)            :: TMP(nao,2)

  real*8                          :: Smo(2,2)
  real*8                          :: Sh(2,2)
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
