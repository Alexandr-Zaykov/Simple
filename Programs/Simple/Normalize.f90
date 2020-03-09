Subroutine Normalize(mol)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(1),intent(in)         :: mol

  if (mol == 'a') call rNormalize(MOA,SaoA,nAOA)
  if (mol == 'b') call rNormalize(MOB,SaoB,nAOB)

End Subroutine Normalize
!------------------------------------------------------------------------------
Subroutine rNormalize(MOs,S,nao)
  implicit none

  integer,parameter               :: nmo=2
  integer,intent(in)              :: nao
  real*8,intent(inout)            :: MOs(nao,nmo)
  real*8,intent(in)               :: S(nao,nao)
  real*8                          :: N
  integer                         :: a,i,j
  real*8                          :: Saa

  do a=1,nmo
    Saa=0.0d0
    do i=1,nao
      do j=1,nao
        Saa=Saa + MOs(i,a)*MOs(j,a)*S(i,j)
      enddo
    enddo
    N=sqrt(1.0d0/Saa)
    do i=1,nao
      MOs(i,a)=N*MOs(i,a)
    enddo
  enddo

  return

End Subroutine rNormalize
