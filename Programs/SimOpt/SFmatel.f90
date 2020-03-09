Subroutine SFmatel
  use Declaration
! version 3.0     July 18, 2019
  implicit none

  real*8                          :: LJ
  real*8,parameter                :: fac=sqrt(1.5d0)

  xB = xBsave
  call TrafoGeom(xB,nB,g,gTMP)
  ! Calculate Lennrd-Jones 6-12 potential or Hard-sphere repulsuion
  EnLJ=LJ()
  if (method==2) then
    EnLJ=EnLJ * 4.33641153d-2 ! kcal to eV
    return
  endif
  
  ! Overlap B in new geometry
  call AOoverlap('bb')
  ! Orbitals of B in new geometry
  MOB = MOBsave
  call RotOrb(MOB,AOsB,nAOB,g)
  ! Normalize B
  call Normalize('b')
  ! Orthogonalize B
  call Orthogonalize('b',moTMPB,SsaveB)
  ! AO intemolecular overlap
  call AOoverlap('ab')
  ! Wolfsberg-Helmhotz approximation
  call Hij(MOA,nAOA,MOB,nAOB,SaoAB,HiiA,HiiB,H)
  ! SF matrix element
  ! Ta = <lA|F|lb><lA|F|hB>-<kA|F|hB><hA|F|lB> = F(2,2)F(2,1)-F(1,1)F(1,2)
  SFa=fac*(H(2,2)*H(2,1)/eCA-H(1,1)*H(1,2)/eAC)
  SFb=fac*(H(2,2)*H(1,2)/eCA-H(1,1)*H(2,1)/eAC)

  if (method==1) SF=SFa
  if (method==2) SF=EnLJ
  if (method==3) SF=MixLJ*EnLJ - SFa*SFa

End Subroutine SFmatel
