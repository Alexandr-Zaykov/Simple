Subroutine Rate
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  real*8,parameter                :: kT=8.6173303d-5 * 298.0d0
  real*8                          :: p1,p2,Bf
  real*8                          :: DESminus,DESplus,TSplussq,TSminussq
  integer                         :: j

  call int2el
  call DSFmatel

  TSplussq=exSF(3,1)*exSF(3,1)
  DESplus=(exSF(3,3)-exSF(1,1))
  TSminussq=exSF(3,2)*exSF(3,2)
  DESminus=(exSF(3,3)-exSF(2,2))

! Boltzmann weigths
  if ((exSF(2,2)-exSF(1,1))>=1) then
    p1 = 1
    p2 = 0
  elseif ((exSF(2,2)-exSF(1,1))<=-1) then
    p1 = 0
    p2 = 1
  else
    Bf=exp((exSF(2,2)-exSF(1,1))/kT)
    p1=Bf/(1.0d0+Bf)
    p2=1.0d0-p1
  endif

! MARCUS THEORY
  call Marcus(lambda_Marcus,p1,p2,DESminus,DESplus,TSplussq,TSminussq,k_rate,endoerg)
  SF=mixLJ*EnLJ-k_rate

End Subroutine Rate
