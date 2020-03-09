Subroutine Rate(prt)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  logical,intent(in)              :: prt

  integer                         :: j

  call int2el
  if (prt) then
    write(6,'(a)')'     Fij(eV)  |hB>      |lB>'
    write(6,'(a,2f10.5)')'     <hA|',(H(1,j),j=1,2)
    write(6,'(a,2f10.5)')'     <lA|',(H(2,j),j=1,2)
    write(6,*)
  endif
  call DSFmatel(prt)

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
  DeltaESF=(p1*(exSF(3,3)-exSF(1,1))+p2*(exSF(3,3)-exSF(2,2)))*1.0d3
  Biexc_binding=((-1)*(exSF(3,3) - H1(5,5)))*1.0d3
  Overall_endoerg=Biexc_binding+DeltaESF
! MARCUS THEORY
  call Marcus(lambda_Marcus, p1, p2, DESminus, DESplus, TSplussq, TSminussq, k_rate, endoerg)

End Subroutine Rate
