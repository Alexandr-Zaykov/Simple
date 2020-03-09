Subroutine Marcus(lambda, p1, p2, DESminus, DESplus, TSpsq, TSmsq, k, DEint)
! version 3.0     July 18, 2019
  implicit none

  real*8,intent(in)               :: lambda,p1,p2,DESminus,DESplus,TSpsq,TSmsq,DEint
  real*8,intent(out)              :: k
  real*8,parameter                :: fourkBT=8.6173303d-5*298.0d0*4,Pi8=4.0d0*ATAN(1.0d0),hbar=6.582119514d-16 ! all in eV
     
  ! Calculate kp and km as the T^2 multiplied by the probability from Boltzmann statistics (population of the state)
  ! Approximate DOS with Marcus theory
  k  = ((p1*TSpsq*EXP(-((DESplus +lambda+DEint)**2)/(fourkBT*lambda))/SQRT(fourkBt*Pi8*lambda)) + &
        (p2*TSmsq*EXP(-((DESminus+lambda+DEint)**2)/(fourkBT*lambda))/SQRT(fourkBt*Pi8*lambda)))*2.0d0*Pi8
        !p2*TSmsq*EXP(-((DESminus+lambda+DEint)**2)/(fourkBT*lambda))/SQRT(fourkBt*Pi8*lambda)))*2.0d0*Pi8/hbar

  return

End Subroutine Marcus
