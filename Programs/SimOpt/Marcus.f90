Subroutine Marcus(lambda,p1,p2,DESminus,DESplus,TSpsq,TSmsq,k,DEint)
! version 3.0     July 18, 2019
  implicit none

  real*8,intent(in)               :: lambda,p1,p2,DESminus,DESplus,TSpsq,TSmsq,DEint
  real*8,intent(out)              :: k
  real*8,parameter                :: fourkT=4.0d0*8.6173303d-5*298.0d0
  real*8,parameter                :: pi8=4.0d0*atan(1.0d0)

  ! Calculate kp and km as the T^2 multiplied by the probability from Boltzmann statistics (population of the state)
  ! Approximate DOS with Marcus theory

  k  = ((p1*TSpsq*exp(-((DESplus +lambda+DEint)**2)/(fourkT*lambda))/sqrt(fourkT*pi8*lambda)) + &
        (p2*TSmsq*exp(-((DESminus+lambda+DEint)**2)/(fourkT*lambda))/sqrt(fourkT*pi8*lambda)))*2.0d0*pi8

  return

End Subroutine Marcus
