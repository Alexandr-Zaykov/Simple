Function LJ
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  real*8                          :: LJ

  integer,parameter               :: nElm=18
  real*8                          :: Re(nElm),eps(nElm)
  real*8                          :: C6(nElm,nElm),C12(nElm,nElm)

  integer                         :: i,j,k
  real*8                          :: d,d6,d12
  real*8,parameter                :: stosest = 106.0d0
  real*8                          :: C6i,C12i,C6j,C12j,Re6
  real*4                          :: x

! Parameters from Gromos (H,C,N,O,F,Na,Si,P,S,Cl) and UFF (He,Li,Be,B,Ne,Mg,Al,Ar)
! Gromos: J. Comput. Chem. 25, 1656-1676, 2004
! UFF:    J. Am. Chem. Soc. 114, 10024-10035, 1992

!                 H                                                                                               He
  data Re /  2.66406d0,                                                                                        2.36200d0, &
!                Li            Be             B             C             N             O             F           Ne
             2.45100d0,    2.74500d0,    4.08300d0,    4.19935d0,    3.52057d0,    3.55332d0,    3.30011d0,    3.24300d0, &
!                Na            Mg            Al            Si             P             S            Cl           Ar
             2.89074d0,    3.02100d0,    4.49900d0,    3.80017d0,    3.80017d0,    3.71276d0,    3.89474d0,    4.49900d0 /

!                 H                                                                                               He
  data eps/  0.02829d0,                                                                                        0.05600d0, &
!                Li            Be             B             C             N             O             F           Ne
             0.02500d0,    0.08500d0,    0.18000d0,    0.12014d0,    0.15291d0,    0.15539d0,    0.10897d0,    0.04200d0, &
!                Na            Mg            Al            Si             P             S            Cl           Ar
             0.01476d0,    0.11100d0,    0.50500d0,    0.58479d0,    0.58479d0,    0.45551d0,    0.30009d0,    0.18500d0 /

  LJ=0.0d0
  if (method>3) then  ! only repulsion for mixed search function
!   Best values for Reduce (1,...,32)
    do i=1,nA
      do j=1,nB
        d=0.0d0
        do k=1,3
          d=d+(xA(k,i)-xB(k,j))*(xA(k,i)-xB(k,j))
        enddo
        d=sqrt(d)
        LJ=LJ+exp(-Reduce*2.0d0)*exp((-244.07402d0-Reduce)*(d/(rWA(i)+rWB(j)))**2)
      enddo
    enddo
    LJ=LJ * 10**(stosest*((45-Reduce)/45)) * ((45-Reduce)/45)**10 ! in eV**2
  else
!   6-12 potential (in kcal mol^-1 A^6, kcal mol^-1 A^12)
    do i=1,nElm
      Re6=Re(i)*Re(i)*Re(i)
      Re6=Re6*Re6
      C6i=2.0d0*eps(i)*Re6
      C12i=eps(i)*Re6*Re6
      do j=1,nElm
        Re6=Re(j)*Re(j)*Re(j)
        Re6=Re6*Re6
        C6j=2.0d0*eps(j)*Re6
        C12j=eps(j)*Re6*Re6
        C6(i,j)=sqrt(C6i*C6j)
        C12(i,j)=sqrt(C12i*C12j)
      enddo
    enddo

    do i=1,nA
      do j=1,nB
        d=0.0d0
        do k=1,3
          d=d+(xA(k,i)-xB(k,j))*(xA(k,i)-xB(k,j))
        enddo
        if (d==0.0d0)d=d+tiny(x)
        d6=d*d*d
        d12=d6*d6
        LJ=LJ+C12(listA(i),listB(j))/d12-C6(listA(i),listB(j))/d6
      enddo
    enddo
  endif
  return

End Function LJ
