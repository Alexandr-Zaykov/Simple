Subroutine int2el
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  real*8                          :: qA(MOsA,2),qB(MOsB,2),tqA(nA),tqB(nB)
  integer                         :: n,i,j,k,l,m
  real*8                          :: normA(2),normB(2)
  real*8                          :: d

  qA = CA
  qB = CB

  normA = 0.0d0
  do i=1,MOsA
    normA(1)=normA(1)+qA(i,1)*qA(i,1)
    normA(2)=normA(2)+qA(i,2)*qA(i,2)
  enddo
  normA(1)=1.0d0/sqrt(normA(1))
  normA(2)=1.0d0/sqrt(normA(2))
  do i=1,MOsA
    qA(i,1)=qA(i,1)*normA(1)
    qA(i,2)=qA(i,2)*normA(2)
  enddo
  normB = 0.0d0
  do i=1,MOsB
    normB(1)=normB(1)+qB(i,1)*qB(i,1)
    normB(2)=normB(2)+qB(i,2)*qB(i,2)
  enddo
  normB(1)=1.0d0/sqrt(normB(1))
  normB(2)=1.0d0/sqrt(normB(2))
  do i=1,MOsB
    qB(i,1)=qB(i,1)*normB(1)
    qB(i,2)=qB(i,2)*normB(2)
  enddo
  
  ! Interaction of S1S0 with S0S1
  Vab=0.0d0
  ! Interaction of between D+ and D- in  D+D- and D-D+
  Jab=0.0d0
  do i=1,MOsA
    do j=1,MOsB
      d=sqrt((xA(1,OnAtmA(i,1))-xB(1,OnAtmB(j,1)))**2+(xA(2,OnAtmA(i,1))-xB(2,OnAtmB(j,1)))**2+(xA(3,OnAtmA(i,1))-xB(3,OnAtmB(j,1)))**2)
      Vab=Vab+qA(i,1)*qA(i,2)*qB(j,1)*qB(j,2)/d
      Jab=Jab-qA(i,1)*qA(i,1)*qB(j,2)*qB(j,2)/d
    enddo
  enddo
! Conversion to eV
  Vab=Vab*14.399645d0
  Jab=Jab*14.399645d0

End Subroutine int2el
