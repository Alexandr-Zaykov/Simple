Subroutine GetGrid()
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  real*8                          :: xminA=999.0d0,xmaxA=-999.0d0
  real*8                          :: yminA=999.0d0,ymaxA=-999.0d0
  real*8                          :: zminA=999.0d0,zmaxA=-999.0d0
  real*8                          :: xminB=999.0d0,xmaxB=-999.0d0
  real*8                          :: yminB=999.0d0,ymaxB=-999.0d0
  real*8                          :: zminB=999.0d0,zmaxB=-999.0d0
  integer                         :: i
  real*8                          :: D,step
  real*8                          :: C
  real*8                          :: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
  integer                         :: Xsteps,Ysteps,Zsteps
  real*8                          :: X,Y,Z

  if (fidelity==1) then     !! fine
    D=2.0d0                    !! distance between molecules+rW
    step=0.25d0                !! step size
  endif  
  if (fidelity==2) then     !! medium
    D=1.0d0
    step=0.5d0
  endif
  if (fidelity==3) then     !! coarse
    D=0.0D0
    step=0.75d0
  endif

  do i=1,nA
    if ((xA(1,i)+rWA(i))<xminA) xminA=xA(1,i)+rWA(i)
    if ((xA(1,i)-rWA(i))<xminA) xminA=xA(1,i)-rWA(i)
    if ((xA(1,i)+rWA(i))>xmaxA) xmaxA=xA(1,i)+rWA(i)
    if ((xA(1,i)-rWA(i))>xmaxA) xmaxA=xA(1,i)-rWA(i)

    if ((xA(2,i)+rWA(i))<yminA) yminA=xA(2,i)+rWA(i)
    if ((xA(2,i)-rWA(i))<yminA) yminA=xA(2,i)-rWA(i)
    if ((xA(2,i)+rWA(i))>ymaxA) ymaxA=xA(2,i)+rWA(i)
    if ((xA(2,i)-rWA(i))>ymaxA) ymaxA=xA(2,i)-rWA(i)

    if ((xA(3,i)+rWA(i))<zminA) zminA=xA(3,i)+rWA(i)
    if ((xA(3,i)-rWA(i))<zminA) zminA=xA(3,i)-rWA(i)
    if ((xA(3,i)+rWA(i))>zmaxA) zmaxA=xA(3,i)+rWA(i)
    if ((xA(3,i)-rWA(i))>zmaxA) zmaxA=xA(3,i)-rWA(i)
  enddo

  do i=1,nB
    if ((xB(1,i)+rWB(i))<xminB) xminB=xB(1,i)+rWB(i)
    if ((xB(1,i)-rWB(i))<xminB) xminB=xB(1,i)-rWB(i)
    if ((xB(1,i)+rWB(i))>xmaxB) xmaxB=xB(1,i)+rWB(i)
    if ((xB(1,i)-rWB(i))>xmaxB) xmaxB=xB(1,i)-rWB(i)

    if ((xB(2,i)+rWB(i))<yminB) yminB=xB(2,i)+rWB(i)
    if ((xB(2,i)-rWB(i))<yminB) yminB=xB(2,i)-rWB(i)
    if ((xB(2,i)+rWB(i))>ymaxB) ymaxB=xB(2,i)+rWB(i)
    if ((xB(2,i)-rWB(i))>ymaxB) ymaxB=xB(2,i)-rWB(i)

    if ((xB(3,i)+rWB(i))<zminB) zminB=xB(3,i)+rWB(i)
    if ((xB(3,i)-rWB(i))<zminB) zminB=xB(3,i)-rWB(i)
    if ((xB(3,i)+rWB(i))>zmaxB) zmaxB=xB(3,i)+rWB(i)
    if ((xB(3,i)-rWB(i))>zmaxB) zmaxB=xB(3,i)-rWB(i)
  enddo

  C=max(abs(xminB-xmaxB),abs(yminB-ymaxB),abs(zminB-zmaxB))/2.0d0

  Xmin=xminA-C-D
  Xmax=xmaxA+C+D
  Ymin=yminA-C-D
  Ymax=ymaxA+C+D
  Zmin=zminA-C-D
  Zmax=zmaxA+C+D
  X=Xmax-Xmin
  Y=Ymax-Ymin
  Z=Zmax-Zmin

  Xsteps=(X+0.5)/step   ! Rounding up
  Ysteps=(Y+0.5)/step   ! Rounding up
  Zsteps=(Z+0.5)/step   ! Rounding up
  Xsteps=nint(Xsteps/2.0)
  Ysteps=nint(Ysteps/2.0)
  Zsteps=nint(Zsteps/2.0)
  Xmin=-step*Xsteps
  Ymin=-step*Ysteps
  Zmin=-step*Zsteps
  Xsteps=Xsteps*2+1
  Ysteps=Ysteps*2+1
  Zsteps=Zsteps*2+1

  if (fidelity==1) then       !! fine
    Sini(1)  =  0.0d0
    Sdelta(1)= 10.0d0
    nScan(1) = 36
    Sini(2)  =  0.0d0
    Sdelta(2)= 10.0d0
    nScan(2) = 36
    Sini(3)  =  0.0d0
    Sdelta(3)= 10.0d0
    nScan(3) = 36
  endif
  if (fidelity==2) then       !! medium
    Sini(1)  =  0.0d0
    Sdelta(1)= 15.0d0
    nScan(1) = 24
    Sini(2)  =  0.0d0
    Sdelta(2)= 15.0d0
    nScan(2) = 24
    Sini(3)  =  0.0d0
    Sdelta(3)= 15.0d0
    nScan(3) = 24
  endif
  if (fidelity==3) then       !! coarse
    Sini(1)  =  0.0d0
    Sdelta(1)= 20.0d0
    nScan(1) = 18
    Sini(2)  =  0.0d0
    Sdelta(2)= 20.0d0
    nScan(2) = 18
    Sini(3)  =  0.0d0
    Sdelta(3)= 20.0d0
    nScan(3) = 18
  endif
  Sini(4)  = Xmin
  Sdelta(4)= step
  nScan(4) = Xsteps
  Sini(5)  = Ymin
  Sdelta(5)= step
  nScan(5) = Ysteps
  Sini(6)  = Zmin
  Sdelta(6)= step
  nScan(6) = Zsteps

End Subroutine GetGrid
