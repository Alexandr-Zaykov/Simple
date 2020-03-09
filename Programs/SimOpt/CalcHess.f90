Subroutine CalcHess(g0,Hinv)
! version 3.0     July 18, 2019
  implicit none

  real*8                          :: g0(6),Hinv(6,6),Hess(6,6)
  real*8                          :: g(6),F0,SF
  real*8                          :: h(6)=(/1.0D-3,1.0D-3,1.0D-3,1.0D-5,1.0D-5,1.0D-5/) ! (3 rotations, 3 translations)

! Hess inversion
  integer                         :: ipiv(6)
  integer                         :: lwork=30
  real*8                          :: work(30)
  integer                         :: info

  integer                         :: i,j

  g = g0
  call MatElm(g,SF)
  F0=SF

  do i=1,6
    do j=i,6
      if (i==j) then
        g = g0
        g(i)=g(i)-h(i)
        call MatElm(g,SF)
        Hess(i,i)=SF-2.0d0*F0

        g = g0
        g(i)=g(i)+h(i)
        call MatElm(g,SF)
        Hess(i,i)=Hess(i,i)+SF
        Hess(i,i)=Hess(i,i)/(h(i)*h(i))
      else
        g = g0
        g(i)=g(i)-h(i)
        g(j)=g(j)-h(j)
        call MatElm(g,SF)
        Hess(i,j)=SF

        g = g0
        g(i)=g(i)+h(i)
        g(j)=g(j)+h(j)
        call MatElm(g,SF)
        Hess(i,j)=Hess(i,j)+SF

        g = g0
        g(i)=g(i)+h(i)
        g(j)=g(j)-h(j)
        call MatElm(g,SF)
        Hess(i,j)=Hess(i,j)-SF

        g = g0
        g(i)=g(i)-h(i)
        g(j)=g(j)+h(j)
        call MatElm(g,SF)
        Hess(i,j)=Hess(i,j)-SF
        Hess(i,j)=Hess(i,j)/(4.0d0*h(i)*h(j))
        Hess(j,i)=Hess(i,j)
      endif
    enddo
  enddo

  Hinv = Hess
  call DGETRF(6,6,Hinv,6,ipiv,info)
  call DGETRI(6,Hinv,6,ipiv,work,lwork,info)

End Subroutine CalcHess
