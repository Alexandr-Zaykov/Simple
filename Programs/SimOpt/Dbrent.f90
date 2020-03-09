Function Dbrent(ax,bx,cx,f,df,TOL,xmin)
! version 3.0     July 18, 2019
  implicit none

  real*8                          :: Dbrent
  real*8,external                 :: df,f
  integer,parameter               :: ITMAX=100
  real*8,parameter                :: ZEPS=1.0d-16
  real*8,parameter                :: ZERO=0.0d0,HALF=0.5d0,ONE=1.0d0,TWO=2.0d0
  real*8                          :: ax,bx,cx,TOL,xmin
  integer                         :: iter
  real*8                          :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde
  real*8                          :: tol1,tol2,u,u1,u2,v,w,x,xm
  logical                         :: ok1,ok2 ! Will be used as flags for whether proposed steps are acceptable or not.

! Given a function f and its derivative function df, and given a bracketing triplet of abscissas
! ax, bx, cx [such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)],
! this routine isolates the minimum to a fractional precision of about TOL using
! a modification of Brent's method that uses derivatives. The abscissa of the minimum is
! returned as xmin, and the minimum function value is returned as dbrent, the returned function value.

  a=min(ax,cx)
  b=max(ax,cx)
  v=bx
  w=v
  x=v
  e=ZERO
  fx=f(x)
  fv=fx
  fw=fx
  dx=df(x) ! All our housekeeping chores are doubled by the necessity of
           ! moving derivative values around as well as function values.
  dv=dx
  dw=dx

  do iter=1,ITMAX
    xm=HALF*(a+b)
    tol1=TOL*abs(x)+ZEPS
    tol2=TWO*tol1
    if(abs(x-xm).le.(tol2-HALF*(b-a))) goto 3
    if(abs(e).gt.tol1) then
      d1=TWO*(b-a) ! Initialize these d's to an out-of-bracket value.
      d2=d1
      if(dw.ne.dx) d1=(w-x)*dx/(dx-dw) ! Secant method with one point.
      if(dv.ne.dx) d2=(v-x)*dx/(dx-dv) ! And the other.
           ! Which of these two estimates of d shall we take? We will insist that they be within
           ! the bracket, and on the side pointed to by the derivative at x:
      u1=x+d1
      u2=x+d2
      ok1=((a-u1)*(u1-b).gt.ZERO).and.(dx*d1.le.ZERO)
      ok2=((a-u2)*(u2-b).gt.ZERO).and.(dx*d2.le.ZERO)
      olde=e ! Movement on the step before last.
      e=d
      if(.not.(ok1.or.ok2))then ! Take only an acceptable d, and if both
                                ! are acceptable, then take the smallest one.
        goto 1
      else if (ok1.and.ok2)then
        if(abs(d1).lt.abs(d2))then
          d=d1
        else
          d=d2
        endif
      else if (ok1)then
        d=d1
      else
        d=d2
      endif
      if(abs(d).gt.abs(HALF*olde))goto 1
      u=x+d
      if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
      goto 2
    endif
1   if(dx.ge.ZERO) then ! Decide which segment by the sign of the derivative.
      e=a-x
    else
      e=b-x
    endif
    d=HALF*e ! Bisect, not golden section.
2   if(abs(d).ge.tol1) then
      u=x+d
      fu=f(u)
    else
      u=x+sign(tol1,d)
      fu=f(u)
      if(fu.gt.fx)goto 3 ! If the minimum step in the downhill direction takes us uphill,
    endif                ! then we are done.

    du=df(u) ! Now all the housekeeping, sigh.
    if(fu.le.fx) then
      if(u.ge.x) then
        a=x
      else
        b=x
      endif
      v=w
      fv=fw
      dv=dw
      w=x
      fw=fx
      dw=dx
      x=u
      fx=fu
      dx=du
    else
      if(u.lt.x) then
        a=u
      else
        b=u
      endif
      if(fu.le.fw .or. w.eq.x) then
        v=w
        fv=fw
        dv=dw
        w=u
        fw=fu
        dw=du
      else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
        v=u
        fv=fu
        dv=du
      endif
    endif
  enddo

  stop 'Dbrent exceeded maximum iterations'

3 xmin=x
  dbrent=fx

  return

End Function Dbrent
