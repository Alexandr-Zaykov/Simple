Subroutine MnBrak(ax,bx,cx,fa,fb,fc,func)
! version 3.0     July 18, 2019
  implicit none

  real*8,external                 :: func
  real*8,parameter                :: GOLD=1.618034d0 ! default ratio by which successive intervals are magnified
  real*8,parameter                :: GLIMIT=100.0d0  ! maximum magnification allowed for a parabolic-fit step
  real*8,parameter                :: TINY=1.0d-20
  real*8                          :: ax,bx,cx,fa,fb,fc
  real*8                          :: dum,fu,q,r,u,ulim

! Given a function func, and given distinct initial points ax and bx, this routine searches
! in the downhill direction (defined by the function as evaluated at the initial points) and
! returns new points ax, bx, cx that bracket a minimum of the function. Also returned are
! the function values at the three points, fa, fb, and fc.

  fa=func(ax)
  fb=func(bx)
  if(fb.gt.fa)then                          ! Switch roles of a and b so that we can go downhill in the direction from a to b.
    dum=ax
    ax=bx
    bx=dum
    dum=fb
    fb=fa
    fa=dum
  endif

  cx=bx+GOLD*(bx-ax)                        ! First guess for c.
  fc=func(cx)

1 if(fb.ge.fc)then                          ! "do while": keep returning here until we bracket.
    r=(bx-ax)*(fb-fc)                       ! Compute u by parabolic extrapolation from a, b, c. TINY
                                            ! is used to prevent any possible division by zero.
    q=(bx-cx)*(fb-fa)
    u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
    ulim=bx+GLIMIT*(cx-bx)                  ! We won't go farther than this. Test various possibilities:
    if((bx-u)*(u-cx).gt.0.)then             ! Parabolic u is between b and c: try it.
      fu=func(u)
      if(fu.lt.fc)then                      ! Got a minimum between b and c.
        ax=bx
        fa=fb
        bx=u
        fb=fu
        return
      else if(fu.gt.fb)then                 ! Got a minimum between between a and u.
        cx=u
        fc=fu
        return
      endif
      u=cx+GOLD*(cx-bx)                     ! Parabolic fit was no use. Use default magnification.
      fu=func(u)
    else if((cx-u)*(u-ulim).gt.0.)then      ! Parabolic fit is between c and its allowed limit.
      fu=func(u)
      if(fu.lt.fc)then
        bx=cx
        cx=u
        u=cx+GOLD*(cx-bx)
        fb=fc
        fc=fu
        fu=func(u)
      endif
    else if((u-ulim)*(ulim-cx).ge.0.)then   ! Limit parabolic u to maximum allowed value.
      u=ulim
      fu=func(u)
    else                                    ! Reject parabolic u, use default magnification.
      u=cx+GOLD*(cx-bx)
      fu=func(u)
    endif
    ax=bx                                   ! Eliminate oldest point and continue.
    bx=cx
    cx=u
    fa=fb
    fb=fc
    fc=fu
    goto 1
  endif

  return

End
