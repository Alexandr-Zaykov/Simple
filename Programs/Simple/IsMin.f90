Function IsMin(Q)
  implicit none

  logical                         :: IsMin
  real*8,intent(in)               :: Q(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1)

  IsMin=.false.

  if ((Q(-1,0,0,0,0,0)>Q(0,0,0,0,0,0)).and.(Q(1,0,0,0,0,0)>Q(0,0,0,0,0,0))) then
    if ((Q(0,-1,0,0,0,0)>Q(0,0,0,0,0,0)).and.(Q(0,1,0,0,0,0)>Q(0,0,0,0,0,0))) then
      if ((Q(0,0,-1,0,0,0)>Q(0,0,0,0,0,0)).and.(Q(0,0,1,0,0,0)>Q(0,0,0,0,0,0))) then
        if ((Q(0,0,0,-1,0,0)>Q(0,0,0,0,0,0)).and.(Q(0,0,0,1,0,0)>Q(0,0,0,0,0,0))) then
          if ((Q(0,0,0,0,-1,0)>Q(0,0,0,0,0,0)).and.(Q(0,0,0,0,1,0)>Q(0,0,0,0,0,0))) then
            if ((Q(0,0,0,0,0,-1)>Q(0,0,0,0,0,0)).and.(Q(0,0,0,0,0,1)>Q(0,0,0,0,0,0))) then
              IsMin=.true.
            endif
          endif
        endif
      endif
    endif
  endif

  return

End Function IsMin
