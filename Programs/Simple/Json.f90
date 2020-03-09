Subroutine Json(fname,const,SF,TA,EBB,dE,k_rate)
! Version 2.1a April 2018
  Implicit none

! Input name
  character(*)                    :: fname
! Constants from input file
  real*8,intent(in)               :: const(6)
! Results
  real*8,intent(in)               :: SF,TA,EBB,dE,k_rate

! open(16,file=fname(1:len_trim(fname)),status='unknown')
  open(16,file=fname,status='unknown')
  write(16,'(a)')'['
  write(16,'(10(a,e15.6),a,f15.6,a)')'{"F":',SF,',"z":',const(6),',"y":',const(5),',"x":',const(4),',"rz":',const(3),',"ry":',const(2),',"rx":',const(1), &
                                     ',"TA":',TA,',"EBB":',EBB,',"ESF":', dE,',"log10_k":', log10(k_rate),'}'
  write(16,'(a)')'['
  close(16)

End Subroutine Json
