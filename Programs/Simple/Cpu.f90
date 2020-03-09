Subroutine Cpu
! version 3.0     July 18, 2019
  implicit none

  real                            :: finish
  real                            :: second,minut,hour,day,mm,hh,dd

  call cpu_time(finish)
  second=modulo(finish,60.0)
  mm=(finish-second)/60.0
  minut=modulo(mm,60.0)
  hh=(mm-minut)/60.0
  hour=modulo(hh,60.0)
  dd=(hh-hour)/24.0
  day=modulo(dd,24.0)
  write(6,'(a,3(i4,a),f8.3,a)')'Elapsed CPU time:',int(day),' days',int(hour),' hours',int(minut),' minutes',second,' seconds'

End Subroutine Cpu
