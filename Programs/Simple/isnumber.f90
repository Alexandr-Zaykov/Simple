function isnumber(s)
! version 3.0     July 18, 2019
  implicit none

  logical                         :: isnumber
  character                       :: s

  isnumber=.false.
  if (ichar(s)<65) isnumber=.true.

  return
End function isnumber
