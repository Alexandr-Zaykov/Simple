Subroutine TrafoS(Smo,S,nao,MOs)
! version 3.0     July 18, 2019
  implicit none

  integer,parameter               :: nmo=2
  integer,intent(in)              :: nao
  real*8,intent(out)              :: Smo(nmo,nmo)
  real*8,intent(in)               :: S(nao,nao)
  real*8,intent(in)               :: MOs(nao,nmo)
  integer                         :: a,b,i,j

  do a=1,nmo
    do b=1,nmo
      Smo(a,b)=0.0d0
      do i=1,nao
        do j=1,nao
          Smo(a,b)=Smo(a,b) + MOs(i,a)*MOs(j,b)*S(i,j)
        enddo
      enddo
    enddo
  enddo

  return

End Subroutine TrafoS
