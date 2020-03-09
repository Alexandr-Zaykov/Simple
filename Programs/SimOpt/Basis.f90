Subroutine Basis
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rBasis('a',nA,labA,nOnAtmA,nNAOA,maxnCA,maxnPA,OnAtmA,DimnContA,DimnPrimA,TypContA,CContA,ngA,alphaA,gA)
  call rBasis('b',nB,labB,nOnAtmB,nNAOB,maxnCB,maxnPB,OnAtmB,DimnContB,DimnPrimB,TypContB,CContB,ngB,alphaB,gB)

End Subroutine Basis
!------------------------------------------------------------------------------
Subroutine rBasis(mol,n,lab,nOnAtm,nNAO,maxnC,maxnP,OnAtm,DimnCont,DimnPrim,TypCont,CCont,ng,alpha,g)
  implicit none

  character(1),intent(in)         :: mol
  integer,intent(in)              :: n
  character(2),intent(in)         :: lab(n)
  integer,intent(in)              :: nOnAtm,nNAO,maxnC,maxnP

  integer,intent(out)             :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(out)             :: DimnCont(nNAO),DimnPrim(nNAO)
  character(1),intent(out)        :: TypCont(nNAO,maxnC)
  real*8,intent(out)              :: CCont(nNAO,maxnC,6)
  integer,intent(out)             :: ng(nNAO,maxnC)
  real*8,intent(out)              :: alpha(nNAO,maxnP)
  real*8,intent(out)              :: g(nNAO,maxnP)

  character(2),allocatable        :: labtmp(:)

  integer                         :: i

  rewind (3)
  allocate (labtmp(n))
  labtmp = lab
  do i=1,n
    call lowercase(labtmp(i))
  enddo
  call ReadOnAtm(mol,OnAtm,nNAO,nOnAtm,n,labtmp)
  deallocate (labtmp)

  rewind (3)
  call ReadBas(mol,n,lab,nNAO,maxnC,maxnP,DimnCont,DimnPrim,TypCont,CCont,ng,alpha,g)

End Subroutine rBasis
