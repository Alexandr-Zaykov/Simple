Subroutine BasDim
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rBasDim('a',nA,labA,nNAOA,nOnAtmA,maxnCA,maxnPA)
  call rBasDim('b',nB,labB,nNAOB,nOnAtmB,maxnCB,maxnPB)

End Subroutine BasDim
!------------------------------------------------------------------------------
Subroutine rBasDim(mol,n,lab,nNAO,nOnAtm,maxnC,maxnP)
  implicit none

  character(1),intent(in)         :: mol
  integer,intent(in)              :: n
  character(2),intent(in)         :: lab(n)

  integer,intent(out)             :: nNAO,nOnAtm,maxnC,maxnP
  character(2),allocatable        :: labtmp(:)

  integer                         :: nOA,mxnC,mxnP

  integer                         :: i

  nOnAtm=0
  maxnC=0
  maxnP=0

! Read number of primitives
  rewind(3)

  allocate (labtmp(n))

  labtmp = lab
  do i=1,n
    call lowercase(labtmp(i))
  enddo

  call getNBas(mol,n,labtmp,nNAO,nOA,mxnC,mxnP)

  nOnAtm=max(nOnAtm,nOA)
  maxnC=max(maxnC,mxnC)
  maxnP=max(maxnP,mxnP)

  deallocate (labtmp)

End Subroutine rBasDim
