module Declaration
! version 3.0     July 18, 2019

! Geometry and input variables
  character(120)                  :: Title
  integer                         :: method
  integer                         :: nA,nB
  character(2),allocatable        :: labA(:),labB(:)
  real*8,allocatable              :: xA(:,:),xB(:,:),xBsave(:,:)
  real*8,allocatable              :: rWA(:),rWB(:)
  integer,allocatable             :: listA(:),listB(:)
  real*8                          :: eAC,eCA,eCT
  real*8                          :: g(6)
  character(3)                    :: gs(6)=(/'R X','R Y','R Z','T X','T Y','T Z'/)

! Repulsion potential
  real*8                          :: MixLJ,EnLJ
  real*8                          :: Reduce

! Matrix elements
  real*8                          :: SF,SFa,SFb
  real*8                          :: H1(5,5),exSF(3,3)

! Wolfsberg-Helmholtz formula
  real*8,allocatable              :: HiiA(:),HiiB(:)
  real*8                          :: H(2,2)

! Optimization
  logical                         :: last=.false.

! Basis variables
  integer                         :: nNAOA,nOnAtmA,maxnCA,maxnPA
  integer                         :: nNAOB,nOnAtmB,maxnCB,maxnPB
  integer,allocatable             :: OnAtmA(:,:),OnAtmB(:,:)
  integer,allocatable             :: DimnContA(:),DimnContB(:)
  integer,allocatable             :: DimnPrimA(:),DimnPrimB(:)
  character(1),allocatable        :: TypContA(:,:),TypContB(:,:)
  real*8,allocatable              :: CContA(:,:,:),CContB(:,:,:)
  integer,allocatable             :: ngA(:,:),ngB(:,:)
  real*8,allocatable              :: alphaA(:,:),alphaB(:,:)
  real*8,allocatable              :: gA(:,:),gB(:,:)
  integer                         :: nAOA,nAOB
  integer                         :: nPGA,nPGB
  character(3),allocatable        :: AOsA(:),AOsB(:)

! Primitive functions
  integer,allocatable             :: aPrimA(:,:),aPrimB(:,:)
  real*8,allocatable              :: alphaPrimA(:),alphaPrimB(:),gPrimA(:),gPrimB(:),CPrimA(:),CPrimB(:)
  integer,allocatable             :: nAOPrimA(:),nAOPrimB(:)
  integer,allocatable             :: nXA(:),nXB(:)
  integer,allocatable             :: a(:,:),b(:,:)

! Overlap
  real*8,allocatable              :: SaoA(:,:),SaoB(:,:),SaoAB(:,:)
  real*8,allocatable              :: GnormA(:),GnormB(:)

! MOs
  real*8,allocatable              :: CA(:,:),CB(:,:)
  real*8,allocatable              :: MOA(:,:),MOB(:,:),MOBsave(:,:)
  integer                         :: MOsA,MOsB

! Temporary arrays
  real*8,allocatable              :: moTMPA(:,:),moTMPB(:,:)
  real*8,allocatable              :: SsaveA(:,:),SsaveB(:,:)
  real*8,allocatable              :: gTMP(:,:)

! 2el integrals
  real*8                          :: Vab,Jab

! Marcus theory
  real*8                          :: lambda_Marcus,endoerg
  real*8                          :: k_rate

end module Declaration
