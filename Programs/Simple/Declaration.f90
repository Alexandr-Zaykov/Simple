module Declaration

! Geometry and input variables
  character(120)                  :: Title
  integer                         :: method
  integer                         :: nA,nB
  character(2),allocatable        :: labA(:),labB(:)
  real*8,allocatable              :: xA(:,:),xB(:,:)
  real*8,allocatable              :: rWA(:),rWB(:)
  integer,allocatable             :: listA(:),listB(:)

! Charge trasfer energies
  real*8                          :: eAC=0.0d0,eCA=0.0d0,eCT=0.0d0

! Repulsion potential
  real*8                          :: MixLJ,EnLJ

! Scan variables
  real*8                          :: g(6),Sini(6),Sdelta(6)
  integer                         :: nScan(6)
  integer*8                       :: nScans
  real*8,allocatable              :: R(:,:,:,:,:,:)
  integer                         :: v1,v2,v3,v4,v5,v6
  character(3)                    :: gs(6)=(/'R X','R Y','R Z','T X','T Y','T Z'/)

! Wolfsberg-Helmholtz formula
  real*8,allocatable              :: HiiA(:),HiiB(:)
  real*8                          :: H(2,2)

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

! 2el integrals
  real*8                          :: Vab,Jab

! Grid
  logical                         :: grid
  integer                         :: fidelity
  real*8                          :: Q(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1)
  integer                         :: q1,q2,q3,q4,q5,q6
  integer                         :: q1v,q2v,q3v,q4v,q5v,q6v
  integer*8                       :: points,Kpoints

! Checking input
  logical                         :: check

! Couplings
  real*8                          :: H1(5,5),exSF(3,3),SFa,SFb,phase(2)

! finding minima
  real*8,allocatable              :: F(:),Tz(:),Ty(:),Tx(:),Rz(:),Ry(:),Rx(:)

! Boltzmann weighting
  real*8,parameter                :: kT=8.6173303d-5 * 298.0d0
  real*8                          :: p1,p2,Bf
  real*8                          :: DeltaESF

! Marcus theory
  real*8                          :: lambda_Marcus,endoerg
  real*8                          :: k_rate
  real*8                          :: DESminus,DESplus,TSplussq,TSminussq
  real*8                          :: Biexc_Binding,Overall_endoerg

! Temporary arrays
  real*8,allocatable              :: moTMPA(:,:),moTMPB(:,:)
  real*8,allocatable              :: SsaveA(:,:),SsaveB(:,:)
  real*8,allocatable              :: gTMP(:,:)

! Reduction of repulsion
  real*8                          :: Reduce

end module Declaration
