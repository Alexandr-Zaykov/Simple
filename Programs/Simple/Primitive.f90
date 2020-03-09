Subroutine Primitive
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rPrimitive(nA,nNAOA,nOnAtmA,OnAtmA,DimnContA,maxnCA,TypContA,CContA,ngA,maxnPA,alphaA,gA,nAOA,a,nPGA,aPrimA,alphaPrimA,gPrimA,CPrimA)
  call rPrimitive(nB,nNAOB,nOnAtmB,OnAtmB,DimnContB,maxnCB,TypContB,CContB,ngB,maxnPB,alphaB,gB,nAOB,b,nPGB,aPrimB,alphaPrimB,gPrimB,CPrimB)

End Subroutine Primitive
!------------------------------------------------------------------------------
Subroutine rPrimitive(n,nNAO,nOnAtm,OnAtm,DimnCont,maxnC,TypCont,CCont,ng,maxnP,alpha,g,nAO,a,nPG,aPrim,alphaPrim,gPrim,CPrim)
  implicit none

  integer,intent(in)              :: n
  integer,intent(in)              :: nNAO
  integer,intent(in)              :: nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(in)              :: DimnCont(nNAO)
  integer,intent(in)              :: maxnC
  character(1),intent(in)         :: TypCont(nNAO,maxnC)
  real*8,intent(in)               :: CCont(nNAO,maxnC,6)
  integer,intent(in)              :: ng(nNAO,maxnC)
  integer,intent(in)              :: maxnP
  real*8,intent(in)               :: alpha(nNAO,maxnP)
  real*8,intent(in)               :: g(nNAO,maxnP)
  integer,intent(in)              :: nAO
  integer,intent(in)              :: a(nAO,3)
  integer,intent(in)              :: nPG

  integer,intent(out)             :: aPrim(nPG,3)
  real*8,intent(out)              :: alphaPrim(nPG),gPrim(nPG),CPrim(nPG)
  

  integer                         :: i,j,k,l,m,ig,im
  integer                         :: iPG,iAO
  integer                         :: aa(3)
  integer                         :: ii,jj
  integer                         :: iCo,nCO
  character(1)                    :: tmp

  iPG=0
  iAO=0
  iCo=1
  im=0

  do i=1,n
    do j=1,nNAO
      do k=1,OnATM(j,0)
        if (i == OnATM(j,k)) then
          ig=1
          do l=1,DimnCont(j)
            iAO=iAO+1
            tmp=TypCont(j,l)
            call lowercase(tmp)
            if (tmp=='s') nCo=1
            if (tmp=='p') nCo=3
            if (tmp=='d') nCo=6
            do ii=1,ng(j,l)
              do jj=1,nCo
                 im=im+1
                CPrim(im)=CCont(j,l,jj)
              enddo
            enddo
            do jj=1,nCo
              do m=ig,ig+ng(j,l)-1
                iPG=iPG+1
                do ii=1,3
                  aa(ii)=a(iCo+jj-1,ii)
                  aPrim(iPG,ii)=aa(ii)
                enddo
                alphaPrim(iPG)=alpha(j,m)
                gPrim(iPG)=g(j,m)
              enddo
            enddo
            iCo=iCo+nCo
            ig=ig+ng(j,l)
          enddo
        endif
      enddo
    enddo
  enddo

End Subroutine rPrimitive
