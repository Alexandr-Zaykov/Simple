Subroutine WriteBasis(mol)
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  character(1),intent(in)         :: mol

  if (mol == 'a') call rWriteBasis(nA,xA,labA,nNAOA,nOnAtmA,OnAtmA,DimnContA,DimnPrimA,maxnCA,TypContA,CContA,ngA,maxnPA,alphaA,gA)
  if (mol == 'b') call rWriteBasis(nB,xB,labB,nNAOB,nOnAtmB,OnAtmB,DimnContB,DimnPrimB,maxnCB,TypContB,CContB,ngB,maxnPB,alphaB,gB)

End Subroutine WriteBasis
!------------------------------------------------------------------------------
Subroutine rWriteBasis(n,x,lab,nNAO,nOnAtm,OnAtm,DimnCont,DimnPrim,maxnC,TypCont,CCont,ng,maxnP,alpha,g)
  implicit none

  integer,intent(in)              :: n
  real*8,intent(in)               :: x(3,n)
  character(2),intent(in)         :: lab(n)
  integer,intent(in)              :: nNAO
  integer,intent(in)              :: nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(in)              :: DimnCont(nNAO)
  integer,intent(in)              :: DimnPrim(nNAO)
  integer,intent(in)              :: maxnC
  character(1),intent(in)         :: TypCont(nNAO,maxnC)
  real*8,intent(in)               :: CCont(nNAO,maxnC,6)
  integer,intent(in)              :: ng(nNAO,maxnC)
  integer,intent(in)              :: maxnP
  real*8,intent(in)               :: alpha(nNAO,maxnP)
  real*8,intent(in)               :: g(nNAO,maxnP)

  integer                         :: i,j,k,l,m,ig,im,ix
  integer                         :: iAO
  character(1)                    :: tmp
  character(1)                    :: S
  character(2)                    :: Px(3)
  character(3)                    :: Dx(6)

  data S  /'s'/
  data Px /'px','py','pz'/
  data Dx /'dxx','dyy','dzz','dxy','dxz','dyz'/

  iAO=0
  do i=1,n
    do j=1,nNAO
      do k=1,OnAtm(j,0)
        if (i == OnAtm(j,k)) then
          write(6,'(i3,2x,a,a,3f12.6,a,i3)')i,lab(i),'[',(x(ix,i),ix=1,3),' ]   NAO:',j
          ig=1
          do l=1,DimnCont(j)
            iAO=iAo+1
            tmp=TypCont(j,l)
            call lowercase(tmp)
            if (tmp=='s') write(6,'(t5,i3,2x,a,6(4x,a,f10.5))')iAO,TypCont(j,l),S,CCont(j,l,1)
            if (tmp=='p') write(6,'(t5,i3,2x,a,6(4x,a,f10.5))')iAO,TypCont(j,l),(Px(m),CCont(j,l,m),m=1,3)
            if (tmp=='d') write(6,'(t5,i3,2x,a,6(4x,a,f10.5))')iAO,TypCont(j,l),(Dx(m),CCont(j,l,m),m=1,6)
            im=0
            do m=ig,ig+ng(j,l)-1
              im=im+1
              write(6,'(t10,i3,2f15.7)')im,alpha(j,m),g(j,m)
            enddo
            ig=ig+ng(j,l)
          enddo
          write(6,*)
        endif
      enddo
    enddo
  enddo

End Subroutine rWriteBasis
