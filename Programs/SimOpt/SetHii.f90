Subroutine SetHii
! version 3.0     July 18, 2019
  use Declaration
  implicit none

  call rSetHii(nA,nNAOA,nOnAtmA,OnAtmA,DimnContA,maxnCA,TypContA,listA,nAOA,HiiA)
  call rSetHii(nB,nNAOB,nOnAtmB,OnAtmB,DimnContB,maxnCB,TypContB,listB,nAOB,HiiB)

End Subroutine SetHii
!------------------------------------------------------------------------------
Subroutine rSetHii(n,nNAO,nOnAtm,OnAtm,DimnCont,maxnC,TypCont,list,nAO,Hii)
  implicit none

  integer,intent(in)              :: n
  integer,intent(in)              :: nNAO
  integer,intent(in)              :: nOnAtm
  integer,intent(in)              :: OnAtm(nNAO,0:nOnAtm)
  integer,intent(in)              :: DimnCont(nNAO)
  integer,intent(in)              :: maxnC
  character(1),intent(in)         :: TypCont(nNAO,maxnC)
  integer,intent(in)              :: list(n)
  integer,intent(in)              :: nAO

  real*8,intent(out)              :: Hii(nAO)

  integer                         :: i,j,k,l,m
  integer                         :: iAO
  integer                         :: ic,nn
  character(1)                    :: tmp
  real*8                          :: Hs(18),Hp(18),Hd(18)
  real*8                          :: Hmm

!               H                                                     He
  data Hs /  -13.6,                                                  0.0, &
!              Li      Be       B       C       N       O       F     Ne
              -5.4, -10.0,  -15.2,  -21.4,  -26.0,  -32.0,  -40.0,   0.0, &
!              Na      Mg      Al      Si       P       S      Cl 
              -5.1,  -9.0,  -12.3,  -17.3,  -18.6,  -20.0,  -30.0,   0.0 /
!               H                                                     He
  data Hp /    0.0,                                                  0.0, &
!              Li      Be       B       C       N       O       F     Ne
              -3.5,  -6.0,   -8.5,  -11.4,  -13.4,  -14.8,  -18.1,   0.0, &
!              Na      Mg      Al      Si       P       S      Cl 
              -3.0,  -4.5,   -6.5,   -9.2,  -14.0,  -13.3,  -15.0,   0.0 /

  data Hd /    0.0,                                                  0.0, &
!              Li      Be       B       C       N       O       F     Ne
               0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0, &
!              Na      Mg      Al      Si       P       S      Cl 
               0.0,   0.0,    0.0,   -6.0,   -7.0,   -8.0,   -9.0,   0.0 /

  Hii = 0.0d0

  do m=1,2
    iAO=0
    do i=1,n
      if (Hs(list(i))==0.0) goto 100
      do j=1,nNAO
        do k=1,OnAtm(j,0)
          if (i == OnAtm(j,k)) then
            do l=1,DimnCont(j)
              tmp=TypCont(j,l)
              call lowercase(tmp)
              if (tmp=='s') then
                nn=1
                Hmm=Hs(list(i))
              endif
              if (tmp=='p') then
                nn=3
                Hmm=Hp(list(i))
              endif
              if (tmp=='d') then
                nn=6
                Hmm=Hd(list(i))
              endif
              do ic=1,nn
                iAO=iAO+1
                Hii(iAO)=Hmm
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
  enddo

  return

100 write(6,'()')'Hii parameters for atom number',list(i),' not defined'
  call exit(8)

End Subroutine rSetHii
