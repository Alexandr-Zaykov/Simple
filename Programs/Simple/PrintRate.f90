Subroutine PrintRate
! version 3.0     March 24, 2020
  use Declaration
  implicit none

  real*8,parameter                :: rad2deg=180.0d0/acos(-1.0d0)

  write(6,'(6x,a)')               'PERTURBATION THEORY - NO MIXING'
  write(6,'(6X,4(a,es17.7),a)')'TA=',SFa,' eV     TA^2=',SFa*SFa,' eV^2       LJ(Rep)=',EnLJ,  &
                             ' eV^2    Mixed LJ-TA^2=',MixLJ*EnLJ-SFa*SFa,' eV^2'
  write(6,'(6X,4(a,es17.7),a)')'TB=',SFb,' eV     TB^2=',SFb*SFb,' eV^2       LJ(Rep)=',EnLJ,  &
                             ' eV^2    Mixed LJ-TB^2=',MixLJ*EnLJ-SFb*SFb,' eV^2'
  
  write(6,'(6x,a)')               '4x4 and 3x3 DIAGONALIZATION - MIXING OF S1S0 and S0S1'
  write(6,'(6X,4(a,es14.7),a)')'T(S+)=',exSF(3,1),' eV     T(S+)^2=',TSplussq,' eV^2'
  write(6,'(6X,4(a,es14.7),a)')'T(S-)=',exSF(3,2),' eV     T(S-)^2=',TSminussq,' eV^2'
  
  write(6,'(6x,a)')               'MARCUS THEORY'
  write(6,'(6X,2(a,es15.7),a)')'k=',k_rate,' s-1    Mixed LJ-k=',MixLJ*EnLJ-k_rate,' s-1'
  
  write(6,*)
  write(6,'(6x,a)')               'DETAILED RESULTS'
  write(6,'(3x,a)')'Excitonic:'
  write(6,'(/,6x,a,f15.6,a)')'(hAlA|hBlB)         =',Vab*1.0d3,' meV'
  write(6,'(/,6x,a,f15.6,a)')'(hAhA|lBlB)         =',Jab*1.0d3,' meV'

  write(6,'(/,3x,a)')'Mixing LE states in |S+>'
  write(6,'(6x,a,f15.6)')    'Mixing(rad)         =',ATAN(1.0d0/phase(1))
  write(6,'(6x,a,f15.6)')    'Mixing(deg)         =',ATAN(1.0d0/phase(1))*rad2deg

  write(6,'(/,3x,a)')'Mixing LE states in |S->'
  write(6,'(6x,a,f15.6)')    'Mixing(rad)         =',ATAN(-1.0d0*phase(2))
  write(6,'(6x,a,f15.6)')    'Mixing(deg)         =',ATAN(-1.0d0*phase(2))*rad2deg

  write(6,'(/,3x,a)')'Davydov Splitting:'
  write(6,'(6x,a,f15.6, a)')    'dE(S- - S+)         =',(exSF(2,2)-exSF(1,1))*1.0d3, ' meV'

  write(6,'(/,6x,a,es19.6,a)') 'T(S+)^2             =',TSplussq*1.0d6, ' meV^2'
  write(6,'(6x,a,es19.6, a)')   'Mixed LJ-Trp(S+)^2  =',(MixLJ*EnLJ-exSF(3,1)*exSF(3,1))*1.0d6, ' meV^2'
  write(6,'(6x,a,f15.6,a)')    'dE(TT - S+)         =',DESplus*1.0d3, ' meV'
  if (phase(1) > 0.0d0) then
    write(6,'(6x,a)')        'S+  phase           =      + (in phase)'
  elseif (phase(1) < 0.0d0) then
    write(6,'(6x,a)')        'S+  phase           =      - (out of phase)'
  endif

  write(6,'(/,6x,a,es19.6,a)') 'T(S-)^2             =',TSminussq*1.0d6, ' meV^2'
  write(6,'(6x,a,es19.6,a)')   'Mixed LJ-Trp(S-)^2  =',(MixLJ*EnLJ-exSF(3,2)*exSF(3,2))*1.0d6, ' meV^2'
  write(6,'(6x,a,f15.6,a)')    'dE(TT - S-)         =',DESminus*1.0d3, ' meV'
  if (phase(2) > 0.0d0) then
    write(6,'(6X,a)')        'S-  phase           =      + (in phase)'
  elseif (phase(2) < 0.0d0) then
    write(6,'(6X,a)')        'S-  phase           =      - (out of phase)'
  endif

  write(6,*)
  write(6,'(/,3x,a)') 'Biexciton Binding Energy:'
  write(6,'(6X,a,f15.6,a)')    'dE(T1T1 - TT)       =', Biexc_binding, ' meV'

  write(6,'(/,3x,a)') 'Endoergicity:'
  write(6,'(6X,a,f15.6,a)')    'dE(process)         =', Overall_endoerg, ' meV'

  write(6,'(/,3x,a)')'Boltzmann weighting:'
  write(6,'(6X,a,f15.6)')    'w(S+)               =',p1
  write(6,'(6X,a,f15.6)')    'w(S-)               =',p2

  write(6,'(/,3x,a)')'Marcus theory:'
  write(6,'(6x,a,f15.6,a)')  'lambda(reorg.)      =',lambda_Marcus,' eV'
  write(6,'(6x,a,f15.6,a)')  'dE(internal)        =',endoerg,' eV'
  write(6,'(6x,a,es19.6,a)') 'Rate const. k       =',k_rate,' s^-1'
  write(6,'(6x,a,f12.3,a)')  'Lifetime tau        =', 1/k_rate*1.0d12, ' ps'
  
  write(6,*)
  write(6,'(6x,a)') 'CSV:'
  write(6,'(13(es16.7,a))') SFa,',', SFb,',', exSF(3,1),',', exSF(3,2),',', k_rate,',', 2*Vab,',', DESplus,',', DESminus,',',&
                          & exSF(2,2)-exSF(1,1),',', Biexc_binding*1.0d-3,',', Overall_endoerg*1.0d-3

End Subroutine PrintRate
