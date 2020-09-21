MODULE nucleation_tstep_mod

CONTAINS

SUBROUTINE nucleation_rate(rhoa,t_seri,pplay,rh,a_xm,b_xm,c_xm,nucl_rate,ntot,x)

  USE aerophys
  USE infotrac
  USE YOMCST, ONLY : RPI, RD

  IMPLICIT NONE

  ! input variables
  REAL rhoa !H2SO4 number density [molecules/cm3]
  REAL t_seri
  REAL pplay
  REAL rh
  REAL a_xm, b_xm, c_xm

  ! output variables
  REAL nucl_rate
  REAL ntot ! total number of molecules in the critical cluster
  REAL x    ! molefraction of H2SO4 in the critical cluster

  ! local variables
  REAL, PARAMETER                               :: k_B=1.3806E-23  ! Boltzmann constant [J/K]
  REAL                                          :: jnuc !nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s)
  REAL                                          :: rc   !radius of the critical cluster in nm 
  REAL VH2SO4mol

  ! call nucleation routine
  CALL binapara(t_seri,rh,rhoa,jnuc,x,ntot,rc)

  IF (ntot < 4.0) THEN
    !set jnuc to collision rate of two H2SO4 molecules (following personal communication of Ulrike Niemeier and Hanna Vehkamäki)
    VH2SO4mol=mH2SO4mol/(1.e-3*(a_xm+t_seri*(b_xm+t_seri*c_xm))) !cm3
    jnuc = rhoa**2. *(3./4.*RPI)**(1./6.) *(12.*k_B*t_seri/mH2SO4mol)**0.5 &
         & *100.*(2.*VH2SO4mol**(1./3.))**2. !1/(cm3s)
    ntot=2.0
    x=1.0
  ENDIF

  ! convert jnuc from particles/cm3/s to kg(H2SO4)/kgA/s
  nucl_rate=jnuc*ntot*x*mH2SO4mol/(pplay/t_seri/RD/1.E6)

END SUBROUTINE nucleation_rate

!--------------------------------------------------------------------------------------------------

SUBROUTINE nucleation_part(nucl_rate,ntot,x,dt,Vbin,tr_seri)

  USE aerophys
  USE infotrac

  IMPLICIT NONE

  ! input variables
  REAL nucl_rate
  REAL ntot ! total number of molecules in the critical cluster
  REAL x    ! molefraction of H2SO4 in the critical cluster
  REAL dt
  REAL Vbin(nbtr_bin)

  ! output variables
  REAL tr_seri(nbtr)

  ! local variables
  INTEGER k
  REAL Vnew
  REAL ff(nbtr_bin)

  ! dry volume of nucleated particle (at reference temperature)
  Vnew=mH2SO4mol*ntot*x/dens_aer_dry

  ! compute distribution factor for particles of intermediate size (from Jacobson 1994, equation 13)
  ff(:)=0.0 
  DO k=1, nbtr_bin 
  ! CK 20160531: bug fix for first bin
    IF (k.LE.(nbtr_bin-1)) THEN 
      IF (Vbin(k).LE.Vnew.AND.Vnew.LT.Vbin(k+1)) THEN
        ff(k)= Vbin(k)/Vnew*(Vbin(k+1)-Vnew)/(Vbin(k+1)-Vbin(k))
      ENDIF
    ENDIF
    IF (k.EQ.1.AND.Vnew.LE.Vbin(k)) THEN
      ff(k)= 1.
    ENDIF
    IF (k.GT.1) THEN 
      IF (Vbin(k-1).LT.Vnew.AND.Vnew.LT.Vbin(k)) THEN
        ff(k)= 1.-ff(k-1)
      ENDIF
    ENDIF
    IF (k.EQ.nbtr_bin.AND.Vnew.GE.Vbin(k)) THEN
      ff(k)= 1.
    ENDIF
  ENDDO
  ! correction of ff for volume conservation
  DO k=1, nbtr_bin         
    ff(k)=ff(k)*Vnew/Vbin(k)
  ENDDO

  ! add nucleated H2SO4 to corresponding aerosol bin
  DO k=1, nbtr_bin
    tr_seri(k+nbtr_sulgas)=tr_seri(k+nbtr_sulgas)+nucl_rate*dt/(mH2SO4mol*ntot*x)*ff(k)
  ENDDO

END SUBROUTINE nucleation_part

!---------------------------------------------------------------------------------------------------

SUBROUTINE binapara(pt,prh,rhoa_in,jnuc,x,ntot,rc)


  !    Fortran 90 subroutine binapara
  !
  !    Calculates parametrized values of nucleation rate,
  !    mole fraction of sulphuric acid
  !    total number of particles, and the radius of the critical cluster
  !    in H2O-H2SO4 system IF temperature, saturatio ratio of water and 
  !    sulfuric acid concentration  are given. 
  !
  !    Copyright (C) 2002 Hanna Vehkamäki
  !
  !    Division of Atmospheric Sciences
  !    Department of Physical Sciences
  !    P.O. Box 64
  !    FIN-00014 University of Helsinki
  !    Finland
  !    
  !    hanna.vehkamaki@helsinki.fi

  IMPLICIT NONE 

  REAL :: pt,t     !temperature in K (190.15-300.15K)
  REAL :: prh,rh    !saturatio ratio of water (0.0001-1)
  REAL,intent(in) :: rhoa_in    !sulfuric acid concentration in 1/cm3 (10^4-10^11 1/cm3)
  REAL,intent(out) :: jnuc    !nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s)
  REAL,intent(out) :: ntot !total number of molecules in the critical cluster (ntot>4)
  REAL,intent(out) :: x    ! molefraction of H2SO4 in the critical cluster       
  REAL,intent(out) :: rc    !radius of the critical cluster in nm     
  REAL :: rhotres    ! treshold concentration of h2so4 (1/cm^3)
                     ! which produces nucleation rate   1/(cm^3 s) as a function of rh and t 
  REAL rhoa

! CK: use intermediate variables to avoid overwriting
  t=pt
  rh=prh
  rhoa=rhoa_in

  IF (t < 190.15) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): temperature < 190.15 K, using 190.15 K (T=',t,'K)'
     t=190.15
  ENDIF

  IF (t > 300.15) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): temperature > 300.15 K, using 300.15 K'
     t=300.15
  ENDIF

  IF (rh < 0.0001) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): saturation ratio of water < 0.0001, using 0.0001'
     rh=0.0001
  ENDIF

  IF (rh > 1.0) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): saturation ratio of water > 1 using 1'
!     print *, 'rh=',rh
     rh=1.0
  ENDIF

  IF (rhoa < 1.e4) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): sulfuric acid concentration < 1e4 1/cm3, using 1e4 1/cm3'
     rhoa=1.e4
  ENDIF

  IF (rhoa > 1.e11) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev,'): sulfuric acid concentration > 1e11 1/cm3, using 1e11 1/cm3'
     rhoa=1.e11
  ENDIF

  x=  0.7409967177282139 - 0.002663785665140117*t + 0.002010478847383187*Log(rh)  &
       & - 0.0001832894131464668*t*Log(rh) + 0.001574072538464286*Log(rh)**2      &
       & - 0.00001790589121766952*t*Log(rh)**2 + 0.0001844027436573778*Log(rh)**3 & 
       & -  1.503452308794887e-6*t*Log(rh)**3 - 0.003499978417957668*Log(rhoa)    &
       & + 0.0000504021689382576*t*Log(rhoa)

  jnuc= 0.1430901615568665 + 2.219563673425199*t - 0.02739106114964264*t**2 +  &
       &  0.00007228107239317088*t**3 + 5.91822263375044/x +                   &
       &  0.1174886643003278*Log(rh) + 0.4625315047693772*t*Log(rh) -          &
       &  0.01180591129059253*t**2*Log(rh) +                                   &
       &  0.0000404196487152575*t**3*Log(rh) + (15.79628615047088*Log(rh))/x - &
       &  0.215553951893509*Log(rh)**2 - 0.0810269192332194*t*Log(rh)**2 +     &
       &  0.001435808434184642*t**2*Log(rh)**2 -                               &
       &  4.775796947178588e-6*t**3*Log(rh)**2 -                               &
       &  (2.912974063702185*Log(rh)**2)/x - 3.588557942822751*Log(rh)**3 +    &
       &  0.04950795302831703*t*Log(rh)**3 -                                   &
       &  0.0002138195118737068*t**2*Log(rh)**3 +                              &
       &  3.108005107949533e-7*t**3*Log(rh)**3 -                               &
       &  (0.02933332747098296*Log(rh)**3)/x +                                 &
       &  1.145983818561277*Log(rhoa) -                                        &
       &  0.6007956227856778*t*Log(rhoa) +                                     &
       &  0.00864244733283759*t**2*Log(rhoa) -                                 &
       &  0.00002289467254710888*t**3*Log(rhoa) -                              &
       &  (8.44984513869014*Log(rhoa))/x +                                     &
       &  2.158548369286559*Log(rh)*Log(rhoa) +                                &
       &  0.0808121412840917*t*Log(rh)*Log(rhoa) -                             &
       &  0.0004073815255395214*t**2*Log(rh)*Log(rhoa) -                       &
       &  4.019572560156515e-7*t**3*Log(rh)*Log(rhoa) +                        &
       &  (0.7213255852557236*Log(rh)*Log(rhoa))/x +                           &
       &  1.62409850488771*Log(rh)**2*Log(rhoa) -                              &
       &  0.01601062035325362*t*Log(rh)**2*Log(rhoa) +                         &
       &  0.00003771238979714162*t**2*Log(rh)**2*Log(rhoa) +                   &
       &  3.217942606371182e-8*t**3*Log(rh)**2*Log(rhoa) -                     &
       &  (0.01132550810022116*Log(rh)**2*Log(rhoa))/x +                       &
       &  9.71681713056504*Log(rhoa)**2 -                                      &
       &  0.1150478558347306*t*Log(rhoa)**2 +                                  &
       &  0.0001570982486038294*t**2*Log(rhoa)**2 +                            &
       &  4.009144680125015e-7*t**3*Log(rhoa)**2 +                             &
       &  (0.7118597859976135*Log(rhoa)**2)/x -                                &
       &  1.056105824379897*Log(rh)*Log(rhoa)**2 +                             &
       &  0.00903377584628419*t*Log(rh)*Log(rhoa)**2 -                         &
       &  0.00001984167387090606*t**2*Log(rh)*Log(rhoa)**2 +                   &
       &  2.460478196482179e-8*t**3*Log(rh)*Log(rhoa)**2 -                     &
       &  (0.05790872906645181*Log(rh)*Log(rhoa)**2)/x -                       &
       &  0.1487119673397459*Log(rhoa)**3 +                                    &
       &  0.002835082097822667*t*Log(rhoa)**3 -                                &
       &  9.24618825471694e-6*t**2*Log(rhoa)**3 +                              &
       &  5.004267665960894e-9*t**3*Log(rhoa)**3 -                             &
       &  (0.01270805101481648*Log(rhoa)**3)/x
  jnuc=exp(jnuc) !1/(cm3s)

  ntot =-0.002954125078716302 - 0.0976834264241286*t + 0.001024847927067835*t**2 - 2.186459697726116e-6*t**3 -    &
       &   0.1017165718716887/x - 0.002050640345231486*Log(rh) - 0.007585041382707174*t*Log(rh) +                 &
       &   0.0001926539658089536*t**2*Log(rh) - 6.70429719683894e-7*t**3*Log(rh) -                                &
       &   (0.2557744774673163*Log(rh))/x + 0.003223076552477191*Log(rh)**2 + 0.000852636632240633*t*Log(rh)**2 - &
       &   0.00001547571354871789*t**2*Log(rh)**2 + 5.666608424980593e-8*t**3*Log(rh)**2 +                        &
       &   (0.03384437400744206*Log(rh)**2)/x + 0.04743226764572505*Log(rh)**3 -                                  &
       &   0.0006251042204583412*t*Log(rh)**3 + 2.650663328519478e-6*t**2*Log(rh)**3 -                            &
       &   3.674710848763778e-9*t**3*Log(rh)**3 - (0.0002672510825259393*Log(rh)**3)/x -                          &
       &   0.01252108546759328*Log(rhoa) + 0.005806550506277202*t*Log(rhoa) -                                     &
       &   0.0001016735312443444*t**2*Log(rhoa) + 2.881946187214505e-7*t**3*Log(rhoa) +                           &
       &   (0.0942243379396279*Log(rhoa))/x - 0.0385459592773097*Log(rh)*Log(rhoa) -                              &
       &   0.0006723156277391984*t*Log(rh)*Log(rhoa) + 2.602884877659698e-6*t**2*Log(rh)*Log(rhoa) +              &
       &   1.194163699688297e-8*t**3*Log(rh)*Log(rhoa) - (0.00851515345806281*Log(rh)*Log(rhoa))/x -              &
       &   0.01837488495738111*Log(rh)**2*Log(rhoa) + 0.0001720723574407498*t*Log(rh)**2*Log(rhoa) -              &
       &   3.717657974086814e-7*t**2*Log(rh)**2*Log(rhoa) -                                                       &
       &   5.148746022615196e-10*t**3*Log(rh)**2*Log(rhoa) +                                                      &
       &   (0.0002686602132926594*Log(rh)**2*Log(rhoa))/x - 0.06199739728812199*Log(rhoa)**2 +                    &
       &   0.000906958053583576*t*Log(rhoa)**2 - 9.11727926129757e-7*t**2*Log(rhoa)**2 -                          &
       &   5.367963396508457e-9*t**3*Log(rhoa)**2 - (0.007742343393937707*Log(rhoa)**2)/x +                       &
       &   0.0121827103101659*Log(rh)*Log(rhoa)**2 - 0.0001066499571188091*t*Log(rh)*Log(rhoa)**2 +               &
       &   2.534598655067518e-7*t**2*Log(rh)*Log(rhoa)**2 -                                                       &
       &   3.635186504599571e-10*t**3*Log(rh)*Log(rhoa)**2 +                                                      &
       &   (0.0006100650851863252*Log(rh)*Log(rhoa)**2)/x + 0.0003201836700403512*Log(rhoa)**3 -                  &
       &   0.0000174761713262546*t*Log(rhoa)**3 + 6.065037668052182e-8*t**2*Log(rhoa)**3 -                        &
       &   1.421771723004557e-11*t**3*Log(rhoa)**3 + (0.0001357509859501723*Log(rhoa)**3)/x
  ntot=exp(ntot)

  rc=exp(-1.6524245+0.42316402*x+0.33466487*log(ntot)) !nm

  IF (jnuc < 1.e-7) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev'): nucleation rate < 1e-7/cm3s, using 0.0/cm3s,'
     jnuc=0.0
  ENDIF

  IF (jnuc > 1.e10) THEN
!     print *,'Warning (ilon=',ilon,'ilev=',ilev, &
!        & '): nucleation rate > 1e10/cm3s, using 1e10/cm3s'
     jnuc=1.e10
  ENDIF

  rhotres=exp( -279.2430007512709 + 11.73439886096903*rh + 22700.92970508331/t &
       & - (1088.644983466801*rh)/t + 1.144362942094912*t                      &
       & - 0.03023314602163684*rh*t - 0.001302541390154324*t**2                &
       & - 6.386965238433532*Log(rh) + (854.980361026715*Log(rh))/t            &
       & + 0.00879662256826497*t*Log(rh)) !1/cm3

  RETURN

END SUBROUTINE binapara

END MODULE nucleation_tstep_mod
