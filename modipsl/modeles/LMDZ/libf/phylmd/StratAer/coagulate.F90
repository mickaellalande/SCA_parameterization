SUBROUTINE COAGULATE(pdtcoag,mdw,tr_seri,t_seri,pplay,dens_aer,is_strato)
  !     -----------------------------------------------------------------------
  !
  !     Author : Christoph Kleinschmitt (with Olivier Boucher)
  !     ------
  !
  !     purpose
  !     -------
  !
  !     interface
  !     ---------
  !      input 
  !       pdtphys        time step duration                 [sec]
  !       tr_seri        tracer mixing ratios               [kg/kg]
  !       mdw             # or mass median diameter          [m]
  !
  !     method
  !     ------
  !
  !     -----------------------------------------------------------------------

  USE dimphy, ONLY : klon,klev
  USE aerophys
  USE infotrac
  USE phys_local_var_mod, ONLY: DENSO4, f_r_wet
  USE YOMCST

  IMPLICIT NONE

  !--------------------------------------------------------

  ! transfer variables when calling this routine
  REAL,INTENT(IN)                               :: pdtcoag ! Time step in coagulation routine [s]
  REAL,DIMENSION(nbtr_bin),INTENT(IN)           :: mdw     ! aerosol particle diameter in each bin [m]
  REAL,DIMENSION(klon,klev,nbtr),INTENT(INOUT)  :: tr_seri ! Concentration Traceur [U/KgA]
  REAL,DIMENSION(klon,klev),INTENT(IN)          :: t_seri  ! Temperature
  REAL,DIMENSION(klon,klev),INTENT(IN)          :: pplay   ! pression pour le mileu de chaque couche (en Pa)
  REAL,DIMENSION(klon,klev)                     :: dens_aer! density of aerosol [kg/m3 aerosol] with default H2SO4 mass
  LOGICAL,DIMENSION(klon,klev),INTENT(IN)       :: is_strato

  ! local variables in coagulation routine
  INTEGER                                       :: i,j,k,nb,ilon,ilev
  REAL, DIMENSION(nbtr_bin)                     :: radius ! aerosol particle radius in each bin [m]
  REAL, DIMENSION(klon,klev,nbtr_bin)           :: tr_t ! Concentration Traceur at time t [U/KgA]
  REAL, DIMENSION(klon,klev,nbtr_bin)           :: tr_tp1 ! Concentration Traceur at time t+1 [U/KgA]
  REAL, DIMENSION(nbtr_bin,nbtr_bin,nbtr_bin)   :: ff   ! Volume fraction of intermediate particles
  REAL, DIMENSION(nbtr_bin)                     :: V    ! Volume of bins
  REAL, DIMENSION(nbtr_bin,nbtr_bin)            :: Vij  ! Volume sum of i and j
  REAL                                          :: eta  ! Dynamic viscosity of air
  REAL, PARAMETER                               :: mair=4.8097E-26 ! Average mass of an air molecule [kg]
  REAL                                          :: zrho ! Density of air
  REAL                                          :: mnfrpth ! Mean free path of air
  REAL, DIMENSION(nbtr_bin)                     :: Kn   ! Knudsen number of particle i
  REAL, DIMENSION(nbtr_bin)                     :: Di   ! Particle diffusion coefficient
  REAL, DIMENSION(nbtr_bin)                     :: m_par   ! Mass of particle i
  REAL, DIMENSION(nbtr_bin)                     :: thvelpar! Thermal velocity of particle i
  REAL, DIMENSION(nbtr_bin)                     :: mfppar  ! Mean free path of particle i
  REAL, DIMENSION(nbtr_bin)                     :: delta! delta of particle i (from equation 21)
  REAL, DIMENSION(nbtr_bin,nbtr_bin)            :: beta ! Coagulation kernel from Brownian diffusion
  REAL                                          :: beta_const ! Constant coagulation kernel (for comparison)
  REAL                                          :: num
  REAL                                          :: numi
  REAL                                          :: denom

  ! Additional variables for coagulation enhancement factor due to van der Waals forces
  ! Taken from Chan and Mozurkewich, Measurement of the coagulation rate constant for sulfuric acid 
  ! particles as a function of particle size using TDMA, Aerosol Science, 32, 321-339, 2001.
  !--ok_vdw is 0 for no vdW forces, 1 for E(0), 2 for E(infinity)
  INTEGER, PARAMETER                            :: ok_vdw = 0
  REAL, PARAMETER                               :: avdW1 = 0.0757
  REAL, PARAMETER                               :: avdW3 = 0.0015
  REAL, PARAMETER                               :: bvdW0 = 0.0151
  REAL, PARAMETER                               :: bvdW1 = -0.186
  REAL, PARAMETER                               :: bvdW3 = -0.0163
  REAL, PARAMETER                               :: AvdW = 6.4e-20 !Hamaker constant in J = 1e7 erg
  REAL                                          :: AvdWi
  REAL                                          :: xvdW
  REAL                                          :: EvdW

  DO i=1, nbtr_bin
   radius(i)=mdw(i)/2.
   V(i)= radius(i)**3.  !neglecting factor 4*RPI/3
  ENDDO

  DO j=1, nbtr_bin
  DO i=1, nbtr_bin
   Vij(i,j)= V(i)+V(j)
  ENDDO
  ENDDO

!--pre-compute the f(i,j,k) from Jacobson equation 13
  ff=0.0 
  DO k=1, nbtr_bin 
  DO j=1, nbtr_bin
  DO i=1, nbtr_bin
    IF (k.EQ.1) THEN
      ff(i,j,k)= 0.0
    ELSEIF (k.GT.1.AND.V(k-1).LT.Vij(i,j).AND.Vij(i,j).LT.V(k)) THEN
      ff(i,j,k)= 1.-ff(i,j,k-1)
    ELSEIF (k.EQ.nbtr_bin) THEN
      IF (Vij(i,j).GE.v(k)) THEN
        ff(i,j,k)= 1.
      ELSE
        ff(i,j,k)= 0.0
      ENDIF
    ELSEIF (k.LE.(nbtr_bin-1).AND.V(k).LE.Vij(i,j).AND.Vij(i,j).LT.V(k+1)) THEN
      ff(i,j,k)= V(k)/Vij(i,j)*(V(k+1)-Vij(i,j))/(V(k+1)-V(k))
    ENDIF
  ENDDO
  ENDDO
  ENDDO

  DO ilon=1, klon
  DO ilev=1, klev 
  !only in the stratosphere
  IF (is_strato(ilon,ilev)) THEN
  !compute actual wet particle radius & volume for every grid box
  DO i=1, nbtr_bin
   radius(i)=f_r_wet(ilon,ilev)*mdw(i)/2.
   V(i)= radius(i)**3.  !neglecting factor 4*RPI/3
  ENDDO

!--Calculations for the coagulation kernel---------------------------------------------------------

  zrho=pplay(ilon,ilev)/t_seri(ilon,ilev)/RD

!--initialize the tracer at time t and convert from [number/KgA] to [number/m3]
  DO i=1, nbtr_bin
  tr_t(ilon,ilev,i) = tr_seri(ilon,ilev,i+nbtr_sulgas) * zrho
  ENDDO

  ! mean free path of air (Pruppacher and Klett, 2010, p.417) [m]
  mnfrpth=6.6E-8*(1.01325E+5/pplay(ilon,ilev))*(t_seri(ilon,ilev)/293.15)
  ! mnfrpth=2.*eta/(zrho*thvelair)
  ! mean free path of air (Prupp. Klett) in [10^-6 m]
  ! ZLAIR = 0.066 *(1.01325E+5/PPLAY)*(T_SERI/293.15)*1.E-06

  ! dynamic viscosity of air (Pruppacher and Klett, 2010, p.417) [kg/(m*s)]
  IF (t_seri(ilon,ilev).GE.273.15) THEN
    eta=(1.718+0.0049*(t_seri(ilon,ilev)-273.15))*1.E-5
  ELSE
    eta=(1.718+0.0049*(t_seri(ilon,ilev)-273.15)-1.2E-5*(t_seri(ilon,ilev)-273.15)**2)*1.E-5
  ENDIF

!--pre-compute the particle diffusion coefficient Di(i) from equation 18
  Di=0.0
  DO i=1, nbtr_bin
   Kn(i)=mnfrpth/radius(i)
   Di(i)=RKBOL*t_seri(ilon,ilev)/(6.*RPI*radius(i)*eta)*(1.+Kn(i)*(1.249+0.42*exp(-0.87/Kn(i))))
  ENDDO

!--pre-compute the thermal velocity of a particle thvelpar(i) from equation 20
  thvelpar=0.0
  DO i=1, nbtr_bin
   m_par(i)=4./3.*RPI*radius(i)**3.*DENSO4(ilon,ilev)*1000.
   thvelpar(i)=sqrt(8.*RKBOL*t_seri(ilon,ilev)/(RPI*m_par(i)))
  ENDDO

!--pre-compute the particle mean free path mfppar(i) from equation 22
  mfppar=0.0
  DO i=1, nbtr_bin
   mfppar(i)=8.*Di(i)/(RPI*thvelpar(i))
  ENDDO

!--pre-compute the mean distance delta(i) from the center of a sphere reached by particles 
!--leaving the surface of the sphere and traveling a distance of particle mfppar(i) from equation 21
  delta=0.0
  DO i=1, nbtr_bin
   delta(i)=((2.*radius(i)+mfppar(i))**3.-(4.*radius(i)**2.+mfppar(i)**2.)**1.5)/ &
           & (6.*radius(i)*mfppar(i))-2.*radius(i)
  ENDDO

!--pre-compute the beta(i,j) from equation 17 in Jacobson
  num=0.0
  DO j=1, nbtr_bin
  DO i=1, nbtr_bin
!
   num=4.*RPI*(radius(i)+radius(j))*(Di(i)+Di(j))
   denom=(radius(i)+radius(j))/(radius(i)+radius(j)+sqrt(delta(i)**2.+delta(j)**2.))+ &
        & 4.*(Di(i)+Di(j))/(sqrt(thvelpar(i)**2.+thvelpar(j)**2.)*(radius(i)+radius(j)))
   beta(i,j)=num/denom
!
!--compute enhancement factor due to van der Waals forces
   IF (ok_vdw .EQ. 0) THEN      !--no enhancement factor
     Evdw=1.0
   ELSEIF (ok_vdw .EQ. 1) THEN  !--E(0) case
     AvdWi = AvdW/(RKBOL*t_seri(ilon,ilev))*(4.*radius(i)*radius(j))/(radius(i)+radius(j))**2.
     xvdW = LOG(1.+AvdWi)
     EvdW = 1. + avdW1*xvdW + avdW3*xvdW**3
   ELSEIF (ok_vdw .EQ. 2) THEN  !--E(infinity) case
     AvdWi = AvdW/(RKBOL*t_seri(ilon,ilev))*(4.*radius(i)*radius(j))/(radius(i)+radius(j))**2.
     xvdW = LOG(1.+AvdWi)
     EvdW = 1. + SQRT(AvdWi/3.)/(1.+bvdW0*SQRT(AvdWi)) + bvdW1*xvdW + bvdW3*xvdW**3.
   ENDIF
!
   beta(i,j)=beta(i,j)*EvdW

  ENDDO
  ENDDO

!--external loop for equation 14
  DO k=1, nbtr_bin

!--calculating denominator sum 
  denom=0.0
  DO j=1, nbtr_bin
  denom=denom+(1.-ff(k,j,k))*beta(k,j)*tr_t(ilon,ilev,j)
  ENDDO

  IF (k.EQ.1) THEN
!--calculate new concentration of smallest bin
    tr_tp1(ilon,ilev,k)=tr_t(ilon,ilev,k)/(1.+pdtcoag*denom)
    ELSE
!--calculating double sum terms in numerator of eq 14
    num=0.0
    DO j=1, k
    numi=0.0
    DO i=1, k-1
    numi=numi+ff(i,j,k)*beta(i,j)*V(i)*tr_tp1(ilon,ilev,i)*tr_t(ilon,ilev,j)
    ENDDO
    num=num+numi
    ENDDO

!--calculate new concentration of other bins
    tr_tp1(ilon,ilev,k)=(V(k)*tr_t(ilon,ilev,k)+pdtcoag*num)/(1.+pdtcoag*denom)/V(k)
  ENDIF

  ENDDO !--end of loop k

!--convert tracer concentration back from [number/m3] to [number/KgA] and write into tr_seri
  DO i=1, nbtr_bin
   tr_seri(ilon,ilev,i+nbtr_sulgas) = tr_tp1(ilon,ilev,i) / zrho
  ENDDO

  ENDIF ! IF IN STRATOSPHERE
  ENDDO !--end of loop klev
  ENDDO !--end of loop klon

END SUBROUTINE COAGULATE
