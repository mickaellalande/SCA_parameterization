!
! $Id: calfis.F 2604 2016-07-26 15:37:18Z emillour $
!
C
C
      SUBROUTINE calfis(lafin,
     $                  jD_cur, jH_cur,
     $                  pucov,
     $                  pvcov,
     $                  pteta,
     $                  pq,
     $                  pmasse,
     $                  pps,
     $                  pp,
     $                  ppk,
     $                  pphis,
     $                  pphi,
     $                  pducov,
     $                  pdvcov,
     $                  pdteta,
     $                  pdq,
     $                  flxw,
     $                  pdufi,
     $                  pdvfi,
     $                  pdhfi,
     $                  pdqfi,
     $                  pdpsfi)
c
c    Auteur :  P. Le Van, F. Hourdin 
c   .........
      USE infotrac, ONLY: nqtot, niadv, tname
      USE control_mod, ONLY: planet_type, nsplit_phys
#ifdef CPP_PHYS
      USE callphysiq_mod, ONLY: call_physiq
#endif 
      USE comconst_mod, ONLY: cpp, daysec, dtphys, dtvr, kappa, pi
      USE comvert_mod, ONLY: preff, presnivs
      
      IMPLICIT NONE
c=======================================================================
c
c   1. rearrangement des tableaux et transformation
c      variables dynamiques  >  variables physiques
c   2. calcul des termes physiques
c   3. retransformation des tendances physiques en tendances dynamiques
c
c   remarques:
c   ----------
c
c    - les vents sont donnes dans la physique par leurs composantes 
c      naturelles.
c    - la variable thermodynamique de la physique est une variable
c      intensive :   T 
c      pour la dynamique on prend    T * ( preff / p(l) ) **kappa
c    - les deux seules variables dependant de la geometrie necessaires
c      pour la physique sont la latitude pour le rayonnement et 
c      l'aire de la maille quand on veut integrer une grandeur 
c      horizontalement.
c    - les points de la physique sont les points scalaires de la 
c      la dynamique; numerotation:
c          1 pour le pole nord
c          (jjm-1)*iim pour l'interieur du domaine
c          ngridmx pour le pole sud
c      ---> ngridmx=2+(jjm-1)*iim
c
c     Input :
c     -------
c       pucov           covariant zonal velocity
c       pvcov           covariant meridional velocity 
c       pteta           potential temperature
c       pps             surface pressure
c       pmasse          masse d'air dans chaque maille
c       pts             surface temperature  (K)
c       callrad         clef d'appel au rayonnement
c
c    Output :
c    --------
c        pdufi          tendency for the natural zonal velocity (ms-1)
c        pdvfi          tendency for the natural meridional velocity 
c        pdhfi          tendency for the potential temperature
c        pdtsfi         tendency for the surface temperature
c
c        pdtrad         radiative tendencies  \  both input
c        pfluxrad       radiative fluxes      /  and output
c
c=======================================================================
c
c-----------------------------------------------------------------------
c
c    0.  Declarations :
c    ------------------

      include "dimensions.h"
      include "paramet.h"

      INTEGER ngridmx
      PARAMETER( ngridmx = 2+(jjm-1)*iim - 1/jjm   )

      include "comgeom2.h"
      include "iniprint.h"

c    Arguments :
c    -----------
      LOGICAL,INTENT(IN) ::  lafin ! .true. for the very last call to physics
      REAL,INTENT(IN):: jD_cur, jH_cur
      REAL,INTENT(IN) :: pvcov(iip1,jjm,llm) ! covariant meridional velocity
      REAL,INTENT(IN) :: pucov(iip1,jjp1,llm) ! covariant zonal velocity
      REAL,INTENT(IN) :: pteta(iip1,jjp1,llm) ! potential temperature
      REAL,INTENT(IN) :: pmasse(iip1,jjp1,llm) ! mass in each cell ! not used
      REAL,INTENT(IN) :: pq(iip1,jjp1,llm,nqtot) ! tracers
      REAL,INTENT(IN) :: pphis(iip1,jjp1) ! surface geopotential
      REAL,INTENT(IN) :: pphi(iip1,jjp1,llm) ! geopotential

      REAL,INTENT(IN) :: pdvcov(iip1,jjm,llm) ! dynamical tendency on vcov
      REAL,INTENT(IN) :: pducov(iip1,jjp1,llm) ! dynamical tendency on ucov
      REAL,INTENT(IN) :: pdteta(iip1,jjp1,llm) ! dynamical tendency on teta
      ! NB: pdteta is used only to compute pcvgt which is in fact not used...
      REAL,INTENT(IN) :: pdq(iip1,jjp1,llm,nqtot) ! dynamical tendency on tracers
      ! NB: pdq is only used to compute pcvgq which is in fact not used...

      REAL,INTENT(IN) :: pps(iip1,jjp1) ! surface pressure (Pa)
      REAL,INTENT(IN) :: pp(iip1,jjp1,llmp1) ! pressure at mesh interfaces (Pa)
      REAL,INTENT(IN) :: ppk(iip1,jjp1,llm) ! Exner at mid-layer
      REAL,INTENT(IN) :: flxw(iip1,jjp1,llm) ! Vertical mass flux on lower mesh interfaces (kg/s) (on llm because flxw(:,:,llm+1)=0)

      ! tendencies (in */s) from the physics
      REAL,INTENT(OUT) :: pdvfi(iip1,jjm,llm) ! tendency on covariant meridional wind
      REAL,INTENT(OUT) :: pdufi(iip1,jjp1,llm) ! tendency on covariant zonal wind
      REAL,INTENT(OUT) :: pdhfi(iip1,jjp1,llm) ! tendency on potential temperature (K/s)
      REAL,INTENT(OUT) :: pdqfi(iip1,jjp1,llm,nqtot) ! tendency on tracers
      REAL,INTENT(OUT) :: pdpsfi(iip1,jjp1) ! tendency on surface pressure (Pa/s)


c    Local variables :
c    -----------------

      INTEGER i,j,l,ig0,ig,iq,iiq
      REAL zpsrf(ngridmx)
      REAL zplev(ngridmx,llm+1),zplay(ngridmx,llm)
      REAL zphi(ngridmx,llm),zphis(ngridmx)
c
      REAL zrot(iip1,jjm,llm) ! AdlC May 2014
      REAL zufi(ngridmx,llm), zvfi(ngridmx,llm)
      REAL zrfi(ngridmx,llm) ! relative wind vorticity
      REAL ztfi(ngridmx,llm),zqfi(ngridmx,llm,nqtot)
      REAL zpk(ngridmx,llm)
c
      REAL pcvgu(ngridmx,llm), pcvgv(ngridmx,llm)
      REAL pcvgt(ngridmx,llm), pcvgq(ngridmx,llm,2)
c
      REAL zdufi(ngridmx,llm),zdvfi(ngridmx,llm)
      REAL zdtfi(ngridmx,llm),zdqfi(ngridmx,llm,nqtot)
      REAL zdpsrf(ngridmx)
c
      REAL zdufic(ngridmx,llm),zdvfic(ngridmx,llm)
      REAL zdtfic(ngridmx,llm),zdqfic(ngridmx,llm,nqtot)
      REAL jH_cur_split,zdt_split
      LOGICAL debut_split,lafin_split
      INTEGER isplit

      REAL zsin(iim),zcos(iim),z1(iim)
      REAL zsinbis(iim),zcosbis(iim),z1bis(iim)
      REAL unskap, pksurcp
c
      REAL flxwfi(ngridmx,llm)  ! Flux de masse verticale sur la grille physiq
c
      
      REAL SSUM

      LOGICAL,SAVE :: firstcal=.true., debut=.true.
!      REAL rdayvrai

c
c-----------------------------------------------------------------------
c
c    1. Initialisations :
c    --------------------
c
c
      IF ( firstcal )  THEN
        debut = .TRUE.
        IF (ngridmx.NE.2+(jjm-1)*iim) THEN
         write(lunout,*) 'STOP dans calfis'
         write(lunout,*)
     &   'La dimension ngridmx doit etre egale a 2 + (jjm-1)*iim'
         write(lunout,*) '  ngridmx  jjm   iim   '
         write(lunout,*) ngridmx,jjm,iim
         STOP
        ENDIF
      ELSE
        debut = .FALSE.
      ENDIF ! of IF (firstcal)

c
c
c-----------------------------------------------------------------------
c   40. transformation des variables dynamiques en variables physiques:
c   ---------------------------------------------------------------

c   41. pressions au sol (en Pascals)
c   ----------------------------------

       
      zpsrf(1) = pps(1,1)

      ig0  = 2
      DO j = 2,jjm
         CALL SCOPY( iim,pps(1,j),1,zpsrf(ig0), 1 )
         ig0 = ig0+iim
      ENDDO

      zpsrf(ngridmx) = pps(1,jjp1)


c   42. pression intercouches et fonction d'Exner:
c
c   -----------------------------------------------------------------
c     .... zplev  definis aux (llm +1) interfaces des couches  ....
c     .... zplay  definis aux (  llm )    milieux des couches  .... 
c   -----------------------------------------------------------------

c    ...    Exner = cp * ( p(l) / preff ) ** kappa     ....
c
       unskap   = 1./ kappa
c
      DO l = 1, llm
        zpk(   1,l ) = ppk(1,1,l)
        zplev( 1,l ) = pp(1,1,l)
        ig0 = 2
          DO j = 2, jjm
             DO i =1, iim
              zpk(   ig0,l ) = ppk(i,j,l)
              zplev( ig0,l ) = pp(i,j,l)
              ig0 = ig0 +1
             ENDDO
          ENDDO
        zpk(   ngridmx,l ) = ppk(1,jjp1,l)
        zplev( ngridmx,l ) = pp(1,jjp1,l)
      ENDDO
        zplev( 1,llmp1 ) = pp(1,1,llmp1)
        ig0 = 2
          DO j = 2, jjm
             DO i =1, iim
              zplev( ig0,llmp1 ) = pp(i,j,llmp1)
              ig0 = ig0 +1
             ENDDO
          ENDDO
        zplev( ngridmx,llmp1 ) = pp(1,jjp1,llmp1)
c
c

c   43. temperature naturelle (en K) et pressions milieux couches .
c   ---------------------------------------------------------------

      DO l=1,llm

         pksurcp     =  ppk(1,1,l) / cpp
         zplay(1,l)  =  preff * pksurcp ** unskap
         ztfi(1,l)   =  pteta(1,1,l) *  pksurcp
         pcvgt(1,l)  =  pdteta(1,1,l) * pksurcp / pmasse(1,1,l)
         ig0         = 2

         DO j = 2, jjm
            DO i = 1, iim
              pksurcp        = ppk(i,j,l) / cpp
              zplay(ig0,l)   = preff * pksurcp ** unskap
              ztfi(ig0,l)    = pteta(i,j,l)  * pksurcp
              pcvgt(ig0,l)   = pdteta(i,j,l) * pksurcp / pmasse(i,j,l)
              ig0            = ig0 + 1
            ENDDO
         ENDDO

         pksurcp       = ppk(1,jjp1,l) / cpp
         zplay(ig0,l)  = preff * pksurcp ** unskap
         ztfi (ig0,l)  = pteta(1,jjp1,l)  * pksurcp
         pcvgt(ig0,l)  = pdteta(1,jjp1,l) * pksurcp/ pmasse(1,jjp1,l)

      ENDDO

c   43.bis traceurs
c   ---------------
c
      DO iq=1,nqtot
          iiq=niadv(iq) 
         DO l=1,llm
            zqfi(1,l,iq) = pq(1,1,l,iiq)
            ig0          = 2
            DO j=2,jjm
               DO i = 1, iim
                  zqfi(ig0,l,iq)  = pq(i,j,l,iiq)
                  ig0             = ig0 + 1
               ENDDO
            ENDDO
            zqfi(ig0,l,iq) = pq(1,jjp1,l,iiq)
         ENDDO
      ENDDO

c   convergence dynamique pour les traceurs "EAU"
! Earth-specific treatment of first 2 tracers (water)
       if (planet_type=="earth") then
        DO iq=1,2
         DO l=1,llm
            pcvgq(1,l,iq)= pdq(1,1,l,iq) / pmasse(1,1,l)
            ig0          = 2
            DO j=2,jjm
               DO i = 1, iim
                  pcvgq(ig0,l,iq) = pdq(i,j,l,iq) / pmasse(i,j,l)
                  ig0             = ig0 + 1
               ENDDO
            ENDDO
            pcvgq(ig0,l,iq)= pdq(1,jjp1,l,iq) / pmasse(1,jjp1,l)
         ENDDO
        ENDDO
       endif ! of if (planet_type=="earth")


c   Geopotentiel calcule par rapport a la surface locale:
c   -----------------------------------------------------

      CALL gr_dyn_fi(llm,iip1,jjp1,ngridmx,pphi,zphi)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,pphis,zphis)
      DO l=1,llm
         DO ig=1,ngridmx
           zphi(ig,l)=zphi(ig,l)-zphis(ig)
         ENDDO
      ENDDO

c   ....  Calcul de la vitesse  verticale  ( en Pa*m*s  ou Kg/s )  ....
c JG : ancien calcule de omega utilise dans physiq.F. Maintenant le flux 
c    de masse est calclue dans advtrac.F  
c      DO l=1,llm
c        pvervel(1,l)=pw(1,1,l) * g /apoln
c        ig0=2
c       DO j=2,jjm
c           DO i = 1, iim
c              pvervel(ig0,l) = pw(i,j,l) * g * unsaire(i,j)
c              ig0 = ig0 + 1
c           ENDDO
c       ENDDO
c        pvervel(ig0,l)=pw(1,jjp1,l) * g /apols
c      ENDDO

c
c   45. champ u:
c   ------------

      DO 50 l=1,llm

         DO 25 j=2,jjm
            ig0 = 1+(j-2)*iim
            zufi(ig0+1,l)= 0.5 * 
     $      ( pucov(iim,j,l)/cu(iim,j) + pucov(1,j,l)/cu(1,j) )
            pcvgu(ig0+1,l)= 0.5 * 
     $      ( pducov(iim,j,l)/cu(iim,j) + pducov(1,j,l)/cu(1,j) )
            DO 10 i=2,iim
               zufi(ig0+i,l)= 0.5 *
     $         ( pucov(i-1,j,l)/cu(i-1,j) + pucov(i,j,l)/cu(i,j) )
               pcvgu(ig0+i,l)= 0.5 *
     $         ( pducov(i-1,j,l)/cu(i-1,j) + pducov(i,j,l)/cu(i,j) )
10         CONTINUE
25      CONTINUE

50    CONTINUE


C  Alvaro de la Camara (May 2014)
C  46.1 Calcul de la vorticite et passage sur la grille physique
C  --------------------------------------------------------------
      DO l=1,llm
        do i=1,iim
          do j=1,jjm
            zrot(i,j,l) = (pvcov(i+1,j,l) - pvcov(i,j,l)
     $                   + pucov(i,j+1,l) - pucov(i,j,l)) 
     $                   / (cu(i,j)+cu(i,j+1)) 
     $                   / (cv(i+1,j)+cv(i,j)) *4
          enddo
        enddo
      ENDDO

c   46.champ v:
c   -----------

      DO l=1,llm
         DO j=2,jjm
            ig0=1+(j-2)*iim
            DO i=1,iim
               zvfi(ig0+i,l)= 0.5 *
     $         ( pvcov(i,j-1,l)/cv(i,j-1) + pvcov(i,j,l)/cv(i,j) )
               pcvgv(ig0+i,l)= 0.5 *
     $         ( pdvcov(i,j-1,l)/cv(i,j-1) + pdvcov(i,j,l)/cv(i,j) )
            ENDDO
               zrfi(ig0 + 1,l)= 0.25 *(zrot(iim,j-1,l)+zrot(iim,j,l)
     &                                +zrot(1,j-1,l)+zrot(1,j,l))
            DO i=2,iim
               zrfi(ig0 + i,l)= 0.25 *(zrot(i-1,j-1,l)+zrot(i-1,j,l)
     $                   +zrot(i,j-1,l)+zrot(i,j,l))   !  AdlC MAY 2014
            ENDDO
         ENDDO
      ENDDO


c   47. champs de vents aux pole nord   
c   ------------------------------
c        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
c        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

      DO l=1,llm

         z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,1,l)/cv(1,1)
         z1bis(1)=(rlonu(1)-rlonu(iim)+2.*pi)*pdvcov(1,1,l)/cv(1,1)
         DO i=2,iim
            z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,1,l)/cv(i,1)
            z1bis(i)=(rlonu(i)-rlonu(i-1))*pdvcov(i,1,l)/cv(i,1)
         ENDDO

         DO i=1,iim
            zcos(i)   = COS(rlonv(i))*z1(i)
            zcosbis(i)= COS(rlonv(i))*z1bis(i)
            zsin(i)   = SIN(rlonv(i))*z1(i)
            zsinbis(i)= SIN(rlonv(i))*z1bis(i)
         ENDDO

         zufi(1,l)  = SSUM(iim,zcos,1)/pi
         pcvgu(1,l) = SSUM(iim,zcosbis,1)/pi
         zvfi(1,l)  = SSUM(iim,zsin,1)/pi
         pcvgv(1,l) = SSUM(iim,zsinbis,1)/pi
         zrfi(1, l) = 0.
      ENDDO


c   48. champs de vents aux pole sud:
c   ---------------------------------
c        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
c        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

      DO l=1,llm

         z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,jjm,l)/cv(1,jjm)
         z1bis(1)=(rlonu(1)-rlonu(iim)+2.*pi)*pdvcov(1,jjm,l)/cv(1,jjm)
         DO i=2,iim
            z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,jjm,l)/cv(i,jjm)
            z1bis(i)=(rlonu(i)-rlonu(i-1))*pdvcov(i,jjm,l)/cv(i,jjm)
         ENDDO

         DO i=1,iim
            zcos(i)    = COS(rlonv(i))*z1(i)
            zcosbis(i) = COS(rlonv(i))*z1bis(i)
            zsin(i)    = SIN(rlonv(i))*z1(i)
            zsinbis(i) = SIN(rlonv(i))*z1bis(i)
         ENDDO

         zufi(ngridmx,l)  = SSUM(iim,zcos,1)/pi
         pcvgu(ngridmx,l) = SSUM(iim,zcosbis,1)/pi
         zvfi(ngridmx,l)  = SSUM(iim,zsin,1)/pi
         pcvgv(ngridmx,l) = SSUM(iim,zsinbis,1)/pi
         zrfi(ngridmx, l) = 0.
      ENDDO
c
c On change de grille, dynamique vers physiq, pour le flux de masse verticale
      CALL gr_dyn_fi(llm,iip1,jjp1,ngridmx,flxw,flxwfi)

c-----------------------------------------------------------------------
c   Appel de la physique:
c   ---------------------



!      write(lunout,*) 'PHYSIQUE AVEC NSPLIT_PHYS=',nsplit_phys
      zdt_split=dtphys/nsplit_phys
      zdufic(:,:)=0.
      zdvfic(:,:)=0.
      zdtfic(:,:)=0.
      zdqfic(:,:,:)=0.

#ifdef CPP_PHYS

       do isplit=1,nsplit_phys

         jH_cur_split=jH_cur+(isplit-1) * dtvr / (daysec *nsplit_phys)
         debut_split=debut.and.isplit==1
         lafin_split=lafin.and.isplit==nsplit_phys

        CALL call_physiq(ngridmx,llm,nqtot,tname,
     &                   debut_split,lafin_split,
     &                   jD_cur,jH_cur_split,zdt_split,
     &                   zplev,zplay,
     &                   zpk,zphi,zphis,
     &                   presnivs,
     &                   zufi,zvfi,zrfi,ztfi,zqfi,
     &                   flxwfi,pducov,
     &                   zdufi,zdvfi,zdtfi,zdqfi,zdpsrf)
                             
!      if (planet_type=="earth") then
!
!         CALL physiq (ngridmx,
!     .             llm,
!     .             debut_split,
!     .             lafin_split,
!     .             jD_cur,
!     .             jH_cur_split,
!     .             zdt_split,
!     .             zplev,
!     .             zplay,
!     .             zphi,
!     .             zphis,
!     .             presnivs,
!     .             zufi,
!     .             zvfi, zrfi,
!     .             ztfi,
!     .             zqfi,
!     .             flxwfi,
!     .             zdufi,
!     .             zdvfi,
!     .             zdtfi,
!     .             zdqfi,
!     .             zdpsrf,
!     .             pducov)
!
!      else if ( planet_type=="generic" ) then
!
!         CALL physiq (ngridmx,     !! ngrid
!     .             llm,            !! nlayer
!     .             nqtot,          !! nq
!     .             tname,          !! tracer names from dynamical core (given in infotrac)
!     .             debut_split,    !! firstcall 
!     .             lafin_split,    !! lastcall
!     .             jD_cur,         !! pday. see leapfrog
!     .             jH_cur_split,   !! ptime "fraction of day"
!     .             zdt_split,      !! ptimestep
!     .             zplev,          !! pplev
!     .             zplay,          !! pplay
!     .             zphi,           !! pphi
!     .             zufi,           !! pu
!     .             zvfi,           !! pv
!     .             ztfi,           !! pt
!     .             zqfi,           !! pq
!     .             flxwfi,         !! pw !! or 0. anyway this is for diagnostic. not used in physiq.
!     .             zdufi,          !! pdu
!     .             zdvfi,          !! pdv
!     .             zdtfi,          !! pdt
!     .             zdqfi,          !! pdq
!     .             zdpsrf,         !! pdpsrf
!     .             tracerdyn)      !! tracerdyn <-- utilite ???
!
!      endif ! of if (planet_type=="earth")

         zufi(:,:)=zufi(:,:)+zdufi(:,:)*zdt_split
         zvfi(:,:)=zvfi(:,:)+zdvfi(:,:)*zdt_split
         ztfi(:,:)=ztfi(:,:)+zdtfi(:,:)*zdt_split
         zqfi(:,:,:)=zqfi(:,:,:)+zdqfi(:,:,:)*zdt_split

         zdufic(:,:)=zdufic(:,:)+zdufi(:,:)
         zdvfic(:,:)=zdvfic(:,:)+zdvfi(:,:)
         zdtfic(:,:)=zdtfic(:,:)+zdtfi(:,:)
         zdqfic(:,:,:)=zdqfic(:,:,:)+zdqfi(:,:,:)

       enddo ! of do isplit=1,nsplit_phys

#endif
! of #ifdef CPP_PHYS

      zdufi(:,:)=zdufic(:,:)/nsplit_phys
      zdvfi(:,:)=zdvfic(:,:)/nsplit_phys
      zdtfi(:,:)=zdtfic(:,:)/nsplit_phys
      zdqfi(:,:,:)=zdqfic(:,:,:)/nsplit_phys


500   CONTINUE

c-----------------------------------------------------------------------
c   transformation des tendances physiques en tendances dynamiques:
c   ---------------------------------------------------------------

c  tendance sur la pression :
c  -----------------------------------

      CALL gr_fi_dyn(1,ngridmx,iip1,jjp1,zdpsrf,pdpsfi)
c
c   62. enthalpie potentielle
c   ---------------------

      DO l=1,llm

         DO i=1,iip1
          pdhfi(i,1,l)    = cpp *  zdtfi(1,l)      / ppk(i, 1  ,l)
          pdhfi(i,jjp1,l) = cpp *  zdtfi(ngridmx,l)/ ppk(i,jjp1,l)
         ENDDO

         DO j=2,jjm
            ig0=1+(j-2)*iim
            DO i=1,iim
               pdhfi(i,j,l) = cpp * zdtfi(ig0+i,l) / ppk(i,j,l)
            ENDDO
               pdhfi(iip1,j,l) =  pdhfi(1,j,l)
         ENDDO

      ENDDO


c   62. humidite specifique
c   ---------------------
! Ehouarn: removed this useless bit: was overwritten at step 63 anyways
!      DO iq=1,nqtot
!         DO l=1,llm
!            DO i=1,iip1
!               pdqfi(i,1,l,iq)    = zdqfi(1,l,iq)
!               pdqfi(i,jjp1,l,iq) = zdqfi(ngridmx,l,iq)
!            ENDDO
!            DO j=2,jjm
!               ig0=1+(j-2)*iim
!               DO i=1,iim
!                  pdqfi(i,j,l,iq) = zdqfi(ig0+i,l,iq)
!               ENDDO
!               pdqfi(iip1,j,l,iq) = pdqfi(1,j,l,iq)
!            ENDDO
!         ENDDO
!      ENDDO

c   63. traceurs
c   ------------
C     initialisation des tendances
      pdqfi(:,:,:,:)=0.
C
      DO iq=1,nqtot
         iiq=niadv(iq)
         DO l=1,llm
            DO i=1,iip1
               pdqfi(i,1,l,iiq)    = zdqfi(1,l,iq)
               pdqfi(i,jjp1,l,iiq) = zdqfi(ngridmx,l,iq)
            ENDDO
            DO j=2,jjm
               ig0=1+(j-2)*iim
               DO i=1,iim
                  pdqfi(i,j,l,iiq) = zdqfi(ig0+i,l,iq)
               ENDDO
               pdqfi(iip1,j,l,iiq) = pdqfi(1,j,l,iq)
            ENDDO
         ENDDO
      ENDDO

c   65. champ u:
c   ------------

      DO l=1,llm

         DO i=1,iip1
            pdufi(i,1,l)    = 0.
            pdufi(i,jjp1,l) = 0.
         ENDDO

         DO j=2,jjm
            ig0=1+(j-2)*iim
            DO i=1,iim-1
               pdufi(i,j,l)=
     $         0.5*(zdufi(ig0+i,l)+zdufi(ig0+i+1,l))*cu(i,j)
            ENDDO
            pdufi(iim,j,l)=
     $      0.5*(zdufi(ig0+1,l)+zdufi(ig0+iim,l))*cu(iim,j)
            pdufi(iip1,j,l)=pdufi(1,j,l)
         ENDDO

      ENDDO


c   67. champ v:
c   ------------

      DO l=1,llm

         DO j=2,jjm-1
            ig0=1+(j-2)*iim
            DO i=1,iim
               pdvfi(i,j,l)=
     $         0.5*(zdvfi(ig0+i,l)+zdvfi(ig0+i+iim,l))*cv(i,j)
            ENDDO
            pdvfi(iip1,j,l) = pdvfi(1,j,l)
         ENDDO
      ENDDO


c   68. champ v pres des poles:
c   ---------------------------
c      v = U * cos(long) + V * SIN(long)

      DO l=1,llm

         DO i=1,iim
            pdvfi(i,1,l)=
     $      zdufi(1,l)*COS(rlonv(i))+zdvfi(1,l)*SIN(rlonv(i))
            pdvfi(i,jjm,l)=zdufi(ngridmx,l)*COS(rlonv(i))
     $      +zdvfi(ngridmx,l)*SIN(rlonv(i))
            pdvfi(i,1,l)=
     $      0.5*(pdvfi(i,1,l)+zdvfi(i+1,l))*cv(i,1)
            pdvfi(i,jjm,l)=
     $      0.5*(pdvfi(i,jjm,l)+zdvfi(ngridmx-iip1+i,l))*cv(i,jjm)
          ENDDO

         pdvfi(iip1,1,l)  = pdvfi(1,1,l)
         pdvfi(iip1,jjm,l)= pdvfi(1,jjm,l)

      ENDDO

c-----------------------------------------------------------------------

700   CONTINUE
 
      firstcal = .FALSE.

      RETURN
      END
