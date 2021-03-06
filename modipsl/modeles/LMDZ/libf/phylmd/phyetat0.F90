! $Id: phyetat0.F90 3328 2018-05-16 16:10:44Z musat $

SUBROUTINE phyetat0 (fichnom, clesphy0, tabcntr0)

  USE dimphy, only: klon, zmasq, klev
  USE iophy, ONLY : init_iophy_new
  USE ocean_cpl_mod,    ONLY : ocean_cpl_init
  USE fonte_neige_mod,  ONLY : fonte_neige_init
  USE pbl_surface_mod,  ONLY : pbl_surface_init
  USE surface_data,     ONLY : type_ocean, version_ocean
  USE phys_state_var_mod, ONLY : ancien_ok, clwcon, detr_therm, dtime, &
       qsol, fevap, z0m, z0h, agesno, &
       du_gwd_rando, du_gwd_front, entr_therm, f0, fm_therm, &
       falb_dir, falb_dif, prw_ancien, prlw_ancien, prsw_ancien, &
       ftsol, pbl_tke, pctsrf, q_ancien, ql_ancien, qs_ancien, radpas, radsol, rain_fall, ratqs, &
       rnebcon, rugoro, sig1, snow_fall, solaire_etat0, sollw, sollwdown, &
       solsw, t_ancien, u_ancien, v_ancien, w01, wake_cstar, wake_deltaq, &
       wake_deltat, wake_delta_pbl_TKE, delta_tsurf, wake_fip, wake_pe, &
       wake_s, wake_dens, zgam, zmax0, zmea, zpic, zsig, &
       zstd, zthe, zval, ale_bl, ale_bl_trig, alp_bl, u10m, v10m, treedrg, &
       ale_wake, ale_bl_stat
!FC
  USE geometry_mod, ONLY : longitude_deg, latitude_deg
  USE iostart, ONLY : close_startphy, get_field, get_var, open_startphy
  USE infotrac_phy, only: nbtr, nqo, type_trac, tname, niadv
  USE traclmdz_mod,    ONLY : traclmdz_from_restart
  USE carbon_cycle_mod, ONLY : carbon_cycle_tr, carbon_cycle_cpl, co2_send
  USE indice_sol_mod, only: nbsrf, is_ter, epsfra, is_lic, is_oce, is_sic
  USE ocean_slab_mod, ONLY: nslay, tslab, seaice, tice, ocean_slab_init
  USE time_phylmdz_mod, ONLY: init_iteration, pdtphys, itau_phy

  IMPLICIT none
  !======================================================================
  ! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: Lecture de l'etat initial pour la physique
  !======================================================================
  include "netcdf.inc"
  include "dimsoil.h"
  include "clesphys.h"
  include "thermcell.h"
  include "compbl.h"
  include "YOMCST.h"
  !======================================================================
  CHARACTER*(*) fichnom

  ! les variables globales lues dans le fichier restart

  REAL tsoil(klon, nsoilmx, nbsrf)
  REAL qsurf(klon, nbsrf)
  REAL snow(klon, nbsrf)
  real fder(klon)
  REAL run_off_lic_0(klon)
  REAL fractint(klon)
  REAL trs(klon, nbtr)
  REAL zts(klon)
  ! pour drag arbres FC
  REAL drg_ter(klon,klev)

  CHARACTER*6 ocean_in
  LOGICAL ok_veget_in

  INTEGER        longcles
  PARAMETER    ( longcles = 20 )
  REAL clesphy0( longcles )

  REAL xmin, xmax

  INTEGER nid, nvarid
  INTEGER ierr, i, nsrf, isoil , k
  INTEGER length
  PARAMETER (length=100)
  INTEGER it, iiq, isw
  REAL tab_cntrl(length), tabcntr0(length)
  CHARACTER*7 str7
  CHARACTER*2 str2
  LOGICAL :: found,phyetat0_get,phyetat0_srf
  REAL :: lon_startphy(klon), lat_startphy(klon)

  ! FH1D
  !     real iolat(jjm+1)
  !real iolat(jjm+1-1/(iim*jjm))

  ! Ouvrir le fichier contenant l'etat initial:

  CALL open_startphy(fichnom)

  ! Lecture des parametres de controle:

  CALL get_var("controle", tab_cntrl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
  ! Les constantes de la physiques sont lues dans la physique seulement.
  ! Les egalites du type
  !             tab_cntrl( 5 )=clesphy0(1)
  ! sont remplacees par
  !             clesphy0(1)=tab_cntrl( 5 )
  ! On inverse aussi la logique.
  ! On remplit les tab_cntrl avec les parametres lus dans les .def
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO i = 1, length
     tabcntr0( i ) = tab_cntrl( i )
  ENDDO

  tab_cntrl(1)=pdtphys
  tab_cntrl(2)=radpas

  ! co2_ppm : value from the previous time step
  IF (carbon_cycle_tr .OR. carbon_cycle_cpl) THEN
     co2_ppm = tab_cntrl(3)
     RCO2    = co2_ppm * 1.0e-06  * 44.011/28.97 
     ! ELSE : keep value from .def
  END IF

  ! co2_ppm0 : initial value of atmospheric CO2 (from create_etat0_limit.e .def)
  co2_ppm0   = tab_cntrl(16)

  solaire_etat0      = tab_cntrl(4)
  tab_cntrl(5)=iflag_con
  tab_cntrl(6)=nbapp_rad

  if (iflag_cycle_diurne.GE.1) tab_cntrl( 7) = iflag_cycle_diurne
  if (soil_model) tab_cntrl( 8) =1.
  if (new_oliq) tab_cntrl( 9) =1.
  if (ok_orodr) tab_cntrl(10) =1.
  if (ok_orolf) tab_cntrl(11) =1.
  if (ok_limitvrai) tab_cntrl(12) =1.

  itau_phy = tab_cntrl(15)

  clesphy0(1)=tab_cntrl( 5 )
  clesphy0(2)=tab_cntrl( 6 )
  clesphy0(3)=tab_cntrl( 7 )
  clesphy0(4)=tab_cntrl( 8 )
  clesphy0(5)=tab_cntrl( 9 )
  clesphy0(6)=tab_cntrl( 10 )
  clesphy0(7)=tab_cntrl( 11 )
  clesphy0(8)=tab_cntrl( 12 )

  ! set time iteration
   CALL init_iteration(itau_phy)

  ! read latitudes and make a sanity check (because already known from dyn)
  CALL get_field("latitude",lat_startphy)
  DO i=1,klon
    IF (ABS(lat_startphy(i)-latitude_deg(i))>=1) THEN
      WRITE(*,*) "phyetat0: Error! Latitude discrepancy wrt startphy file:",&
                 " i=",i," lat_startphy(i)=",lat_startphy(i),&
                 " latitude_deg(i)=",latitude_deg(i)
      ! This is presumably serious enough to abort run
      CALL abort_physic("phyetat0","discrepancy in latitudes!",1)
    ENDIF
    IF (ABS(lat_startphy(i)-latitude_deg(i))>=0.0001) THEN
      WRITE(*,*) "phyetat0: Warning! Latitude discrepancy wrt startphy file:",&
                 " i=",i," lat_startphy(i)=",lat_startphy(i),&
                 " latitude_deg(i)=",latitude_deg(i)
    ENDIF
  ENDDO

  ! read longitudes and make a sanity check (because already known from dyn)
  CALL get_field("longitude",lon_startphy)
  DO i=1,klon
    IF (ABS(lon_startphy(i)-longitude_deg(i))>=1) THEN
      WRITE(*,*) "phyetat0: Error! Longitude discrepancy wrt startphy file:",&
                 " i=",i," lon_startphy(i)=",lon_startphy(i),&
                 " longitude_deg(i)=",longitude_deg(i)
      ! This is presumably serious enough to abort run
      CALL abort_physic("phyetat0","discrepancy in longitudes!",1)
    ENDIF
    IF (ABS(lon_startphy(i)-longitude_deg(i))>=0.0001) THEN
      WRITE(*,*) "phyetat0: Warning! Longitude discrepancy wrt startphy file:",&
                 " i=",i," lon_startphy(i)=",lon_startphy(i),&
                 " longitude_deg(i)=",longitude_deg(i)
    ENDIF
  ENDDO

  ! Lecture du masque terre mer

  CALL get_field("masque", zmasq, found)
  IF (.NOT. found) THEN
     PRINT*, 'phyetat0: Le champ <masque> est absent'
     PRINT *, 'fichier startphy non compatible avec phyetat0'
  ENDIF

  ! Lecture des fractions pour chaque sous-surface

  ! initialisation des sous-surfaces

  pctsrf = 0.

  ! fraction de terre

  CALL get_field("FTER", pctsrf(:, is_ter), found)
  IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FTER> est absent'

  ! fraction de glace de terre

  CALL get_field("FLIC", pctsrf(:, is_lic), found)
  IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FLIC> est absent'

  ! fraction d'ocean

  CALL get_field("FOCE", pctsrf(:, is_oce), found)
  IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FOCE> est absent'

  ! fraction glace de mer

  CALL get_field("FSIC", pctsrf(:, is_sic), found)
  IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FSIC> est absent'

  !  Verification de l'adequation entre le masque et les sous-surfaces

  fractint( 1 : klon) = pctsrf(1 : klon, is_ter)  &
       + pctsrf(1 : klon, is_lic)
  DO i = 1 , klon
     IF ( abs(fractint(i) - zmasq(i) ) .GT. EPSFRA ) THEN
        WRITE(*, *) 'phyetat0: attention fraction terre pas ',  &
             'coherente ', i, zmasq(i), pctsrf(i, is_ter) &
             , pctsrf(i, is_lic)
        WRITE(*, *) 'Je force la coherence zmasq=fractint'
        zmasq(i) = fractint(i)
     ENDIF
  END DO
  fractint (1 : klon) =  pctsrf(1 : klon, is_oce)  &
       + pctsrf(1 : klon, is_sic)
  DO i = 1 , klon
     IF ( abs( fractint(i) - (1. - zmasq(i))) .GT. EPSFRA ) THEN
        WRITE(*, *) 'phyetat0 attention fraction ocean pas ',  &
             'coherente ', i, zmasq(i) , pctsrf(i, is_oce) &
             , pctsrf(i, is_sic)
        WRITE(*, *) 'Je force la coherence zmasq=1.-fractint'
        zmasq(i) = 1. - fractint(i)
     ENDIF
  END DO

!===================================================================
! Lecture des temperatures du sol:
!===================================================================

  found=phyetat0_get(1,ftsol(:,1),"TS","Surface temperature",283.)
  IF (found) THEN
     DO nsrf=2,nbsrf
        ftsol(:,nsrf)=ftsol(:,1)
     ENDDO
  ELSE
     found=phyetat0_srf(1,ftsol,"TS","Surface temperature",283.)
  ENDIF

!===================================================================
  ! Lecture des albedo difus et direct
!===================================================================

  DO nsrf = 1, nbsrf
     DO isw=1, nsw
        IF (isw.GT.99) THEN
           PRINT*, "Trop de bandes SW"
           call abort_physic("phyetat0", "", 1)
        ENDIF
        WRITE(str2, '(i2.2)') isw
        found=phyetat0_srf(1,falb_dir(:, isw,:),"A_dir_SW"//str2//"srf","Direct Albedo",0.2)
        found=phyetat0_srf(1,falb_dif(:, isw,:),"A_dif_SW"//str2//"srf","Direct Albedo",0.2)
     ENDDO
  ENDDO

  found=phyetat0_srf(1,u10m,"U10M","u a 10m",0.)
  found=phyetat0_srf(1,v10m,"V10M","v a 10m",0.)

!===================================================================
  ! Lecture des temperatures du sol profond:
!===================================================================

   DO isoil=1, nsoilmx
        IF (isoil.GT.99) THEN
           PRINT*, "Trop de couches "
           call abort_physic("phyetat0", "", 1)
        ENDIF
        WRITE(str2,'(i2.2)') isoil
        found=phyetat0_srf(1,tsoil(:, isoil,:),"Tsoil"//str2//"srf","Temp soil",0.)
        IF (.NOT. found) THEN
           PRINT*, "phyetat0: Le champ <Tsoil"//str7//"> est absent"
           PRINT*, "          Il prend donc la valeur de surface"
           tsoil(:, isoil, :)=ftsol(:, :)
        ENDIF
   ENDDO

!=======================================================================
! Lecture precipitation/evaporation
!=======================================================================

  found=phyetat0_srf(1,qsurf,"QS","Near surface hmidity",0.)
  found=phyetat0_get(1,qsol,"QSOL","Surface hmidity / bucket",0.)
  found=phyetat0_srf(1,snow,"SNOW","Surface snow",0.)
  found=phyetat0_srf(1,fevap,"EVAP","evaporation",0.)
  found=phyetat0_get(1,snow_fall,"snow_f","snow fall",0.)
  found=phyetat0_get(1,rain_fall,"rain_f","rain fall",0.)

!=======================================================================
! Radiation
!=======================================================================

  found=phyetat0_get(1,solsw,"solsw","net SW radiation surf",0.)
  found=phyetat0_get(1,sollw,"sollw","net LW radiation surf",0.)
  found=phyetat0_get(1,sollwdown,"sollwdown","down LW radiation surf",0.)
  IF (.NOT. found) THEN
     sollwdown = 0. ;  zts=0.
     do nsrf=1,nbsrf
        zts(:)=zts(:)+ftsol(:,nsrf)*pctsrf(:,nsrf)
     enddo
     sollwdown(:)=sollw(:)+RSIGMA*zts(:)**4
  ENDIF

  found=phyetat0_get(1,radsol,"RADS","Solar radiation",0.)
  found=phyetat0_get(1,fder,"fder","Flux derivative",0.) 


  ! Lecture de la longueur de rugosite 
  found=phyetat0_srf(1,z0m,"RUG","Z0m ancien",0.001)
  IF (found) THEN
     z0h(:,1:nbsrf)=z0m(:,1:nbsrf)
  ELSE
     found=phyetat0_srf(1,z0m,"Z0m","Roughness length, momentum ",0.001)
     found=phyetat0_srf(1,z0h,"Z0h","Roughness length, enthalpy ",0.001)
  ENDIF
!FC
  IF (ifl_pbltree>0) THEN
!CALL get_field("FTER", pctsrf(:, is_ter), found)
    treedrg(:,1:klev,1:nbsrf)= 0.0
    CALL get_field("treedrg_ter", drg_ter(:,:), found)
!  found=phyetat0_srf(1,treedrg,"treedrg","drag from vegetation" , 0.)
    !lecture du profile de freinage des arbres
    IF (.not. found ) THEN
      treedrg(:,1:klev,1:nbsrf)= 0.0
    ELSE
      treedrg(:,1:klev,is_ter)= drg_ter(:,:)
!     found=phyetat0_srf(klev,treedrg,"treedrg","freinage arbres",0.)
    ENDIF
  ELSE
    ! initialize treedrg(), because it will be written in restartphy.nc
    treedrg(:,:,:) = 0.0 
  ENDIF

  ! Lecture de l'age de la neige:
  found=phyetat0_srf(1,agesno,"AGESNO","SNOW AGE",0.001)

  ancien_ok=.true.
  ancien_ok=ancien_ok.AND.phyetat0_get(klev,t_ancien,"TANCIEN","TANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(klev,q_ancien,"QANCIEN","QANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(klev,ql_ancien,"QLANCIEN","QLANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(klev,qs_ancien,"QSANCIEN","QSANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(klev,u_ancien,"UANCIEN","UANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(klev,v_ancien,"VANCIEN","VANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(1,prw_ancien,"PRWANCIEN","PRWANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(1,prlw_ancien,"PRLWANCIEN","PRLWANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(1,prsw_ancien,"PRSWANCIEN","PRSWANCIEN",0.)

  ! Ehouarn: addtional tests to check if t_ancien, q_ancien contain
  !          dummy values (as is the case when generated by ce0l,
  !          or by iniaqua)
  if ( (maxval(q_ancien).eq.minval(q_ancien))       .or. &
       (maxval(ql_ancien).eq.minval(ql_ancien))     .or. &
       (maxval(qs_ancien).eq.minval(qs_ancien))     .or. &
       (maxval(prw_ancien).eq.minval(prw_ancien))   .or. &
       (maxval(prlw_ancien).eq.minval(prlw_ancien)) .or. &
       (maxval(prsw_ancien).eq.minval(prsw_ancien)) .or. &
       (maxval(t_ancien).eq.minval(t_ancien)) ) then
    ancien_ok=.false.
  endif

  found=phyetat0_get(klev,clwcon,"CLWCON","CLWCON",0.)
  found=phyetat0_get(klev,rnebcon,"RNEBCON","RNEBCON",0.)
  found=phyetat0_get(klev,ratqs,"RATQS","RATQS",0.)

  found=phyetat0_get(1,run_off_lic_0,"RUNOFFLIC0","RUNOFFLIC0",0.)

!==================================
!  TKE
!==================================
!
  IF (iflag_pbl>1) then
     found=phyetat0_srf(klev+1,pbl_tke,"TKE","Turb. Kinetic. Energ. ",1.e-8)
  ENDIF

  IF (iflag_pbl>1 .AND. iflag_wake>=1  .AND. iflag_pbl_split >=1 ) then
    found=phyetat0_srf(klev+1,wake_delta_pbl_tke,"DELTATKE","Del TKE wk/env",0.)
    found=phyetat0_srf(1,delta_tsurf,"DELTA_TSURF","Delta Ts wk/env ",0.)
  ENDIF   !(iflag_pbl>1 .AND. iflag_wake>=1 .AND. iflag_pbl_split >=1 )

!==================================
!  thermiques, poches, convection
!==================================

! Emanuel
  found=phyetat0_get(klev,sig1,"sig1","sig1",0.)
  found=phyetat0_get(klev,w01,"w01","w01",0.)

! Wake
  found=phyetat0_get(klev,wake_deltat,"WAKE_DELTAT","Delta T wake/env",0.)
  found=phyetat0_get(klev,wake_deltaq,"WAKE_DELTAQ","Delta hum. wake/env",0.)
  found=phyetat0_get(1,wake_s,"WAKE_S","Wake frac. area",0.)
!jyg<
!  Set wake_dens to -1000. when there is no restart so that the actual
!  initialization is made in calwake.
!!  found=phyetat0_get(1,wake_dens,"WAKE_DENS","Wake num. /unit area",0.)
  found=phyetat0_get(1,wake_dens,"WAKE_DENS","Wake num. /unit area",-1000.)
!>jyg
  found=phyetat0_get(1,wake_cstar,"WAKE_CSTAR","WAKE_CSTAR",0.)
  found=phyetat0_get(1,wake_pe,"WAKE_PE","WAKE_PE",0.)
  found=phyetat0_get(1,wake_fip,"WAKE_FIP","WAKE_FIP",0.)

! Thermiques
  found=phyetat0_get(1,zmax0,"ZMAX0","ZMAX0",40.)
  found=phyetat0_get(1,f0,"F0","F0",1.e-5)
  found=phyetat0_get(klev+1,fm_therm,"FM_THERM","Thermals mass flux",0.)
  found=phyetat0_get(klev,entr_therm,"ENTR_THERM","Thermals Entrain.",0.)
  found=phyetat0_get(klev,detr_therm,"DETR_THERM","Thermals Detrain.",0.)

! ALE/ALP
  found=phyetat0_get(1,ale_bl,"ALE_BL","ALE BL",0.)
  found=phyetat0_get(1,ale_bl_trig,"ALE_BL_TRIG","ALE BL_TRIG",0.)
  found=phyetat0_get(1,alp_bl,"ALP_BL","ALP BL",0.)
  found=phyetat0_get(1,ale_wake,"ALE_WAKE","ALE_WAKE",0.)
  found=phyetat0_get(1,ale_bl_stat,"ALE_BL_STAT","ALE_BL_STAT",0.)

!===========================================
  ! Read and send field trs to traclmdz
!===========================================

  IF (type_trac == 'lmdz') THEN
     DO it=1, nbtr                                                                  
!!        iiq=niadv(it+2)                                                           ! jyg
        iiq=niadv(it+nqo)                                                           ! jyg
        found=phyetat0_get(1,trs(:,it),"trs_"//tname(iiq), &
              "Surf trac"//tname(iiq),0.)
     END DO
     CALL traclmdz_from_restart(trs)

     IF (carbon_cycle_cpl) THEN
        ALLOCATE(co2_send(klon), stat=ierr)
        IF (ierr /= 0) CALL abort_physic('phyetat0', 'pb allocation co2_send', 1)
        found=phyetat0_get(1,co2_send,"co2_send","co2 send",0.)
     END IF
  END IF

!===========================================
!  ondes de gravite / relief
!===========================================

!  ondes de gravite non orographiques
  if (ok_gwd_rando) found = &
       phyetat0_get(klev,du_gwd_rando,"du_gwd_rando","du_gwd_rando",0.)
  IF (.not. ok_hines .and. ok_gwd_rando) found &
       = phyetat0_get(klev,du_gwd_front,"du_gwd_front","du_gwd_front",0.)

!  prise en compte du relief sous-maille
  found=phyetat0_get(1,zmea,"ZMEA","sub grid orography",0.)
  found=phyetat0_get(1,zstd,"ZSTD","sub grid orography",0.)
  found=phyetat0_get(1,zsig,"ZSIG","sub grid orography",0.)
  found=phyetat0_get(1,zgam,"ZGAM","sub grid orography",0.)
  found=phyetat0_get(1,zthe,"ZTHE","sub grid orography",0.)
  found=phyetat0_get(1,zpic,"ZPIC","sub grid orography",0.)
  found=phyetat0_get(1,zval,"ZVAL","sub grid orography",0.)
  found=phyetat0_get(1,zmea,"ZMEA","sub grid orography",0.)
  found=phyetat0_get(1,rugoro,"RUGSREL","sub grid orography",0.)

!===========================================
! Initialize ocean
!===========================================

  IF ( type_ocean == 'slab' ) THEN
      CALL ocean_slab_init(dtime, pctsrf)
      IF (nslay.EQ.1) THEN
        found=phyetat0_get(1,tslab,"tslab01","tslab",0.)
        IF (.NOT. found) THEN
            found=phyetat0_get(1,tslab,"tslab","tslab",0.)
        END IF
      ELSE
          DO i=1,nslay
            WRITE(str2,'(i2.2)') i
            found=phyetat0_get(1,tslab(:,i),"tslab"//str2,"tslab",0.)  
          END DO
      END IF
      IF (.NOT. found) THEN 
          PRINT*, "phyetat0: Le champ <tslab> est absent"
          PRINT*, "Initialisation a tsol_oce"
          DO i=1,nslay
              tslab(:,i)=MAX(ftsol(:,is_oce),271.35)
          END DO
      END IF 

      ! Sea ice variables
      IF (version_ocean == 'sicINT') THEN
          found=phyetat0_get(1,tice,"slab_tice","slab_tice",0.)
          IF (.NOT. found) THEN 
              PRINT*, "phyetat0: Le champ <tice> est absent"
              PRINT*, "Initialisation a tsol_sic"
                  tice(:)=ftsol(:,is_sic)
          END IF 
          found=phyetat0_get(1,seaice,"seaice","seaice",0.)
          IF (.NOT. found) THEN
              PRINT*, "phyetat0: Le champ <seaice> est absent"
              PRINT*, "Initialisation a 0/1m suivant fraction glace"
              seaice(:)=0.
              WHERE (pctsrf(:,is_sic).GT.EPSFRA)
                  seaice=917.
              END WHERE
          END IF
      END IF !sea ice INT
  END IF ! Slab        

  ! on ferme le fichier
  CALL close_startphy

  ! Initialize module pbl_surface_mod 

  CALL pbl_surface_init(fder, snow, qsurf, tsoil)

  ! Initialize module ocean_cpl_mod for the case of coupled ocean
  IF ( type_ocean == 'couple' ) THEN
     CALL ocean_cpl_init(dtime, longitude_deg, latitude_deg)
  ENDIF

  CALL init_iophy_new(latitude_deg, longitude_deg)

  ! Initilialize module fonte_neige_mod      
  CALL fonte_neige_init(run_off_lic_0)

END SUBROUTINE phyetat0

!===================================================================
FUNCTION phyetat0_get(nlev,field,name,descr,default)
!===================================================================
! Lecture d'un champ avec contrôle
! Function logique dont le resultat indique si la lecture
! s'est bien passée
! On donne une valeur par defaut dans le cas contraire
!===================================================================

USE iostart, ONLY : get_field
USE dimphy, only: klon 
USE print_control_mod, ONLY: lunout

IMPLICIT NONE

LOGICAL phyetat0_get

! arguments
INTEGER,INTENT(IN) :: nlev
CHARACTER*(*),INTENT(IN) :: name,descr
REAL,INTENT(IN) :: default
REAL,DIMENSION(klon,nlev),INTENT(INOUT) :: field

! Local variables
LOGICAL found

   CALL get_field(name, field, found)
   IF (.NOT. found) THEN
     WRITE(lunout,*) "phyetat0: Le champ <",name,"> est absent"
     WRITE(lunout,*) "Depart legerement fausse. Mais je continue"
     field(:,:)=default
   ENDIF
   WRITE(lunout,*) name, descr, MINval(field),MAXval(field)
   phyetat0_get=found

RETURN
END FUNCTION phyetat0_get

!================================================================
FUNCTION phyetat0_srf(nlev,field,name,descr,default)
!===================================================================
! Lecture d'un champ par sous-surface avec contrôle
! Function logique dont le resultat indique si la lecture
! s'est bien passée
! On donne une valeur par defaut dans le cas contraire
!===================================================================

USE iostart, ONLY : get_field
USE dimphy, only: klon 
USE indice_sol_mod, only: nbsrf
USE print_control_mod, ONLY: lunout

IMPLICIT NONE

LOGICAL phyetat0_srf
! arguments
INTEGER,INTENT(IN) :: nlev
CHARACTER*(*),INTENT(IN) :: name,descr
REAL,INTENT(IN) :: default
REAL,DIMENSION(klon,nlev,nbsrf),INTENT(INOUT) :: field

! Local variables
LOGICAL found,phyetat0_get
INTEGER nsrf
CHARACTER*2 str2
  
     IF (nbsrf.GT.99) THEN
        WRITE(lunout,*) "Trop de sous-mailles"
        call abort_physic("phyetat0", "", 1)
     ENDIF

     DO nsrf = 1, nbsrf
        WRITE(str2, '(i2.2)') nsrf
        found= phyetat0_get(nlev,field(:,:, nsrf), &
        name//str2,descr//" srf:"//str2,default)
     ENDDO

     phyetat0_srf=found

RETURN
END FUNCTION phyetat0_srf

