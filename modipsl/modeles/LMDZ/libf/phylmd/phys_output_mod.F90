! $Id: phys_output_mod.F90 3284 2018-03-17 08:58:45Z fairhead $
!

MODULE phys_output_mod 
  USE indice_sol_mod
  USE phys_output_var_mod
  USE phys_output_write_mod, ONLY : phys_output_write
  REAL, DIMENSION(nfiles),SAVE :: ecrit_files

! Abderrahmane 12 2007
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Ecreture des Sorties du modele dans les fichiers Netcdf :
! histmth.nc : moyennes mensuelles
! histday.nc : moyennes journalieres
! histhf.nc  : moyennes toutes les 3 heures
! histins.nc : valeurs instantanees
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Ouverture des fichier et definition des variable de sortie !!!!!!!!
  !! histbeg, histvert et histdef
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  SUBROUTINE phys_output_open(rlon,rlat,pim,tabij,ipt,jpt,plon,plat, &
       jjmp1,nlevSTD,clevSTD,rlevSTD, dtime, ok_veget, &
       type_ocean, iflag_pbl,iflag_pbl_split,ok_mensuel,ok_journe, &
       ok_hf,ok_instan,ok_LES,ok_ade,ok_aie, read_climoz, &
       phys_out_filestations, &
       new_aod, aerosol_couple, flag_aerosol_strat, &
       pdtphys, paprs, pphis, pplay, lmax_th, ptconv, ptconvth, ivap, &
       d_u, d_t, qx, d_qx, zmasse, ok_sync)   

    USE iophy 
    USE dimphy
    USE infotrac_phy, ONLY: nqtot, nqo, niadv, tname, ttext, type_trac
    USE ioipsl
    USE phys_cal_mod, only : hour, calend
    USE mod_phys_lmdz_para
    !Martin
    USE surface_data, ONLY : ok_snow
    USE phys_output_ctrlout_mod
    USE mod_grid_phy_lmdz, only: klon_glo,nbp_lon,nbp_lat
    USE print_control_mod, ONLY: prt_level,lunout
    USE vertical_layers_mod, ONLY: ap,bp,preff,presnivs, aps, bps, pseudoalt
    USE time_phylmdz_mod, ONLY: day_ini, itau_phy, start_time, annee_ref, day_ref
#ifdef REPROBUS
    USE chem_rep, ONLY: nbnas, tnamenas, ttextnas
#endif
#ifdef CPP_XIOS
    ! ug Pour les sorties XIOS
    USE wxios
#endif

    IMPLICIT NONE
    include "clesphys.h"
    include "thermcell.h"
    include "YOMCST.h"

    ! ug Nouveaux arguments n\'ecessaires au histwrite_mod:
    INTEGER, INTENT(IN)                         :: ivap
    INTEGER, DIMENSION(klon), INTENT(IN)        :: lmax_th
    LOGICAL, INTENT(IN)                         :: ok_sync
    LOGICAL, DIMENSION(klon, klev), INTENT(IN)  :: ptconv, ptconvth
    REAL, INTENT(IN)                            :: pdtphys
    REAL, DIMENSION(klon), INTENT(IN)           :: pphis
    REAL, DIMENSION(klon, klev), INTENT(IN)     :: pplay, d_u, d_t
    REAL, DIMENSION(klon, klev+1), INTENT(IN)   :: paprs
    REAL, DIMENSION(klon,klev,nqtot), INTENT(IN):: qx, d_qx
    REAL, DIMENSION(klon, klev), INTENT(IN)      :: zmasse


    REAL,DIMENSION(klon),INTENT(IN) :: rlon
    REAL,DIMENSION(klon),INTENT(IN) :: rlat
    INTEGER, INTENT(IN)             :: pim
    INTEGER, DIMENSION(pim)            :: tabij
    INTEGER,DIMENSION(pim), INTENT(IN) :: ipt, jpt
    REAL,DIMENSION(pim), INTENT(IN) :: plat, plon
    REAL,DIMENSION(pim,2) :: plat_bounds, plon_bounds

    INTEGER                               :: jjmp1
    INTEGER                               :: nlevSTD, radpas
    LOGICAL                               :: ok_mensuel, ok_journe, ok_hf, ok_instan
    LOGICAL                               :: ok_LES,ok_ade,ok_aie
    INTEGER                               :: flag_aerosol_strat
    LOGICAL                               :: new_aod, aerosol_couple
    INTEGER, INTENT(IN)::  read_climoz ! read ozone climatology
    !     Allowed values are 0, 1 and 2
    !     0: do not read an ozone climatology
    !     1: read a single ozone climatology that will be used day and night
    !     2: read two ozone climatologies, the average day and night
    !     climatology and the daylight climatology

    REAL                                  :: dtime
    INTEGER                               :: idayref
    REAL                                  :: zjulian_start, zjulian
    CHARACTER(LEN=4), DIMENSION(nlevSTD)  :: clevSTD
    REAL, DIMENSION(nlevSTD)              :: rlevSTD
    INTEGER                               :: nsrf, k, iq, iiq, iff, i, j, ilev
    INTEGER                               :: naero
    LOGICAL                               :: ok_veget
    INTEGER                               :: iflag_pbl
    INTEGER                               :: iflag_pbl_split
    CHARACTER(LEN=4)                      :: bb2
    CHARACTER(LEN=2)                      :: bb3
    CHARACTER(LEN=6)                      :: type_ocean
    INTEGER, DIMENSION(nbp_lon*jjmp1)         ::  ndex2d
    INTEGER, DIMENSION(nbp_lon*jjmp1*klev)    :: ndex3d
    INTEGER                               :: imin_ins, imax_ins
    INTEGER                               :: jmin_ins, jmax_ins
    INTEGER, DIMENSION(nfiles)            :: phys_out_levmin, phys_out_levmax
    INTEGER, DIMENSION(nfiles)            :: phys_out_filelevels
    CHARACTER(LEN=20), DIMENSION(nfiles)  :: chtimestep = (/ 'Default', 'Default', 'Default', 'Default', 'Default', &
                                                             'Default', 'Default', 'Default', 'Default', 'Default' /)
    LOGICAL, DIMENSION(nfiles)            :: phys_out_filekeys
    LOGICAL, DIMENSION(nfiles)            :: phys_out_filestations

!!!!!!!!!! stockage dans une region limitee pour chaque fichier !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                 entre [phys_out_lonmin,phys_out_lonmax] et [phys_out_latmin,phys_out_latmax]

    LOGICAL, DIMENSION(nfiles), SAVE  :: phys_out_regfkey       = (/ .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
                                                                     .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /)
    REAL, DIMENSION(nfiles), SAVE     :: phys_out_lonmin        = (/ -180., -180., -180., -180., -180., &
                                                                     -180., -180., -180., -180., -180. /)
    REAL, DIMENSION(nfiles), SAVE     :: phys_out_lonmax        = (/  180.,  180.,  180.,  180.,  180., &
                                                                      180.,  180.,  180.,  180.,  180. /)
    REAL, DIMENSION(nfiles), SAVE     :: phys_out_latmin        = (/  -90.,  -90.,  -90.,  -90.,  -90., &
                                                                      -90.,  -90.,  -90.,  -90.,  -90. /)
    REAL, DIMENSION(nfiles), SAVE     :: phys_out_latmax        = (/   90.,   90.,   90.,   90.,   90., &
                                                                       90.,   90.,   90.,   90.,   90. /)
    REAL, DIMENSION(klev,2) :: Ahyb_bounds, Bhyb_bounds
    REAL, DIMENSION(klev+1)   :: lev_index
                
#ifdef CPP_XIOS
    ! ug Variables utilis\'ees pour r\'ecup\'erer le calendrier pour xios
    INTEGER :: x_an, x_mois, x_jour
    REAL :: x_heure
    INTEGER :: ini_an, ini_mois, ini_jour
    REAL :: ini_heure
#endif
    INTEGER                         :: ISW
    REAL, DIMENSION(NSW)            :: wl1_sun, wl2_sun !wavelength bounds (in um) for SW
    REAL, DIMENSION(NSW)            :: wn1_sun, wn2_sun !wavenumber bounds (in m-1) for SW 
    REAL, DIMENSION(NSW)            :: spectband  !mean wavenumb. of each sp.band
    REAL, DIMENSION(NSW,2)          :: spbnds_sun !bounds of spectband

    WRITE(lunout,*) 'Debut phys_output_mod.F90'
    ! Initialisations (Valeurs par defaut

    DO ilev=1,klev
      Ahyb_bounds(ilev,1) = ap(ilev)
      Ahyb_bounds(ilev,2) = ap(ilev+1)
      Bhyb_bounds(ilev,1) = bp(ilev)
      Bhyb_bounds(ilev,2) = bp(ilev+1)
      lev_index(ilev) = REAL(ilev)
    END DO
      lev_index(klev+1) = REAL(klev+1)

    IF (.NOT. ALLOCATED(o_trac)) ALLOCATE(o_trac(nqtot))
    IF (.NOT. ALLOCATED(o_trac_cum)) ALLOCATE(o_trac_cum(nqtot))
#ifdef REPROBUS
    IF (.NOT. ALLOCATED(o_nas)) ALLOCATE(o_nas(nbnas))
#endif
    ALLOCATE(o_dtr_the(nqtot),o_dtr_con(nqtot),o_dtr_lessi_impa(nqtot))
    ALLOCATE(o_dtr_lessi_nucl(nqtot),o_dtr_insc(nqtot),o_dtr_bcscav(nqtot))
    ALLOCATE(o_dtr_evapls(nqtot),o_dtr_ls(nqtot),o_dtr_trsp(nqtot))
    ALLOCATE(o_dtr_sscav(nqtot),o_dtr_sat(nqtot),o_dtr_uscav(nqtot))
    ALLOCATE(o_dtr_dry(nqtot),o_dtr_vdf(nqtot))

    levmax = (/ klev, klev, klev, klev, klev, klev, nlevSTD, nlevSTD, nlevSTD, klev /)

    phys_out_filenames(1) = 'histmth'
    phys_out_filenames(2) = 'histday'
    phys_out_filenames(3) = 'histhf6h'
    phys_out_filenames(4) = 'histhf3h'
    phys_out_filenames(5) = 'histhf3hm'
    phys_out_filenames(6) = 'histstn'
    phys_out_filenames(7) = 'histmthNMC'
    phys_out_filenames(8) = 'histdayNMC'
    phys_out_filenames(9) = 'histhfNMC'
    phys_out_filenames(10)= 'histstrataer'

    type_ecri(1) = 'ave(X)'
    type_ecri(2) = 'ave(X)'
    type_ecri(3) = 'inst(X)'
    type_ecri(4) = 'inst(X)'
    type_ecri(5) = 'ave(X)'
    type_ecri(6) = 'inst(X)'
    type_ecri(7) = 'inst(X)'
    type_ecri(8) = 'inst(X)'
    type_ecri(9) = 'inst(X)'
    type_ecri(10)= 'ave(X)'

    clef_files(1) = ok_mensuel
    clef_files(2) = ok_journe
    clef_files(3) = ok_hf
    clef_files(4) = ok_instan
    clef_files(5) = ok_LES
    clef_files(6) = ok_instan
    clef_files(7) = ok_histNMC(1)
    clef_files(8) = ok_histNMC(2)
    clef_files(9) = ok_histNMC(3)
#ifdef CPP_StratAer
    clef_files(10)= .TRUE.
#else
    clef_files(10)= .FALSE.
#endif

    !sortir des fichiers "stations" si clef_stations(:)=.TRUE.
    clef_stations(1) = .FALSE.
    clef_stations(2) = .FALSE.
    clef_stations(3) = .FALSE.
    clef_stations(4) = .FALSE.
    clef_stations(5) = .FALSE.
    clef_stations(6) = .FALSE.
    clef_stations(7) = .FALSE.
    clef_stations(8) = .FALSE.
    clef_stations(9) = .FALSE.
    clef_stations(10)= .FALSE.

    lev_files(1) = lev_histmth
    lev_files(2) = lev_histday
    lev_files(3) = lev_histhf
    lev_files(4) = lev_histins
    lev_files(5) = lev_histLES
    lev_files(6) = lev_histins
    lev_files(7) = levout_histNMC(1)
    lev_files(8) = levout_histNMC(2)
    lev_files(9) = levout_histNMC(3)
    lev_files(10)= 5

    ecrit_files(1) = ecrit_mth
    ecrit_files(2) = ecrit_day
    ecrit_files(3) = ecrit_hf
    ecrit_files(4) = ecrit_ins
    ecrit_files(5) = ecrit_LES
    ecrit_files(6) = ecrit_ins
    ecrit_files(7) = freq_outNMC(1)
    ecrit_files(8) = freq_outNMC(2)
    ecrit_files(9) = freq_outNMC(3)
    ecrit_files(10)= ecrit_mth

    !! Lectures des parametres de sorties dans physiq.def

    CALL getin('phys_out_regfkey',phys_out_regfkey)
    CALL getin('phys_out_lonmin',phys_out_lonmin)
    CALL getin('phys_out_lonmax',phys_out_lonmax)
    CALL getin('phys_out_latmin',phys_out_latmin)
    CALL getin('phys_out_latmax',phys_out_latmax)
    phys_out_levmin(:)=levmin(:)
    CALL getin('phys_out_levmin',levmin)
    phys_out_levmax(:)=levmax(:)
    CALL getin('phys_out_levmax',levmax)
    CALL getin('phys_out_filenames',phys_out_filenames)
    phys_out_filekeys(:)=clef_files(:)
    CALL getin('phys_out_filekeys',clef_files)
    phys_out_filestations(:)=clef_stations(:)
    CALL getin('phys_out_filestations',clef_stations)
    phys_out_filelevels(:)=lev_files(:)
    CALL getin('phys_out_filelevels',lev_files)
    CALL getin('phys_out_filetimesteps',chtimestep)
    phys_out_filetypes(:)=type_ecri(:)
    CALL getin('phys_out_filetypes',type_ecri)

    type_ecri_files(:)=type_ecri(:)

!    if (ok_all_xml) phys_out_filelevels = 999

    WRITE(lunout,*)'phys_out_lonmin=',phys_out_lonmin
    WRITE(lunout,*)'phys_out_lonmax=',phys_out_lonmax
    WRITE(lunout,*)'phys_out_latmin=',phys_out_latmin
    WRITE(lunout,*)'phys_out_latmax=',phys_out_latmax
    WRITE(lunout,*)'phys_out_filenames=',phys_out_filenames
    WRITE(lunout,*)'phys_out_filetypes=',type_ecri
    WRITE(lunout,*)'phys_out_filekeys=',clef_files
    WRITE(lunout,*)'phys_out_filestations=',clef_stations
    WRITE(lunout,*)'phys_out_filelevels=',lev_files
    WRITE(lunout,*)'phys_out_regfkey=',phys_out_regfkey

! A noter pour
! l heure initiale - dans les fichiers histoire hist* - on met comme  
! heure de debut soit la vraie heure (pour le 1D) soit 0h (pour le 3D) 
! afin d avoir une seule sortie mensuelle par mois lorsque l on tourne 
! par annee (IM).
!
     idayref = day_ref
     IF (klon_glo==1) THEN
       ! current_time (used to compute hour) is updated at the begining of
       ! the physics; to set the correct outputs "initial time" we thus 
       ! have to use (hour-dtphys).
         CALL ymds2ju(annee_ref, 1, idayref, hour-pdtphys, zjulian)
         print *,'phys_output_mod: annee,iday,hour,zjulian=',annee_ref,idayref, hour, zjulian
     ELSE
         CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
         CALL ymds2ju(annee_ref, 1, day_ini, start_time*rday, zjulian_start)
     ENDIF

#ifdef CPP_XIOS
    ! ug R\'eglage du calendrier xios
    !Temps julian => an, mois, jour, heure
    CALL ju2ymds(zjulian, x_an, x_mois, x_jour, x_heure)
    CALL ju2ymds(zjulian_start, ini_an, ini_mois, ini_jour, ini_heure)
    CALL wxios_set_cal(dtime, calend, x_an, x_mois, x_jour, x_heure, ini_an, &
                       ini_mois, ini_jour, ini_heure )
#endif

!!!!!!!!!!!!!!!!!!!!!!! Boucle sur les fichiers !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Appel de histbeg et histvert pour creer le fichier et les niveaux verticaux !!
    ! Appel des histbeg pour definir les variables (nom, moy ou inst, freq de sortie ..
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    zdtime_moy = dtime         ! Frequence ou l on moyenne


  ecrit_files(7) = ecrit_files(1)
  ecrit_files(8) = ecrit_files(2)
  ecrit_files(9) = ecrit_files(3)

  DO iff=1,nfiles

       ! Calculate ecrit_files for all files
      IF ( chtimestep(iff).eq.'Default' ) THEN
          ! Par defaut ecrit_files = (ecrit_mensuel ecrit_jour ecrit_hf
          ! ...)*86400.
          ecrit_files(iff)=ecrit_files(iff)*86400.
      ELSE IF (chtimestep(iff).eq.'-1') THEN
          PRINT*,'ecrit_files(',iff,') < 0 so IOIPSL work on different'
          PRINT*,'months length'
          ecrit_files(iff)=-1.
      ELSE
       CALL convers_timesteps(chtimestep(iff),dtime,ecrit_files(iff)) 
      ENDIF

       WRITE(lunout,*)'ecrit_files(',iff,')= ',ecrit_files(iff)
       zoutm(iff) = ecrit_files(iff) ! Frequence ou l on ecrit en seconde


#ifdef CPP_XIOS
!!! Ouverture de chaque fichier XIOS !!!!!!!!!!!
    IF (.not. ok_all_xml) THEN
      IF (prt_level >= 10) THEN
        print*,'phys_output_open: call wxios_add_file with phys_out_filenames(iff)=',trim(phys_out_filenames(iff))    
      ENDIF
      CALL wxios_add_file(phys_out_filenames(iff),chtimestep(iff),lev_files(iff))  
    ENDIF

!!! Declaration des axes verticaux de chaque fichier:
    IF (prt_level >= 10) THEN
      print*,'phys_output_open: Declare vertical axes for each file'
    ENDIF

   IF (iff.LE.6.OR.iff.EQ.10) THEN
    CALL wxios_add_vaxis("presnivs", &
            levmax(iff) - levmin(iff) + 1, presnivs(levmin(iff):levmax(iff)))
    CALL wxios_add_vaxis("Ahyb", &
            levmax(iff) - levmin(iff) + 1, aps(levmin(iff):levmax(iff)), positif='down', &
            bnds=Ahyb_bounds(levmin(iff):levmax(iff),:))
    CALL wxios_add_vaxis("Bhyb", &
            levmax(iff) - levmin(iff) + 1, bps(levmin(iff):levmax(iff)), positif='down', &
            bnds=Bhyb_bounds(levmin(iff):levmax(iff),:))
    CALL wxios_add_vaxis("klev", levmax(iff) - levmin(iff) + 1, &
                          lev_index(levmin(iff):levmax(iff)))
    CALL wxios_add_vaxis("klevp1", klev+1, &
                          lev_index(1:klev+1))
    CALL wxios_add_vaxis("bnds", 2, (/1.,2./))

     CALL wxios_add_vaxis("Alt", &
            levmax(iff) - levmin(iff) + 1, pseudoalt)

    IF (NSW.EQ.6) THEN
!
!wl1_sun: minimum bound of wavelength (in um)
!
      wl1_sun(1)=0.180
      wl1_sun(2)=0.250
      wl1_sun(3)=0.440
      wl1_sun(4)=0.690
      wl1_sun(5)=1.190
      wl1_sun(6)=2.380
!
!wl2_sun: maximum bound of wavelength (in um)
!
      wl2_sun(1)=0.250
      wl2_sun(2)=0.440
      wl2_sun(3)=0.690
      wl2_sun(4)=1.190
      wl2_sun(5)=2.380
      wl2_sun(6)=4.000
!
    ELSE IF(NSW.EQ.2) THEN
!
!wl1_sun: minimum bound of wavelength (in um)
!
      wl1_sun(1)=0.250
      wl1_sun(2)=0.690
!
!wl2_sun: maximum bound of wavelength (in um)
!
      wl2_sun(1)=0.690
      wl2_sun(2)=4.000
    ENDIF

    DO ISW=1, NSW
     wn1_sun(ISW)=1.e+6/wl1_sun(ISW) 
     wn2_sun(ISW)=1.e+6/wl2_sun(ISW) 
     spbnds_sun(ISW,1)=wn2_sun(ISW)
     spbnds_sun(ISW,2)=wn1_sun(ISW)
     spectband(ISW)=(wn1_sun(ISW)+wn2_sun(ISW))/2
    ENDDO
!
!!! ajout axe vertical spectband : solar band number
    CALL wxios_add_vaxis("spectband", NSW, spectband, positif='down')
   ELSE
    ! NMC files
    CALL wxios_add_vaxis("plev", &
            levmax(iff) - levmin(iff) + 1, rlevSTD(levmin(iff):levmax(iff)))
   ENDIF
#endif 

        IF (clef_files(iff)) THEN
!!!!!!!!!!!!!!!!! Traitement dans le cas ou l'on veut stocker sur un domaine limite !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (phys_out_regfkey(iff)) THEN
             imin_ins=1
             imax_ins=nbp_lon
             jmin_ins=1
             jmax_ins=jjmp1

             ! correction abderr        
             DO i=1,nbp_lon
                WRITE(lunout,*)'io_lon(i)=',io_lon(i)
                IF (io_lon(i).le.phys_out_lonmin(iff)) imin_ins=i
                IF (io_lon(i).le.phys_out_lonmax(iff)) imax_ins=i+1
             ENDDO

             DO j=1,jjmp1
                WRITE(lunout,*)'io_lat(j)=',io_lat(j)
                IF (io_lat(j).ge.phys_out_latmin(iff)) jmax_ins=j+1
                IF (io_lat(j).ge.phys_out_latmax(iff)) jmin_ins=j
             ENDDO

             WRITE(lunout,*)'On stoke le fichier histoire numero ',iff,' sur ', &
                  imin_ins,imax_ins,jmin_ins,jmax_ins
             WRITE(lunout,*)'longitudes : ', &
                  io_lon(imin_ins),io_lon(imax_ins), &
                  'latitudes : ', &
                  io_lat(jmax_ins),io_lat(jmin_ins)

             CALL histbeg(phys_out_filenames(iff),nbp_lon,io_lon,jjmp1,io_lat, &
                  imin_ins,imax_ins-imin_ins+1, &
                  jmin_ins,jmax_ins-jmin_ins+1, &
                  itau_phy,zjulian,dtime,nhorim(iff),nid_files(iff))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !IM fichiers stations
          ELSE IF (clef_stations(iff)) THEN

             IF (prt_level >= 10) THEN
             WRITE(lunout,*)'phys_output_open: iff=',iff,'  phys_out_filenames(iff)=',phys_out_filenames(iff)
             ENDIF
             
             CALL histbeg_phy_all(rlon,rlat,pim,tabij,ipt,jpt,plon,plat,plon_bounds,plat_bounds, &
                  phys_out_filenames(iff), &
                  itau_phy,zjulian,dtime,nhorim(iff),nid_files(iff))
          ELSE

             IF (prt_level >= 10) THEN
             WRITE(lunout,*)'phys_output_open: iff=',iff,'  phys_out_filenames(iff)=',phys_out_filenames(iff)
             ENDIF

             CALL histbeg_phy_all(phys_out_filenames(iff),itau_phy,zjulian,&
                 dtime,nhorim(iff),nid_files(iff))
          ENDIF

#ifndef CPP_IOIPSL_NO_OUTPUT 
          IF (iff.LE.6.OR.iff.EQ.10) THEN
             CALL histvert(nid_files(iff), "presnivs", "Vertical levels", "Pa", &  
               levmax(iff) - levmin(iff) + 1, &
               presnivs(levmin(iff):levmax(iff)), nvertm(iff),"down")
!!!! Composantes de la coordonnee sigma-hybride 
          CALL histvert(nid_files(iff), "Ahyb","Ahyb comp of Hyb Cord ", "Pa", &
               levmax(iff) - levmin(iff) + 1,aps,nvertap(iff))

          CALL histvert(nid_files(iff), "Bhyb","Bhyb comp of Hyb Cord", " ", &
               levmax(iff) - levmin(iff) + 1,bps,nvertbp(iff))

          CALL histvert(nid_files(iff), "Alt","Height approx for scale heigh of 8km at levels", "Km", &                       
               levmax(iff) - levmin(iff) + 1,pseudoalt,nvertAlt(iff))

          ELSE
          ! NMC files
             CALL histvert(nid_files(iff), "plev", "pressure", "Pa", &
               levmax(iff) - levmin(iff) + 1, &
              rlevSTD(levmin(iff):levmax(iff)), nvertm(iff), "down")
          ENDIF
#endif

     ENDIF ! clef_files

       IF (nqtot>=nqo+1) THEN
!
            DO iq=nqo+1,nqtot 
            iiq=niadv(iq)
            o_trac(iq-nqo) = ctrl_out((/ 1, 5, 5, 5, 10, 10, 11, 11, 11, 11 /), &
                           tname(iiq),'Tracer '//ttext(iiq), "-",  &
                           (/ '', '', '', '', '', '', '', '', '', '' /))
            o_dtr_vdf(iq-nqo) = ctrl_out((/ 4, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                              'd'//trim(tname(iq))//'_vdf',  &
                              'Tendance tracer '//ttext(iiq), "-" , &
                              (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_the(iq-nqo) = ctrl_out((/ 5, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                              'd'//trim(tname(iq))//'_the', &
                              'Tendance tracer '//ttext(iiq), "-", &
                              (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_con(iq-nqo) = ctrl_out((/ 5, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                              'd'//trim(tname(iq))//'_con', &
                              'Tendance tracer '//ttext(iiq), "-", &
                              (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_lessi_impa(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                                     'd'//trim(tname(iq))//'_lessi_impa', &
                                     'Tendance tracer '//ttext(iiq), "-", &
                                     (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_lessi_nucl(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                                     'd'//trim(tname(iq))//'_lessi_nucl', &
                                     'Tendance tracer '//ttext(iiq), "-", &
                                     (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_insc(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                               'd'//trim(tname(iq))//'_insc', &
                               'Tendance tracer '//ttext(iiq), "-", &
                               (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_bcscav(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                                 'd'//trim(tname(iq))//'_bcscav', &
                                 'Tendance tracer '//ttext(iiq), "-", &
                                 (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_evapls(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                                 'd'//trim(tname(iq))//'_evapls', &
                                 'Tendance tracer '//ttext(iiq), "-", &
                                 (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_ls(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                             'd'//trim(tname(iq))//'_ls', &
                             'Tendance tracer '//ttext(iiq), "-", &
                             (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_trsp(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                               'd'//trim(tname(iq))//'_trsp', &
                               'Tendance tracer '//ttext(iiq), "-", &
                               (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_sscav(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                                'd'//trim(tname(iq))//'_sscav', &
                                'Tendance tracer '//ttext(iiq), "-", &
                                (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_sat(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                               'd'//trim(tname(iq))//'_sat', &
                               'Tendance tracer '//ttext(iiq), "-", &
                               (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_uscav(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                                'd'//trim(tname(iq))//'_uscav', &
                                'Tendance tracer '//ttext(iiq), "-", &
                                 (/ '', '', '', '', '', '', '', '', '', '' /))

            o_dtr_dry(iq-nqo) = ctrl_out((/ 7, 7, 7, 7, 10, 10, 11, 11, 11, 11 /), &
                              'cum'//'d'//trim(tname(iq))//'_dry', &
                              'tracer tendency dry deposition'//ttext(iiq), "-", &
                              (/ '', '', '', '', '', '', '', '', '', '' /))

            o_trac_cum(iq-nqo) = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11 /), &
                               'cum'//tname(iiq),&
                               'Cumulated tracer '//ttext(iiq), "-", &
                               (/ '', '', '', '', '', '', '', '', '', '' /))
            ENDDO
      ENDIF
      IF (type_trac=='repr') THEN
#ifdef REPROBUS
         DO iiq=1,nbnas
            o_nas(iiq) = ctrl_out((/ 4, 5, 5, 5, 10, 10, 11, 11, 11, 11 /), &
                 tnamenas(iiq),ttextnas(iiq), "-", &
                 (/ '', '', '', '', '', '', '', '', '', '' /))
         ENDDO
#endif
      ENDIF

   ENDDO !  iff

    ! Updated write frequencies due to phys_out_filetimesteps. 
    ! Write frequencies are now in seconds.  
    ecrit_mth = ecrit_files(1)
    ecrit_day = ecrit_files(2)
    ecrit_hf  = ecrit_files(3)
    ecrit_ins = ecrit_files(4)
    ecrit_LES = ecrit_files(5)
    ecrit_ins = ecrit_files(6)

    IF (prt_level >= 10) THEN
      WRITE(lunout,*)'swaerofree_diag=',swaerofree_diag
      WRITE(lunout,*)'swaero_diag=',swaero_diag
      WRITE(lunout,*)'dryaod_diag=',dryaod_diag
      WRITE(lunout,*)'ok_4xCO2atm=',ok_4xCO2atm
      WRITE(lunout,*)'phys_output_open: ends here'
    ENDIF

  END SUBROUTINE phys_output_open

  SUBROUTINE convers_timesteps(str,dtime,timestep)

    use ioipsl
    USE phys_cal_mod
    USE time_phylmdz_mod, ONLY: day_ref, annee_ref
    USE print_control_mod, ONLY: lunout

    IMPLICIT NONE

    CHARACTER(LEN=20)   :: str
    CHARACTER(LEN=10)   :: type
    INTEGER             :: ipos,il
    real                :: ttt,xxx,timestep,dayseconde,dtime
    parameter (dayseconde=86400.)

    ipos=scan(str,'0123456789.',.TRUE.)
    !  
    il=len_trim(str)
    WRITE(lunout,*) "ipos = ", ipos
    WRITE(lunout,*) "il = ", il
    IF (ipos == 0) CALL abort_physic("convers_timesteps", "bad str", 1)
    read(str(1:ipos),*) ttt
    WRITE(lunout,*)ttt
    type=str(ipos+1:il)

    IF ( il == ipos ) THEN
       type='day'
    ENDIF

    IF ( type == 'day'.or.type == 'days'.or.type == 'jours'.or.type == 'jour' ) timestep = ttt * dayseconde
    IF ( type == 'mounths'.or.type == 'mth'.or.type == 'mois' ) THEN
       WRITE(lunout,*)'annee_ref,day_ref mon_len',annee_ref,day_ref,mth_len
       timestep = ttt * dayseconde * mth_len
    ENDIF
    IF ( type == 'hours'.or.type == 'hr'.or.type == 'heurs') timestep = ttt * dayseconde / 24.
    IF ( type == 'mn'.or.type == 'minutes'  ) timestep = ttt * 60.
    IF ( type == 's'.or.type == 'sec'.or.type == 'secondes'   ) timestep = ttt
    IF ( type == 'TS' ) timestep = ttt * dtime

    WRITE(lunout,*)'type =      ',type
    WRITE(lunout,*)'nb j/h/m =  ',ttt
    WRITE(lunout,*)'timestep(s)=',timestep

  END SUBROUTINE convers_timesteps

END MODULE phys_output_mod
