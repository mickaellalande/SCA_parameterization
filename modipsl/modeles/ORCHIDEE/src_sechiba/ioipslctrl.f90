! ================================================================================================================================
!  MODULE       : ioipslctrl
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF          This module contains subroutine for initialisation of IOIPSL history files and restart files
!!
!!\n DESCRIPTION: This module contains subroutine for initialisation of IOIPSL history files and restart files. The subroutines
!!                ioipslctrl_history, ioipslctrl_histstom, ioipslctrl_histstomipcc, ioipslctrl_restini where previously stored in
!!                intersurf module. 
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_sechiba/ioipslctrl.f90 $
!! $Date: 2015-02-19 18:42:48 +0100 (jeu. 19 févr. 2015) $
!! $Revision: 2548 $
!! \n
!_ ================================================================================================================================

MODULE ioipslctrl

  USE IOIPSL
  USE ioipsl_para 
  USE defprec
  USE constantes
  USE time, ONLY : one_day, dt_sechiba
  USE constantes_soil
  USE pft_parameters
  USE thermosoilc, ONLY : thermosoilc_levels
  USE grid 

  IMPLICIT NONE


  LOGICAL, SAVE                    :: ok_histsync             !! Flag activate syncronization of IOIPSL output
  !$OMP THREADPRIVATE(ok_histsync)
   REAL(r_std), SAVE               :: dw                      !! Frequency of history write (sec.)
!$OMP THREADPRIVATE(dw)
  INTEGER(i_std),PARAMETER         :: max_hist_level = 11     !! 

  PRIVATE
  PUBLIC :: ioipslctrl_history, ioipslctrl_histstom, ioipslctrl_histstomipcc, ioipslctrl_restini
  PUBLIC :: dw, max_hist_level, ok_histsync

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE    : ioipslctrl_history
!!
!>\BRIEF	 This subroutine initialize the IOIPSL output files
!! 
!! DESCRIPTION   : This subroutine initialize the IOIPSL output files sechiab_history.nc and sechiba_out_2.nc. It also calls the
!!                 the subroutines ioipslctrl_histstom and ioipslctrl_histstomipcc for initialization of the IOIPSL stomate output files.
!!                 This subroutine was previously called intsurf_history and located in module intersurf.
!!
!! RECENT CHANGE(S): None
!!
!! \n
!_ ================================================================================================================================
  SUBROUTINE ioipslctrl_history(iim, jjm, lon, lat, kindex, kjpindex, istp_old, date0, dt, hist_id, hist2_id, &
       hist_id_stom, hist_id_stom_IPCC)
    
    USE mod_orchidee_para
    !   
    !  This subroutine initialized the history files for the land-surface scheme
    !
    IMPLICIT NONE
    
    INTEGER(i_std), INTENT(in)                  :: iim, jjm  !! Size in x and y of the data to be handeled
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in) :: lon, lat  !! Longitude and latitude of the data points
    INTEGER(i_std),INTENT (in)                            :: kjpindex
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: kindex
    
    INTEGER(i_std), INTENT(in)                  :: istp_old  !! Time step counter
    REAL(r_std), INTENT(in)                     :: date0     !! Julian day at which istp=0
    REAL(r_std), INTENT(in)                     :: dt        !! Time step of the counter in seconds

    INTEGER(i_std), INTENT(out)                 :: hist_id   !! History file identification for SECHIBA
    INTEGER(i_std), INTENT(out)                 :: hist2_id  !! History file 2 identification for SECHIBA
    !! History file identification for STOMATE and IPCC
    INTEGER(i_std), INTENT(out)                 :: hist_id_stom, hist_id_stom_IPCC 
    !
    !  LOCAL
    !
    CHARACTER(LEN=80) :: histname,histname2                    !! Name of history files for SECHIBA
    CHARACTER(LEN=80) :: stom_histname, stom_ipcc_histname     !! Name of history files for STOMATE
    LOGICAL           :: ok_histfile2                 !! Flag to switch on histfile 2 for SECHIBA
    REAL(r_std)       :: dw2                          !! frequency of history write (sec.)
    CHARACTER(LEN=30)   :: flux_op                    !! Operations to be performed on fluxes
    CHARACTER(LEN=40)   :: flux_insec, flux_scinsec   !! Operation in seconds
    INTEGER(i_std)     :: hist_level, hist2_level     !! history output level (default is 10 => maximum output)
    CHARACTER(LEN=40),DIMENSION(max_hist_level) :: &
         & ave, avecels, avescatter, fluxop, &
         & fluxop_scinsec, tmincels, tmaxcels, once, sumscatter  !! The various operation to be performed
    CHARACTER(LEN=40),DIMENSION(max_hist_level) :: &
         & ave2, avecels2, avescatter2, fluxop2, &
         & fluxop_scinsec2, tmincels2, tmaxcels2, once2, sumscatter2  !! The various operation to be performed
    CHARACTER(LEN=80) :: global_attribute              !! for writing attributes in the output files
    INTEGER(i_std)     :: i, jst
    ! SECHIBA AXIS
    INTEGER(i_std)     :: hori_id                      !! ID of the default horizontal longitude and latitude map.
    INTEGER(i_std)     :: vegax_id, laiax_id, solax_id, soltax_id, nobioax_id !! ID's for two vertical coordinates
    INTEGER(i_std)     :: soildiagax_id                !! ID for diagnostic soil levels
    INTEGER(i_std)     :: solayax_id                   !! ID for the vertical axis of the CWRR hydrology 
    INTEGER(i_std)     :: hori_id2                      !! ID of the default horizontal longitude and latitude map.
    INTEGER(i_std)     :: vegax_id2, laiax_id2, solax_id2, soltax_id2, nobioax_id2, albax_id2 !! ID's for two vertical coordinates
    INTEGER(i_std)     :: solayax_id2                   !! ID for the vertical axis of the CWRR hydrology 
    INTEGER(i_std)     :: snowax_id                     !! ID for snow level axis

    ! STOMATE AXIS
    INTEGER(i_std)     :: hist_PFTaxis_id
! deforestation
    INTEGER(i_std)     :: hist_pool_10axis_id
    INTEGER(i_std)     :: hist_pool_100axis_id
    INTEGER(i_std)     :: hist_pool_11axis_id
    INTEGER(i_std)     :: hist_pool_101axis_id
    ! STOMATE IPCC AXIS
    INTEGER(i_std)     :: hist_IPCC_PFTaxis_id
    !
    INTEGER(i_std)                         :: ier
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lon_rect, lat_rect
    !
    REAL(r_std),DIMENSION(nvm)   :: veg
    REAL(r_std),DIMENSION(nlai+1):: indlai
    REAL(r_std),DIMENSION(ngrnd) :: sol
    REAL(r_std),DIMENSION(nstm)  :: soltyp
    REAL(r_std),DIMENSION(nnobio):: nobiotyp
    REAL(r_std),DIMENSION(2)     :: albtyp
    REAL(r_std),DIMENSION(nslm)  :: solay
    REAL(r_std),DIMENSION(nsnow) :: snowlev           !! Layers for snow axis
    !
    CHARACTER(LEN=80)           :: var_name           !! To store variables names
    !
    ! STOMATE history file
    REAL(r_std)                  :: hist_days_stom     !!- GK time step in days for this history file
    REAL(r_std)                  :: hist_dt_stom       !!- GK time step in seconds for this history file
    REAL(r_std)                  :: dt_stomate_loc     !!  for test : time step of slow processes and STOMATE
    REAL(r_std),DIMENSION(nvm)   :: hist_PFTaxis       !!- GK An axis we need for the history files
!
    REAL(r_std),DIMENSION(10)  :: hist_pool_10axis     !! Deforestation axis
    REAL(r_std),DIMENSION(100)  :: hist_pool_100axis     !! Deforestation axis
    REAL(r_std),DIMENSION(11)  :: hist_pool_11axis     !! Deforestation axis
    REAL(r_std),DIMENSION(101)  :: hist_pool_101axis     !! Deforestation axis
    !
    ! IPCC history file
    REAL(r_std)                  :: hist_days_stom_ipcc     !!- GK time step in days for this history file
    REAL(r_std)                  :: hist_dt_stom_ipcc       !!- GK time step in seconds for this history file
!
    !
    !=====================================================================
    !- 3.0 Setting up the history files
    !=====================================================================
    !- 3.1 SECHIBA
    !=====================================================================
    !Config Key   = ALMA_OUTPUT
    !Config Desc  = Should the output follow the ALMA convention
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = If this logical flag is set to true the model
    !Config         will output all its data according to the ALMA 
    !Config         convention. It is the recommended way to write
    !Config         data out of ORCHIDEE.
    !Config Units = [FLAG]
    CALL getin_p('ALMA_OUTPUT', almaoutput)    
    IF (printlev>=2) WRITE(numout,*) 'ALMA_OUTPUT', almaoutput
    !-
    !Config Key   = OUTPUT_FILE
    !Config Desc  = Name of file in which the output is going to be written
    !Config If    = OK_SECHIBA
    !Config Def   = sechiba_history.nc
    !Config Help  = This file is going to be created by the model
    !Config         and will contain the output from the model.
    !Config         This file is a truly COADS compliant netCDF file.
    !Config         It will be generated by the hist software from
    !Config         the IOIPSL package.
    !Config Units = [FILE]
    !-
    histname='sechiba_history.nc'
    CALL getin_p('OUTPUT_FILE', histname)
    IF (printlev>=2) WRITE(numout,*) 'OUTPUT_FILE', histname
    !-
    !Config Key   = WRITE_STEP
    !Config Desc  = Frequency in seconds for sechiba_history.nc file with IOIPSL
    !Config If    = OK_SECHIBA, NOT XIOS_ORCHIDEE_OK
    !Config Def   = one_day
    !Config Help  = This variables gives the frequency in the output
    !Config         file sechiba_history.nc if using IOIPSL.
    !Config         This variable is not read if XIOS is activated.
    !Config Units = [seconds]
    !-
    dw = one_day
    IF (xios_orchidee_ok) THEN
      dw=0
      IF (printlev>=2) WRITE(numout,*) 'All IOIPSL output are deactivated because this run uses XIOS'
    ELSE
      CALL getin_p('WRITE_STEP', dw)
      IF ( dw == 0 .AND. printlev>=1) WRITE(numout,*) 'sechiba_history file will not be created'
    END IF
    
    veg(1:nvm)   = (/ (REAL(i,r_std),i=1,nvm) /)
    indlai(1:nlai+1) = (/ (REAL(i,r_std),i=1,nlai+1) /)
    soltyp(1:nstm) = (/ (REAL(i,r_std),i=1,nstm) /)
    nobiotyp(1:nnobio) = (/ (REAL(i,r_std),i=1,nnobio) /)
    albtyp(1:2) = (/ (REAL(i,r_std),i=1,2) /)
    solay(1:nslm) = (/ (REAL(i,r_std),i=1,nslm) /)
    snowlev =  (/ (REAL(i,r_std),i=1,nsnow) /)

    ! Get the vertical soil levels for the thermal scheme
    IF (hydrol_cwrr) THEN
       sol(1:ngrnd) = znt(:)
    ELSE
       sol(1:ngrnd) = thermosoilc_levels()
    END IF

    !
    !- We need to flux averaging operation as when the data is written
    !- from within SECHIBA a scatter is needed. In the driver on the other
    !- hand the data is 2D and can be written is it is.
    !-
    WRITE(flux_op,'("ave(scatter(X*",F8.1,"))")') one_day/dt
    ! WRITE(flux_op,'("(ave(scatter(X))*",F8.1,")")') one_day/dt
!    WRITE(flux_sc,'("ave(X*",F8.1,")")') one_day/dt
!    WRITE(flux_insec,'("ave(X*",F8.6,")")') un/dt
!    WRITE(flux_insec,'("ave(X*",F12.10,")")') un/dt
    WRITE(flux_scinsec,'("ave(scatter(X*",F12.10,"))")') un/dt
    IF (printlev>=2) WRITE(numout,*) 'flux_op=',flux_op,' one_day/dt=', one_day/dt, ' dt=',dt,' dw=', dw
    !-
    !Config Key   = SECHIBA_HISTLEVEL
    !Config Desc  = SECHIBA history output level (0..10)
    !Config If    = OK_SECHIBA and HF
    !Config Def   = 5
    !Config Help  = Chooses the list of variables in the history file. 
    !Config         Values between 0: nothing is written; 10: everything is 
    !Config         written are available More details can be found on the web under documentation.
    !Config Units = [-]
    !-
    hist_level = 5
    CALL getin_p('SECHIBA_HISTLEVEL', hist_level)
    !-
    IF (printlev>=2) WRITE(numout,*) 'SECHIBA history level: ',hist_level
    IF ( (hist_level > max_hist_level).OR.(hist_level < 0) ) THEN
       STOP 'This history level is not allowed'
    ENDIF
    !-
    !- define operations as a function of history level.
    !- Above hist_level, operation='never'
    !-
    ave(1:max_hist_level) = 'ave(scatter(X))'
    IF (hist_level < max_hist_level) THEN
       ave(hist_level+1:max_hist_level) = 'never'
    ENDIF
    sumscatter(1:max_hist_level) = 't_sum(scatter(X))'
    IF (hist_level < max_hist_level) THEN
       sumscatter(hist_level+1:max_hist_level) = 'never'
    ENDIF

    avecels(1:max_hist_level) = 'ave(cels(scatter(X)))'
    IF (hist_level < max_hist_level) THEN
       avecels(hist_level+1:max_hist_level) = 'never'
    ENDIF

    avescatter(1:max_hist_level) = 'ave(scatter(X))'
    IF (hist_level < max_hist_level) THEN
       avescatter(hist_level+1:max_hist_level) = 'never'
    ENDIF
    tmincels(1:max_hist_level) = 't_min(cels(scatter(X)))'
    IF (hist_level < max_hist_level) THEN
       tmincels(hist_level+1:max_hist_level) = 'never'
    ENDIF
    tmaxcels(1:max_hist_level) = 't_max(cels(scatter(X)))'
    IF (hist_level < max_hist_level) THEN
       tmaxcels(hist_level+1:max_hist_level) = 'never'
    ENDIF

    fluxop(1:max_hist_level) = flux_op
    IF (hist_level < max_hist_level) THEN
       fluxop(hist_level+1:max_hist_level) = 'never'
    ENDIF

    fluxop_scinsec(1:max_hist_level) = flux_scinsec
    IF (hist_level < max_hist_level) THEN
       fluxop_scinsec(hist_level+1:max_hist_level) = 'never'
    ENDIF
    once(1:max_hist_level) = 'once(scatter(X))'
    IF (hist_level < max_hist_level) THEN
       once(hist_level+1:max_hist_level) = 'never'
    ENDIF


    !- Initialize sechiba_history output file
    !- 
    IF ( dw == 0 ) THEN
       ! sechiba_history file will not be created.
       hist_id = -1

    ELSE
       ! sechiba_history file will be created

       ! If running in parallel (mpi_size>1), test if there are at least 2 latitude bands(jj_nb) for current MPI process. 
       ! The model can work with 1 latitude band but the rebuild fails. Therefor exit if this is the cas. 
       IF ( jj_nb < 2 .AND. mpi_size > 1) THEN
          CALL ipslerr_p(3,"ioipslctrl_history","The current MPI process has jj_nb=1 (1 band of latitude) but", &
               "the IOIPSL rebuild tool can not work if jj_nb is less than 2 per MPI process.", &
               "Change to a lower number of MPI processors or make the region bigger in latitudes.")
       END IF

       !- Calculation necessary for initialization of sechiba_history file
       !- Check if we have by any change a rectilinear grid. This would allow us to 
       !- simplify the output files.
    IF (is_omp_root) THEN
       !
       IF ( GridType == "RegLonLat" ) THEN
          ALLOCATE(lon_rect(iim),stat=ier)
          IF (ier .NE. 0) THEN
             WRITE (numout,*) ' error in lon_rect allocation. We stop. We need iim words = ',iim
             STOP 'intersurf_history'
          ENDIF
          ALLOCATE(lat_rect(jjm),stat=ier)
          IF (ier .NE. 0) THEN
             WRITE (numout,*) ' error in lat_rect allocation. We stop. We need jjm words = ',jjm
             STOP 'intersurf_history'
          ENDIF
          lon_rect(:) = lon(:,1)
          lat_rect(:) = lat(1,:)
       ENDIF
       !-
       !-
       !-
       ! Initialize sechiba_history file
       IF ( .NOT. almaoutput ) THEN
          !- 
          IF ( GridType == "RegLonLat" ) THEN
#ifdef CPP_PARA
             CALL histbeg(histname, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id,orch_domain_id)
#else
             CALL histbeg(histname, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id)
#endif
             IF (printlev >= 2) WRITE(numout,*)  'HISTBEG --->',istp_old,date0,dt,dw,hist_id
          ELSE
#ifdef CPP_PARA
             CALL histbeg(histname, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id,domain_id=orch_domain_id)
#else
             CALL histbeg(histname, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id)
#endif
          ENDIF
          !-
          CALL histvert(hist_id, 'veget', 'Vegetation types', '1', &
               &    nvm,   veg, vegax_id)
          CALL histvert(hist_id, 'laiax', 'Nb LAI', 'm', &
               &   nlai+1,indlai, laiax_id)
          CALL histvert(hist_id, 'solth', 'Soil levels',      'm', &
               &    ngrnd, sol, solax_id)
          CALL histvert(hist_id, 'soiltyp', 'Soil types',      '1', &
               &    nstm, soltyp, soltax_id)
          CALL histvert(hist_id, 'nobio', 'Other surface types',      '1', &
               &    nnobio, nobiotyp, nobioax_id)
          IF (  hydrol_cwrr ) THEN
             CALL histvert(hist_id, 'solay', 'Hydrol soil levels',      'm', &
                  &    nslm, diaglev(1:nslm), solayax_id)
             CALL histvert(hist_id, 'soildiag', 'Diagnostic soil levels', 'm', &
                  &    nslm, diaglev(1:nslm), soildiagax_id)
          ENDIF

          CALL histvert(hist_id, 'snowlev', 'Snow levels',      'm', &
               &    nsnow, snowlev, snowax_id)

          !-
          !- SECHIBA_HISTLEVEL = 1
          !-
          CALL histdef(hist_id, 'evap', 'Evaporation', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(1), dt,dw)
          CALL histdef(hist_id, 'coastalflow', 'Diffuse coastal flow', 'm^3/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'riverflow', 'River flow to the oceans', 'm^3/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw) 
          CALL histdef(hist_id, 'temp_sol', 'Surface Temperature', 'C', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avecels(1), dt,dw)
          CALL histdef(hist_id, 'rain', 'Rainfall', 'mm/d',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(1), dt,dw)
          CALL histdef(hist_id, 'snowf', 'Snowfall', 'mm/d',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(1), dt,dw)
          CALL histdef(hist_id, 'netrad', 'Net radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'lai', 'Leaf Area Index', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(1), dt,dw)
          !
          IF (  hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'reinf_slope', 'Slope index for each grid box', '1', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(1),  dt,dw)
             CALL histdef(hist_id, 'soilindex', 'Soil index', '1', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, once(1),  dt,dw)
          ENDIF
          !
          IF ( river_routing ) THEN
             CALL histdef(hist_id, 'basinmap', 'Aproximate map of the river basins', ' ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw) 
             CALL histdef(hist_id, 'nbrivers', 'Number or rivers in the outflow grid box', ' ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw)  
          ENDIF
          !-
          !- SECHIBA_HISTLEVEL = 2
          !-
          CALL histdef(hist_id, 'subli', 'Sublimation', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
          CALL histdef(hist_id, 'evapnu', 'Bare soil evaporation', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
          CALL histdef(hist_id, 'runoff', 'Surface runoff', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
          CALL histdef(hist_id, 'drainage', 'Deep drainage', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
          IF ( river_routing ) THEN
             CALL histdef(hist_id, 'riversret', 'Return from endorheic rivers', 'mm/d', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
             CALL histdef(hist_id, 'hydrographs', 'Hydrographs of gridbox outflow', 'm^3/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(2), dt,dw)
          ENDIF
          IF ( hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'evapnu_soil', 'Bare soil evap for soil type', 'mm/d', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, fluxop(2), dt,dw)
             CALL histdef(hist_id, 'drainage_soil', 'Drainage for soil type', 'mm/d', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, fluxop(2), dt,dw)
             CALL histdef(hist_id, 'transpir_soil', 'Transpir for soil type', 'mm/d', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, fluxop(2), dt,dw)
             CALL histdef(hist_id, 'runoff_soil', 'Runoff for soil type', 'mm/d', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, fluxop(2), dt,dw)
          ENDIF
          !
          CALL histdef(hist_id, 'tair', 'Air Temperature', 'K',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'qair', 'Air humidity', 'g/g',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'q2m', '2m Air humidity', 'g/g',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 't2m', '2m Air Temperature', 'K',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'alb_vis', 'Albedo visible', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'alb_nir', 'Albedo near infrared', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'soilalb_vis', 'Soil Albedo visible', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'soilalb_nir', 'Soil Albedo near infrared', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'vegalb_vis', 'Vegetation Albedo visible', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'vegalb_nir', 'Vegetation Albedo near infrared', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'z0m', 'Surface roughness for momentum', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'z0h', 'Surface roughness for heat', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'roughheight', 'Effective roughness height', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'transpir', 'Transpiration', 'mm/d', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop(2), dt,dw)
          CALL histdef(hist_id, 'inter', 'Interception loss', 'mm/d', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop(2), dt,dw)
          !-
          !- SECHIBA_HISTLEVEL = 3
          !-
          CALL histdef(hist_id, 'tsol_max', 'Maximum Surface Temperature',&
               & 'C', iim,jjm, hori_id, 1,1,1, -99, 32, tmaxcels(3), dt,dw)
          CALL histdef(hist_id, 'tsol_min', 'Minimum Surface Temperature',&
               & 'C', iim,jjm, hori_id, 1,1,1, -99, 32, tmincels(3), dt,dw)
          CALL histdef(hist_id, 'fluxsens', 'Sensible Heat Flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(3), dt,dw)
          CALL histdef(hist_id, 'fluxlat', 'Latent Heat Flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(3), dt,dw)
          CALL histdef(hist_id, 'snow', 'Snow mass', 'kg/m^2', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'snowage', 'Snow age', '?', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'snownobio', 'Snow on other surfaces', 'kg/m^2', &
               & iim,jjm, hori_id, nnobio,1, nnobio, nobioax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'snownobioage', 'Snow age on other surfaces', 'd', &
               & iim,jjm, hori_id, nnobio,1, nnobio, nobioax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'vegetfrac', 'Fraction of vegetation', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'maxvegetfrac', 'Maximum fraction of vegetation', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'tot_bare_soil', "Total Bare Soil Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'nobiofrac', 'Fraction of other surface types', '1', &
               & iim,jjm, hori_id, nnobio, 1, nnobio, nobioax_id, 32, avescatter(3), dt,dw)
          IF ( do_floodplains ) THEN
             CALL histdef(hist_id, 'flood_frac', 'Flooded fraction', '1', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(3), dt,dw)
             CALL histdef(hist_id, 'reinfiltration', 'Reinfiltration from floodplains', 'mm/d', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(3), dt,dw)
          ENDIF
          IF ( hydrol_cwrr ) THEN
             DO jst=1,nstm
             
                ! var_name= "mc_1" ... "mc_3"
                WRITE (var_name,"('moistc_',i1)") jst
                CALL histdef(hist_id, var_name, 'Soil Moisture profile for soil type', 'm3/m3', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(3),  dt,dw)
                
                ! var_name= "vegetsoil_1" ... "vegetsoil_3"
                WRITE (var_name,"('vegetsoil_',i1)") jst
                CALL histdef(hist_id, var_name, 'Fraction of vegetation on soil types', '%', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3),  dt,dw)
                
                ! var_name= "kfact_root_1" ... "kfact_root_3"
                WRITE (var_name,"('kfactroot_',i1)") jst
                CALL histdef(hist_id, var_name, 'Root fraction profile for soil type', '%', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(3), dt,dw)
                
             ENDDO

             IF (ok_freeze_cwrr) THEN
                CALL histdef(hist_id, 'profil_froz_hydro', 'Frozen fraction for each hydrological soil layer', '-', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1),  dt,dw)
                CALL histdef(hist_id, 'temp_hydro', 'Soil temperature interpolated on hydrological layers', 'K', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1), dt,dw)
             END IF

             CALL histdef(hist_id, 'kk_moy', 'Mean hydrological conductivity', 'mm/d', &
                  & iim,jjm,hori_id, nslm,1,nslm, solayax_id, 32, avescatter(1),  dt,dw)
             DO jst=1,nstm
                WRITE (var_name,"('profil_froz_hydro_',i1)") jst
                CALL histdef(hist_id, trim(var_name), 'Frozen fraction for each hydrological soil layer and soiltile', '-', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1),  dt,dw)
             ENDDO

             CALL histdef(hist_id, 'ptn_beg', 'Soil temperature from previous time step', 'K', &
                  & iim,jjm,hori_id, ngrnd,1,ngrnd, solax_id, 32, avescatter(1),  dt,dw)

             IF ( ok_freeze_thermix ) THEN
                CALL histdef(hist_id, 'pcappa_supp', 'Additional heat capacity due to soil freezing for each soil layer', 'J/K', &
                     iim,jjm,hori_id, ngrnd,1,ngrnd, solax_id, 32, avescatter(1),  dt,dw)
             END IF

             CALL histdef(hist_id, 'shum_ngrnd_perma', 'Saturation degree on the thermal axes', '-', &
                  & iim,jjm,hori_id,ngrnd,1,ngrnd, solax_id, 32, avescatter(1),  dt,dw)
             CALL histdef(hist_id, 'shumdiag_perma', 'Saturation degree of the soil', '-', &
                  & iim,jjm,hori_id,nslm,1,nslm, soildiagax_id, 32, avescatter(1),  dt,dw)
          ENDIF
          !
          CALL histdef(hist_id, 'frac_bare', 'Bare soil fraction for each tile', '-', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'soiltile', 'Fraction of soil tiles', '%', &
               & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, once(3),  dt,dw)
          !-
          !- SECHIBA_HISTLEVEL = 4
          !-
          IF ( .NOT. hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'dss', 'Up-reservoir Height', 'm',  &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'gqsb', 'Upper Soil Moisture', 'Kg/m^2',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'bqsb', 'Lower Soil Moisture', 'Kg/m^2',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
          ELSE
             CALL histdef(hist_id, 'humtot', 'Total Soil Moisture', 'Kg/m^2', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'humtot_soil', 'Soil Moisture for soil type', 'Kg/m^2', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'njsc', 'Soil class used for hydrology', '-', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, once(4), dt,dw)
          ENDIF
          CALL histdef(hist_id, 'qsintveg', 'Water on canopy', 'Kg/m^2', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(4), dt,dw)
          CALL histdef(hist_id, 'rstruct', 'Structural resistance', 's/m', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(4), dt,dw)
          IF ( ok_co2 ) THEN
             CALL histdef(hist_id, 'gpp', 'Net assimilation of carbon by the vegetation', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
          ENDIF
          IF ( ok_stomate ) THEN
             CALL histdef(hist_id, 'nee', 'Net Ecosystem Exchange', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
             CALL histdef(hist_id, 'maint_resp', 'Maintenance respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
             CALL histdef(hist_id, 'hetero_resp', 'Heterotrophic respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
             CALL histdef(hist_id, 'growth_resp', 'Growth respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
             CALL histdef(hist_id, 'npp', 'Net Primary Production', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt, dw)
          ENDIF
          CALL histdef(hist_id, 'precisol', 'Throughfall', 'mm/d',  &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop(4), dt,dw)
          CALL histdef(hist_id, 'drysoil_frac', 'Fraction of visibly dry soil', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(4), dt,dw)
          CALL histdef(hist_id, 'evapot', 'Potential evaporation', 'mm/d',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(4), dt,dw)
          CALL histdef(hist_id, 'evapot_corr', 'Potential evaporation', 'mm/d',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(4), dt,dw)
          !-
          !- SECHIBA_HISTLEVEL = 5
          !-
          CALL histdef(hist_id, 'swnet', 'Net solar radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(5), dt,dw)
          CALL histdef(hist_id, 'swdown', 'Incident solar radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(5), dt,dw)
          CALL histdef(hist_id, 'lwdown', 'Absorbed downward longwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(5), dt,dw)
          CALL histdef(hist_id, 'lwnet', 'Net surface longwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(5), dt,dw)
          !-
          !- SECHIBA_HISTLEVEL = 6
          !-
          CALL histdef(hist_id, 'ptn', 'Deep ground temperature', 'K', &
               & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32, avescatter(6),  dt,dw)
          CALL histdef(hist_id, 'snowmelt', 'snow melt', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(6), dt,dw)
          CALL histdef(hist_id, 'frac_snow_veg', 'snow fraction on vegeted area','-', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(6), dt,dw)
          CALL histdef(hist_id, 'frac_snow_nobio', 'snow fraction on non-vegeted area', '-', &
               & iim,jjm, hori_id, nnobio, 1,nnobio, nobioax_id, 32, avescatter(6), dt,dw)
          CALL histdef(hist_id, 'pgflux', 'extra energy used for melting top snow layer', '-', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(6), dt,dw)

          IF (ok_explicitsnow) THEN
             CALL histdef(hist_id, 'grndflux', 'ground heat flux', 'W/m2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(6), dt,dw)
             CALL histdef(hist_id, 'sfcfrz', 'surface frozen fraction', '-', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(6), dt,dw)
             CALL histdef(hist_id, 'snowrho', 'Snow density profile', 'kg/m3', & 
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6), dt,dw)
             CALL histdef(hist_id, 'snowtemp','Snow temperature profile','K', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'snowdz','Snow depth profile','m', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'snowliq','Snow liquid content profile','m', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'snowgrain','Snow grain profile','m', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'snowheat','Snow Heat profile','J/m2', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'radsink','Solar Radiation profile','W/m2', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
          END IF
          CALL histdef(hist_id, 'Qg', 'Ground heat flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)

          IF (hydrol_cwrr) THEN
             IF (ok_freeze_thermix) THEN
                CALL histdef(hist_id, 'profil_froz', 'Frozen fraction of the soil', '-', &
                     & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32, avescatter(6),  dt,dw)
             END IF
             CALL histdef(hist_id, 'pkappa', 'Soil thermal conductivity', 'W/m/K', &
                  & iim,jjm,hori_id, ngrnd, 1, ngrnd, solax_id, 32, avescatter(1),  dt,dw)
             CALL histdef(hist_id, 'pcapa', 'Apparent heat capacity', 'J/m3/K', &
                  & iim,jjm,hori_id, ngrnd, 1, ngrnd, solax_id, 32, avescatter(1),  dt,dw)
          END IF

          !-
          !- SECHIBA_HISTLEVEL = 7
          !-
          IF ( river_routing ) THEN
             CALL histdef(hist_id, 'fastr', 'Fast flow reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(7), dt,dw)
             CALL histdef(hist_id, 'slowr', 'Slow flow reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(7), dt,dw)
             CALL histdef(hist_id, 'streamr', 'Stream flow reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(7), dt,dw)
             CALL histdef(hist_id, 'lakevol', 'Volume in lake reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(7), dt,dw)
             
             !-
             !- SECHIBA_HISTLEVEL = 8
             !-
             CALL histdef(hist_id, 'pondr', 'Ponds reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
             CALL histdef(hist_id, 'swampmap', 'Map of swamps', 'm^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(8), dt,dw)
             !
             IF ( do_irrigation ) THEN
                CALL histdef(hist_id, 'irrigation', 'Net irrigation', 'mm/d', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(8), dt,dw)
                CALL histdef(hist_id, 'netirrig', 'Net irrigation requirement', 'mm/d', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(8), dt,dw)
                CALL histdef(hist_id, 'irrigmap', 'Map of irrigated surfaces', 'm^2', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, once(8), dt,dw)
             ENDIF
             IF ( do_floodplains ) THEN
                CALL histdef(hist_id, 'floodmap', 'Map of floodplains', 'm^2', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, once(8), dt,dw)
                CALL histdef(hist_id, 'floodh', 'Floodplains height', 'mm', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
                CALL histdef(hist_id, 'floodr', 'Floodplains reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
                CALL histdef(hist_id, 'floodout', 'Flow out of floodplains', 'mm/d', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(8), dt,dw)
                CALL histdef(hist_id, 'evapflo', 'Floodplains evaporation', 'mm/d', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(8), dt,dw)
             ENDIF
             !
          ENDIF

          IF ( hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'k_litt', 'Litter cond', 'mm/d', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          ENDIF
          CALL histdef(hist_id, 'beta', 'Beta Function', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'raero', 'Aerodynamic resistance', 's/m',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'cdrag', 'Drag coefficient for LE and SH', '?',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'Wind', 'Wind speed', 'm/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'qsatt' , 'Surface saturated humidity', 'g/g', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'vbeta1', 'Beta for sublimation', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'vbeta4', 'Beta for bare soil', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'vbeta5', 'Beta for floodplains', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          IF ( ok_co2 ) THEN
             CALL histdef(hist_id, 'gsmean', 'mean stomatal conductance', 'mol/m2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(8), dt,dw)
          ENDIF
          CALL histdef(hist_id, 'humrel', 'Soil moisture stress', '-',  &
               & iim,jjm, hori_id, nvm,1,nvm, vegax_id, 32, avescatter(8), dt,dw)
          !-
          !- SECHIBA_HISTLEVEL = 9
          !-
          !-
          !- SECHIBA_HISTLEVEL = 10
          !-
          IF ( ok_co2 ) THEN
             CALL histdef(hist_id, 'cimean', 'Stomatal CO2 concentation', 'mmole/m2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'cim', 'cim', 'ppm', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'leafci', 'leaf ci', 'ppm', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'gs', 'gs', 'mol m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'assimi', 'assimi', 'mol m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Rd', 'Rd', 'mol m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Cc', 'Cc', 'ppm', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'limitphoto', 'limitphoto', '-', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Vc', 'Vc', 'mol m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Vj', 'Vj', 'mol e- m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'gm', 'gm', 'mol m-2 s-1', &
                  & iim, jjm, hori_id, 1, 1, 1, -99, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'gammastar', 'gammastar', 'ppm', &
                  & iim, jjm, hori_id, 1, 1, 1, -99, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Kmo', 'Kmo', 'ppm', &
                  & iim, jjm, hori_id, 1, 1, 1, -99, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Kmc', 'Kmc', 'ppm', &
                  & iim, jjm, hori_id, 1, 1, 1, -99, 32, avescatter(10), dt,dw)
          ENDIF
          CALL histdef(hist_id, 'vbeta3', 'Beta for Transpiration', 'mm/d', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
          CALL histdef(hist_id, 'rveget', 'Canopy resistance', 's/m', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
          IF ( .NOT. hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'rsol', 'Soil resistance', 's/m',  &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(10), dt,dw)
          ENDIF
          CALL histdef(hist_id,'vbeta2','Beta for Interception loss','mm/d', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
          CALL histdef(hist_id, 'qsintmax', 'Maximum Interception capacity', 'Kg/m^2', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)

          !- SECHIBA_HISTLEVEL = 11
          !-

          CALL histdef(hist_id, 'mrsos', "Moisture in Upper 0.1 m of Soil Column", "kg m-2", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'mrso', "Total Soil Moisture Content", "kg m-2", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'mrros', "Surface Runoff", "kg m-2 s-1", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)
          
          CALL histdef(hist_id, 'mrro', "Total Runoff", "kg m-2 s-1", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)

          CALL histdef(hist_id, 'prveg', "Precipitation onto Canopy", "kg m-2 s-1", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)


          CALL histdef(hist_id, 'evspsblveg', "Evaporation from Canopy", "kg m-2 s-1", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)
          
          CALL histdef(hist_id, 'evspsblsoi', "Water Evaporation from Soil", "kg m-2 s-1", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)
          
          CALL histdef(hist_id, 'tran', "Transpiration", "kg m-2 s-1", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)
          
          CALL histdef(hist_id, 'treeFrac', "Tree Cover Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'grassFrac', "Natural Grass Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'cropFrac', "Crop Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'baresoilFrac', "Bare Soil Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'residualFrac', &
               & "Fraction of Grid Cell that is Land but Neither Vegetation-Covered nor Bare Soil", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          IF ( ok_bvoc ) THEN
             CALL histdef(hist_id, 'PAR', 'PAR', 'umol phot/m^2/s',  &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
             IF ( ok_radcanopy ) THEN
                CALL histdef(hist_id, 'PARsun', 'Sunlit Leaf PAR', 'umol phot/m^2/s', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'PARsh', 'Shaded Leaf Area PAR', 'umol phot/m^2/s', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'laisun', 'Sunlit Leaf Area Index', '1', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'laish', 'Shaded Leaf Area Index', '1', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'Fdf', 'Fdf', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                IF ( ok_multilayer ) then 
                   CALL histdef(hist_id, 'PARsuntab', 'Sunlit Leaf PAR', 'umol phot/m^2/s', &
                        & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(11), dt,dw)
                   CALL histdef(hist_id, 'PARshtab', 'Shaded Leaf Area PAR', 'umol phot/m^2/s', &
                        & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(11), dt,dw)
                ENDIF
                CALL histdef(hist_id, 'coszang', 'coszang', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'PARdf', 'PARdf', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'PARdr', 'PARdr', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'Trans', 'Trans', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
             END IF
             
             CALL histdef(hist_id, 'flx_fertil_no', 'flx_fertil_no', 'ngN/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'CRF', 'CRF', '1', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_co2_bbg_year', 'flx_co2_bbg_year', 'kgC/m^2/yr ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(11), dt,dw)  
             CALL histdef(hist_id, 'N_qt_WRICE_year', 'N_qt_WRICE_year', 'kgN/yr ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(11), dt,dw)  
             CALL histdef(hist_id, 'N_qt_OTHER_year', 'N_qt_OTHER_year', 'kgN/yr ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(11), dt,dw)  
             CALL histdef(hist_id, 'flx_iso', 'flx_iso', 'kgC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_mono', 'flx_mono', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_apinen', 'flx_apinen', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_bpinen', 'flx_bpinen', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_limonen', 'flx_limonen', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_myrcen', 'flx_myrcen', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_sabinen', 'flx_sabinen', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_camphen', 'flx_camphen', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_3caren', 'flx_3caren', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_tbocimen', 'flx_tbocimen', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_othermono', 'flx_othermono', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_sesquiter', 'flx_sesquiter', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_ORVOC', 'flx_ORVOC', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_MBO', 'flx_MBO', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_methanol', 'flx_methanol', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_acetone', 'flx_acetone', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_acetal', 'flx_acetal', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_formal', 'flx_formal', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_acetic', 'flx_acetic', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_formic', 'flx_formic', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_no_soil', 'flx_no_soil', 'ngN/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_no', 'flx_no', 'ngN/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'fco2', 'fco2', '-', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
          ENDIF

       ELSE 
          !-
          !- This is the ALMA convention output now
          !-
          !- 
          IF ( GridType == "RegLonLat" ) THEN
#ifdef CPP_PARA
             CALL histbeg(histname, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id,orch_domain_id)
#else
             CALL histbeg(histname, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id)
#endif
          ELSE
#ifdef CPP_PARA
             CALL histbeg(histname, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id,domain_id=orch_domain_id)
#else
             CALL histbeg(histname, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id)
#endif
          ENDIF
          !-
          CALL histvert(hist_id, 'veget', 'Vegetation types', '1', &
               &    nvm,   veg, vegax_id)
          CALL histvert(hist_id, 'laiax', 'Nb LAI', 'm', &
               &   nlai+1,indlai, laiax_id)
          CALL histvert(hist_id, 'solth', 'Soil levels',      'm', &
               &    ngrnd, sol, solax_id)
          CALL histvert(hist_id, 'soiltyp', 'Soil types',      '1', &
               &    nstm, soltyp, soltax_id)
          CALL histvert(hist_id, 'nobio', 'Other surface types',      '1', &
               &    nnobio, nobiotyp, nobioax_id)
          IF (  hydrol_cwrr ) THEN
             CALL histvert(hist_id, 'solay', 'Hydrol soil levels',      'm', &
                  &    nslm, diaglev(1:nslm), solayax_id)
          ENDIF
          !-
          !-  Vegetation
          !-
          CALL histdef(hist_id, 'vegetfrac', 'Fraction of vegetation', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'maxvegetfrac', 'Maximum fraction of vegetation', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'nobiofrac', 'Fraction of other surface types', '1', &
               & iim,jjm, hori_id, nnobio, 1, nnobio, nobioax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'lai', 'Leaf Area Index', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(1), dt,dw)
          !-
          !- Forcing variables
          !-
          CALL histdef(hist_id, 'SinAng', 'Net shortwave radiation', '-',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(9), dt,dw)
          CALL histdef(hist_id, 'LWdown', 'Downward longwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(9), dt,dw)
          CALL histdef(hist_id, 'SWdown', 'Downward shortwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(9), dt,dw)
          CALL histdef(hist_id, 'Tair', 'Near surface air temperature at forcing level', 'K',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'Qair', 'Near surface specific humidity at forcing level', 'g/g',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'SurfP', 'Surface Pressure', 'hPa',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(9), dt,dw)
          CALL histdef(hist_id, 'Windu', 'Eastward wind', 'm/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(9), dt,dw)
          CALL histdef(hist_id, 'Windv', 'Northward wind', 'm/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(9), dt,dw)
          !- 
          !-  General energy balance
          !-
          CALL histdef(hist_id, 'SWnet', 'Net shortwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'LWnet', 'Net longwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'Qh', 'Sensible heat flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'Qle', 'Latent heat flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'Qg', 'Ground heat flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'Qf', 'Energy of fusion', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'Qv', 'Energy of sublimation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'DelSurfHeat', 'Change in surface layer heat', 'J/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, sumscatter(1), dt,dw)
          CALL histdef(hist_id, 'DelColdCont', 'Change in snow surface layer cold content', 'J/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, sumscatter(1), dt,dw)
          !-
          !- General water balance
          !-
          CALL histdef(hist_id, 'Snowf', 'Snowfall rate', 'kg/m^2/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Rainf', 'Rainfall rate', 'kg/m^2/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Evap', 'Total Evapotranspiration', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Qs', 'Surface runoff', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Qsb', 'Sub-surface runoff', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Qrec', 'Recharge', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Qsm', 'Snowmelt', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'DelSoilMoist', 'Change in soil moisture', 'kg/m^2',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, sumscatter(1), dt,dw)
          CALL histdef(hist_id, 'DelSurfStor', 'Change in Surface Water Storage','kg/m^2',&
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, sumscatter(1), dt,dw)
          CALL histdef(hist_id, 'DelIntercept', 'Change in interception storage', 'kg/m^2',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, sumscatter(1), dt,dw)
          CALL histdef(hist_id, 'DelSWE', 'Change in Snow Water Equivalent', 'kg/m^2',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, sumscatter(1), dt,dw)
          IF ( do_irrigation ) THEN
             CALL histdef(hist_id, 'Qirrig', 'Irrigation', 'kg/m^2/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'Qirrig_req', 'Irrigation requirement', 'kg/m^2/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          ENDIF
          !-
          !- Surface state
          !-
          CALL histdef(hist_id, 'AvgSurfT', 'Average surface temperature', 'K', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'PotSurfT', 'Potential (Unstressed) surface temperature', 'K', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'RadT', 'Surface radiative temperature', 'K', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'Albedo', 'Albedo', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SurfStor', 'Surface Water Storage','kg/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SWE', 'Snow Water Equivalent', 'kg/m^2', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'InterceptVeg', 'Intercepted Water on Canopy', 'Kg/m^2', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(1), dt,dw)
          !!-
          !-  Sub-surface state
          !- 
          IF ( .NOT. hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'SoilMoist', '3D average layer soil moisture', 'kg/m^2',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(1), dt,dw)
          ELSE
             CALL histdef(hist_id, 'SoilMoist', '3D average layer soil moisture', 'kg/m^2',  &
                  & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1), dt,dw)

             IF (ok_freeze_cwrr) THEN
                CALL histdef(hist_id, 'profil_froz_hydro', 'Frozen fraction for each hydrological soil layer', '-', &
                     & iim,jjm, hori_id, nslm, 1, nslm,solayax_id, 32, avescatter(1),  dt,dw)
                DO jst=1,nstm
                   WRITE (var_name,"('profil_froz_hydro_',i1)") jst
                   CALL histdef(hist_id, trim(var_name), 'Frozen fraction for each hydrological soil layer and soiltile', '-', &
                        & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1),  dt,dw)
                ENDDO

                CALL histdef(hist_id, 'temp_hydro', 'Soil temperature interpolated on hydrological layers', 'K', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1), dt,dw)
                CALL histdef(hist_id, 'kk_moy', 'Mean hydrological conductivity', 'mm/d', &
                     & iim,jjm,hori_id, nslm,1,nslm, solayax_id, 32, avescatter(1),  dt,dw)
             ENDIF
          END IF
          CALL histdef(hist_id, 'SoilWet', 'Total soil wetness', '-',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SoilTemp', 'Soil temperature profile', 'K', &
               & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32, avescatter(1),  dt,dw)
          !- 
          !-  Evaporation components
          !-
          CALL histdef(hist_id, 'PotEvap', 'Potential evapotranspiration', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'PotEvapOld', 'Potential evapotranspiration old method', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'ECanop', 'Interception evaporation', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'TVeg', 'Transpiration', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'ESoil', 'Bare soil evaporation', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'EWater', 'Open water evaporation', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'RootMoist','Root zone soil water', 'kg/m^2',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SubSnow','Snow sublimation','kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'ACond', 'Aerodynamic conductance', 'm/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          IF ( do_floodplains ) THEN
             CALL histdef(hist_id, 'Qflood', 'Floodplain Evaporation', 'kg/m^2/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          ENDIF
          !-
          !- Surface turbulence
          !-
          CALL histdef(hist_id, 'Z0m', 'Roughness height for momentum', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'Z0h', 'Roughness height for heat', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'EffectHeight', 'Effective displacement height (h-d)', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          !-
          !-
          !-  Cold Season Processes
          !-
          CALL histdef(hist_id, 'SnowFrac', 'Snow cover fraction', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SAlbedo', 'Snow albedo', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SnowDepth', '3D snow depth', 'm', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          !-
          !- Hydrologic variables
          !-
          IF ( river_routing ) THEN
             CALL histdef(hist_id, 'slowr', 'Slow flow reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
             CALL histdef(hist_id, 'fastr', 'Fast flow reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
             CALL histdef(hist_id, 'streamr', 'Stream flow reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
             CALL histdef(hist_id, 'lakevol', 'Volume in lake reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
             CALL histdef(hist_id, 'pondr', 'Ponds reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
             !-
             !-
             !-
             CALL histdef(hist_id, 'SwampMap', 'Map of swamp areas', 'm^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw)
             CALL histdef(hist_id, 'Dis', 'Simulated River Discharge', 'm^3/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(2), dt,dw)
             CALL histdef(hist_id, 'CoastalFlow', 'Diffuse coastal flow', 'm^3/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(2), dt,dw)
             CALL histdef(hist_id, 'RiverFlow', 'River flow to the oceans', 'm^3/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(2), dt,dw)
             IF ( do_irrigation ) THEN
                CALL histdef(hist_id, 'IrrigationMap', 'Map of irrigated areas', 'm^2', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw)
             ENDIF
             !
             !
             IF ( do_floodplains ) THEN
                CALL histdef(hist_id, 'FloodplainsMap', 'Map of flooded areas', 'm^2', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw)
                CALL histdef(hist_id, 'FloodFrac', 'Floodplain Fraction', '-', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
             ENDIF
          ENDIF
          !-
          CALL histdef(hist_id, 'humrel', 'Soil moisture stress', '-',  &
               & iim,jjm, hori_id, nvm,1,nvm, vegax_id, 32, avescatter(8), dt,dw)
          !-
          !-  The carbon budget
          !-
          IF ( ok_co2 ) THEN
            CALL histdef(hist_id, 'cimean', 'Stomatal CO2 concentation', 'mmole/m2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
            CALL histdef(hist_id, 'GPP', 'Net assimilation of carbon by the vegetation', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
            CALL histdef(hist_id, 'gsmean', 'mean stomatal conductance', 'mol/m2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(8), dt,dw)
             CALL histdef(hist_id, 'cim', 'cim', 'ppm', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'leafci', 'leaf Ci', 'ppm', &
                  & iim,jjm, hori_id,nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'gs', 'gs', 'mol m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'assimi', 'assimi', 'mol m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Rd', 'Rd', 'mol m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Cc', 'Cc', 'ppm', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'limitphoto', 'limitphoto', '-', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Vc', 'Vc', 'mol m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Vj', 'Vj', 'mol e- m-2 s-1', &
                  & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'gm', 'gm', 'mol m-2 s-1', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'gammastar', 'gammastar', 'ppm', &
                  & iim, jjm, hori_id, 1, 1, 1, -99, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Kmo', 'Kmo', 'ppm', &
                  & iim, jjm, hori_id, 1, 1, 1, -99, 32, avescatter(10), dt,dw)
             CALL histdef(hist_id, 'Kmc', 'Kmc', 'ppm', &
                  & iim, jjm, hori_id, 1, 1, 1, -99, 32, avescatter(10), dt,dw)
          ENDIF
          IF ( ok_stomate ) THEN
             CALL histdef(hist_id, 'NEE', 'Net Ecosystem Exchange', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'maint_resp', 'Maintenance respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'hetero_resp', 'Heterotrophic respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'growth_resp', 'Growth respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'npp', 'Net Primary Production', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
          ENDIF
          !
      ENDIF
       !-
       !- Forcing and grid information
       !-
       CALL histdef(hist_id, 'LandPoints', 'Land Points', '1', &
            & iim,jjm, hori_id, 1,1,1, -99, 32, once(10), dt,dw)  
       CALL histdef(hist_id, 'Areas', 'Mesh areas', 'm2', &
            & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt, dw)
       CALL histdef(hist_id, 'Contfrac', 'Continental fraction', '1', &
            & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt, dw)
       !-
       ! Write the names of the pfts in the history files
       global_attribute="PFT_name"
       DO i=1,nvm
          WRITE(global_attribute(9:10),"(I2.2)") i
          CALL histglobal_attr(hist_id, global_attribute, PFT_name(i))
       ENDDO
       !-
       CALL histend(hist_id)
    ENDIF ! IF (is_omp_root)
 
    END IF !IF ( dw == 0 )
    !
    !
    ! Second SECHIBA hist file
    !
    !-
    !Config Key   = SECHIBA_HISTFILE2
    !Config Desc  = Flag to switch on histfile 2 for SECHIBA (hi-frequency ?)
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = This Flag switch on the second SECHIBA writing for hi (or low) 
    !Config         frequency writing. This second output is optional and not written
    !Config         by default.
    !Config Units = [FLAG]
    !-
    ok_histfile2=.FALSE.
    CALL getin_p('SECHIBA_HISTFILE2', ok_histfile2)
    IF (printlev >= 2) WRITE(numout,*) 'SECHIBA_HISTFILE2 ', ok_histfile2
    !
    !-
    !Config Key   = WRITE_STEP2
    !Config Desc  = Frequency in seconds at which to WRITE output
    !Config If    = SECHIBA_HISTFILE2
    !Config Def   = 1800.0
    !Config Help  = This variables gives the frequency the output 2 of
    !Config         the model should be written into the netCDF file.
    !Config         It does not affect the frequency at which the
    !Config         operations such as averaging are done.
    !Config         That is IF the coding of the calls to histdef
    !Config         are correct !
    !Config Units = [seconds]
    !-
    dw2 = 1800.0
    CALL getin_p('WRITE_STEP2', dw2)
    
    ! Deactivate sechiba_histfile2 if the frequency is set to zero
    IF ( dw2 == 0 ) THEN
       ok_histfile2=.FALSE.
       IF (printlev >= 2) WRITE(numout,*) 'WRITE_STEP2 was set to zero and therfore SECHIBA_HISTFILE2 is deactivated.'
    ELSE IF ( hist_id < 0 ) THEN
       ! Deactivate all history files if sechiba_history file is deactivated
       ok_histfile2=.FALSE.
       IF (printlev >= 2) WRITE(numout,*) 'SECHIBA_HISTFILE2 will not be created because sechiba_history file is deactivated.'
    END IF

    hist2_id = -1
    !
    IF (ok_histfile2) THEN
       !-
       !Config Key   = SECHIBA_OUTPUT_FILE2
       !Config Desc  = Name of file in which the output number 2 is going to be written
       !Config If    = SECHIBA_HISTFILE2
       !Config Def   = sechiba_out_2.nc
       !Config Help  = This file is going to be created by the model
       !Config         and will contain the output 2 from the model.
       !Config Units = [FILE]
       !-
       histname2='sechiba_out_2.nc'
       CALL getin_p('SECHIBA_OUTPUT_FILE2', histname2)
       IF (printlev >= 2) WRITE(numout,*) 'SECHIBA_OUTPUT_FILE2 ', histname2
       !-
       !Config Key   = SECHIBA_HISTLEVEL2
       !Config Desc  = SECHIBA history 2 output level (0..10)
       !Config If    = SECHIBA_HISTFILE2
       !Config Def   = 1
       !Config Help  = Chooses the list of variables in the history file. 
       !Config         Values between 0: nothing is written; 10: everything is 
       !Config         written are available More details can be found on the web under documentation.
       !Config         web under documentation.
       !Config         First level contains all ORCHIDEE outputs.
       !Config Units = [-] 
       !-
       hist2_level = 1
       CALL getin_p('SECHIBA_HISTLEVEL2', hist2_level)
       !-
       IF (printlev >= 2) WRITE(numout,*) 'SECHIBA history level 2 : ',hist2_level
       IF ( (hist2_level > max_hist_level).OR.(hist2_level < 0) ) THEN
          STOP 'This history level 2 is not allowed'
       ENDIF
       !
       !-
       !- define operations as a function of history level.
       !- Above hist2_level, operation='never'
       !-
       ave2(1:max_hist_level) = 'ave(scatter(X))'
       IF (hist2_level < max_hist_level) THEN
          ave2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       sumscatter2(1:max_hist_level) = 't_sum(scatter(X))'
       IF (hist2_level < max_hist_level) THEN
          sumscatter2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       avecels2(1:max_hist_level) = 'ave(cels(scatter(X)))'
       IF (hist2_level < max_hist_level) THEN
          avecels2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       avescatter2(1:max_hist_level) = 'ave(scatter(X))'
       IF (hist2_level < max_hist_level) THEN
          avescatter2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       tmincels2(1:max_hist_level) = 't_min(cels(scatter(X)))'
       IF (hist2_level < max_hist_level) THEN
          tmincels2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       tmaxcels2(1:max_hist_level) = 't_max(cels(scatter(X)))'
       IF (hist2_level < max_hist_level) THEN
          tmaxcels2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       fluxop2(1:max_hist_level) = flux_op
       IF (hist2_level < max_hist_level) THEN
          fluxop2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       fluxop_scinsec2(1:max_hist_level) = flux_scinsec
       IF (hist2_level < max_hist_level) THEN
          fluxop_scinsec2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       once2(1:max_hist_level) = 'once(scatter(X))'
       IF (hist2_level < max_hist_level) THEN
          once2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       ! 
       IF (is_omp_root) THEN
          IF ( .NOT. almaoutput ) THEN
             !- 
             IF ( GridType == "RegLonLat" ) THEN
#ifdef CPP_PARA
                CALL histbeg(histname2, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id,orch_domain_id)
#else
                CALL histbeg(histname2, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id)
#endif
                IF (printlev >= 2) WRITE(numout,*)  'HISTBEG2 --->',istp_old,date0,dt,dw2,hist2_id
             ELSE
#ifdef CPP_PARA
                CALL histbeg(histname2, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id,domain_id=orch_domain_id)
#else
                CALL histbeg(histname2, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id)
#endif
             ENDIF
             !-
             CALL histvert(hist2_id, 'veget', 'Vegetation types', '1', &
                  &    nvm,   veg, vegax_id2)
             CALL histvert(hist2_id, 'laiax', 'Nb LAI', '1', &
                  &    nlai+1,   indlai, laiax_id2)
             CALL histvert(hist2_id, 'solth', 'Soil levels',      'm', &
                  &    ngrnd, sol, solax_id2)
             CALL histvert(hist2_id, 'soiltyp', 'Soil types',      '1', &
                  &    nstm, soltyp, soltax_id2)
             CALL histvert(hist2_id, 'nobio', 'Other surface types',      '1', &
                  &    nnobio, nobiotyp, nobioax_id2)
             CALL histvert(hist2_id, 'albtyp', 'Albedo Types',     '1', &
                  &    2, albtyp, albax_id2)
             IF (  hydrol_cwrr ) THEN
                CALL histvert(hist2_id, 'solay', 'Hydrol soil levels',      'm', &
                     &    nslm, solay, solayax_id2)
             ENDIF
             !-
             !- SECHIBA_HISTLEVEL2 = 1
             !-
             CALL histdef(hist2_id, 'ptn', 'Deep ground temperature', 'K', &
                  & iim,jjm, hori_id2, ngrnd, 1, ngrnd, solax_id2, 32, avescatter2(2),  dt, dw2)  

             CALL histdef(hist2_id, 'mrsos', "Moisture in Upper 0.1 m of Soil Column", "kg m-2", &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(2), dt,dw2)            

             CALL histdef(hist2_id, 'mrso', "Total Soil Moisture Content", "kg m-2", &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(11), dt,dw2)
             
             CALL histdef(hist2_id, 'mrros', "Surface Runoff", "kg m-2 s-1", &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, fluxop_scinsec2(11), dt,dw2)

             CALL histdef(hist2_id, 'mrro', "Total Runoff", "kg m-2 s-1", &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, fluxop_scinsec2(2), dt,dw2)      

             !-
             !- SECHIBA_HISTLEVEL2 = 2
             !-
             CALL histdef(hist2_id, 'cdrag', 'Drag coefficient for LE and SH', '?',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt, dw2)
             CALL histdef(hist2_id, 'soilalb_vis', 'Soil Albedo visible', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
             CALL histdef(hist2_id, 'soilalb_nir', 'Soil Albedo near infrared', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
             CALL histdef(hist2_id, 'vegalb_vis', 'Vegetation Albedo visible', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
             CALL histdef(hist2_id, 'vegalb_nir', 'Vegetation Albedo near infrared', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
             CALL histdef(hist2_id, 'z0m', 'Surface roughness for momentum', 'm',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt, dw2)
             CALL histdef(hist2_id, 'z0h', 'Surface roughness for heat', 'm',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt, dw2)
             CALL histdef(hist2_id, 'coastalflow', 'Diffuse coastal flow', 'm^3/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(2), dt, dw2)
             CALL histdef(hist2_id, 'riverflow', 'River flow to the oceans', 'm^3/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(2), dt, dw2) 
             CALL histdef(hist2_id, 'tsol_rad', 'Radiative surface temperature', 'C', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avecels2(2), dt, dw2)
             IF ( do_floodplains ) THEN
                CALL histdef(hist2_id, 'floodout', 'Flow out of floodplains', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(2), dt,dw2)
                CALL histdef(hist2_id, 'vevapflo', 'Floodplains evaporation', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(2), dt, dw2)
                CALL histdef(hist2_id, 'flood_frac', 'Flooded fraction', '1', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
                CALL histdef(hist2_id, 'reinfiltration', 'Reinfiltration from floodplains', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(2), dt,dw2)
             ENDIF
             CALL histdef(hist2_id, 'vevapnu', 'Bare soil evaporation', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(2), dt, dw2)
             CALL histdef(hist2_id, 'temp_sol', 'New Surface Temperature', 'C', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avecels2(2), dt, dw2)
             CALL histdef(hist2_id, 'qsurf', 'Near surface specific humidity', 'g/g',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(2), dt, dw2)
             CALL histdef(hist2_id, 'albedo', 'Albedo', '1', &
                  & iim,jjm, hori_id2, 2,1,2, albax_id2, 32, avescatter2(2), dt, dw2)
             CALL histdef(hist2_id, 'fluxsens', 'Sensible Heat Flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(2), dt, dw2)
             CALL histdef(hist2_id, 'fluxlat', 'Latent Heat Flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(2), dt, dw2)
             CALL histdef(hist2_id, 'emis', 'Surface emissivity', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(2), dt, dw2)
             !-
             !- SECHIBA_HISTLEVEL2 = 3
             !-
             CALL histdef(hist2_id, 'evap', 'Evaporation', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(3), dt, dw2)
             CALL histdef(hist2_id, 'rain', 'Rainfall', 'mm/d',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(3), dt, dw2)
             CALL histdef(hist2_id, 'snowf', 'Snowfall', 'mm/d',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(3), dt, dw2)
             CALL histdef(hist2_id, 'netrad', 'Net radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(3), dt, dw2)
             CALL histdef(hist2_id, 'lai', 'Leaf Area Index', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(3), dt, dw2)
             IF ( river_routing ) THEN
                CALL histdef(hist2_id, 'basinmap', 'Aproximate map of the river basins', ' ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(3), dt, dw2) 
                CALL histdef(hist2_id, 'nbrivers', 'Number or rivers in the outflow grid box', ' ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(3), dt, dw2)  
             ENDIF
             IF (check_waterbal) THEN
                CALL histdef(hist2_id, 'TotWater', 'Total amount of water at end of time step', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(3), dt, dw2)
                CALL histdef(hist2_id, 'TotWaterFlux', 'Total water flux', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(3), dt, dw2)
             ENDIF

             !-
             !- SECHIBA_HISTLEVEL2 = 4
             !-
             CALL histdef(hist2_id, 'subli', 'Sublimation', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(4), dt, dw2)
             CALL histdef(hist2_id, 'runoff', 'Surface runoff', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(4), dt, dw2)
             CALL histdef(hist2_id, 'drainage', 'Deep drainage', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(4), dt, dw2)
             IF ( river_routing ) THEN
                CALL histdef(hist2_id, 'riversret', 'Return from endorheic rivers', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(4), dt, dw2)
                CALL histdef(hist2_id, 'hydrographs', 'Hydrographs of gridbox outflow', 'm^3/s', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(4), dt, dw2)
             ENDIF
             IF ( hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'evapnu_soil', 'Bare soil evap for soil type', 'mm/d', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, fluxop2(4), dt, dw2)
                CALL histdef(hist2_id, 'drainage_soil', 'Drainage for soil type', 'mm/d', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, fluxop2(4), dt, dw2)
                CALL histdef(hist2_id, 'transpir_soil', 'Transpir for soil type', 'mm/d', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, fluxop2(4), dt, dw2)
                CALL histdef(hist2_id, 'runoff_soil', 'Runoff for soil type', 'mm/d', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, fluxop2(4), dt, dw2)
             ENDIF
             !
             CALL histdef(hist2_id, 'tair', 'Air Temperature', 'K',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             CALL histdef(hist2_id, 'qair', 'Air humidity', 'g/g',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             CALL histdef(hist2_id, 'q2m', '2m Air humidity', 'g/g',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             CALL histdef(hist2_id, 't2m', '2m Air Temperature', 'K',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
            CALL histdef(hist2_id, 'alb_vis', 'Albedo visible', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             CALL histdef(hist2_id, 'alb_nir', 'Albedo near infrared', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             CALL histdef(hist2_id, 'roughheight', 'Effective roughness height', 'm',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(4), dt, dw2)
             CALL histdef(hist2_id, 'transpir', 'Transpiration', 'mm/d', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop2(4), dt, dw2)
             CALL histdef(hist2_id, 'inter', 'Interception loss', 'mm/d', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop2(4), dt, dw2)
             !-
             !- SECHIBA_HISTLEVEL2 = 5
             !-
             CALL histdef(hist2_id, 'tsol_max', 'Maximum Surface Temperature',&
                  & 'C', iim,jjm, hori_id2, 1,1,1, -99, 32, tmaxcels2(5), dt, dw2)
             CALL histdef(hist2_id, 'tsol_min', 'Minimum Surface Temperature',&
                  & 'C', iim,jjm, hori_id2, 1,1,1, -99, 32, tmincels2(5), dt, dw2)
             CALL histdef(hist2_id, 'snow', 'Snow mass', 'kg/m^2', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'snowage', 'Snow age', '?', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'snownobio', 'Snow on other surfaces', 'kg/m^2', &
                  & iim,jjm, hori_id2, nnobio,1, nnobio, nobioax_id2, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'snownobioage', 'Snow age on other surfaces', 'd', &
                  & iim,jjm, hori_id2, nnobio,1, nnobio, nobioax_id2, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'vegetfrac', 'Fraction of vegetation', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'maxvegetfrac', 'Maximum fraction of vegetation', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'nobiofrac', 'Fraction of other surface types', '1', &
                  & iim,jjm, hori_id2, nnobio, 1, nnobio, nobioax_id2, 32, avescatter2(5), dt, dw2)
             IF ( hydrol_cwrr ) THEN
                DO jst=1,nstm
                   
                   ! var_name= "mc_1" ... "mc_3"
                   WRITE (var_name,"('moistc_',i1)") jst
                   CALL histdef(hist2_id, var_name, 'Soil Moisture profile for soil type', 'm3/m3', &
                        & iim,jjm, hori_id2, nslm, 1, nslm, solayax_id2, 32, avescatter2(5), dt, dw2)
                   
                   ! var_name= "vegetsoil_1" ... "vegetsoil_3"
                   WRITE (var_name,"('vegetsoil_',i1)") jst
                   CALL histdef(hist2_id, var_name, 'Fraction of vegetation on soil types', '%', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(5), dt, dw2)
                   
                   ! var_name= "kfact_root_1" ... "kfact_root_3"
                   WRITE (var_name,"('kfactroot_',i1)") jst
                   CALL histdef(hist2_id, var_name, 'Root fraction profile for soil type', '%', &
                        & iim,jjm, hori_id2, nslm, 1, nslm, solayax_id2, 32, avescatter2(5), dt,dw2)
                ENDDO

             ENDIF
             !-
             !- SECHIBA_HISTLEVEL2 = 6
             !-
             IF ( .NOT. hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'dss', 'Up-reservoir Height', 'm',  &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id, 32, avescatter2(6), dt,dw2)
                CALL histdef(hist2_id, 'gqsb', 'Upper Soil Moisture', 'Kg/m^2',  &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(6), dt, dw2)
                CALL histdef(hist2_id, 'bqsb', 'Lower Soil Moisture', 'Kg/m^2',  &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(6), dt, dw2)
             ELSE
                CALL histdef(hist2_id, 'humtot', 'Total Soil Moisture', 'Kg/m^2', &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(6), dt, dw2)
                CALL histdef(hist2_id, 'humtot_soil', 'Soil Moisture for soil type', 'Kg/m^2', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, avescatter2(6), dt, dw2)
             ENDIF
             CALL histdef(hist2_id, 'qsintveg', 'Water on canopy', 'Kg/m^2', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(6), dt, dw2)
             CALL histdef(hist2_id, 'rstruct', 'Structural resistance', 's/m', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(6), dt, dw2)
             IF ( ok_co2 ) THEN
                CALL histdef(hist2_id, 'gpp', 'Net assimilation of carbon by the vegetation', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt,dw2)
             ENDIF
             IF ( ok_stomate ) THEN
                CALL histdef(hist2_id, 'nee', 'Net Ecosystem Exchange', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt,dw2)
                CALL histdef(hist2_id, 'maint_resp', 'Maintenance respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt,dw2)
                CALL histdef(hist2_id, 'hetero_resp', 'Heterotrophic respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt,dw2)
                CALL histdef(hist2_id, 'growth_resp', 'Growth respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt, dw2)
                CALL histdef(hist2_id, 'npp', 'Net Primary Production', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt, dw2)
             ENDIF
             CALL histdef(hist2_id, 'precisol', 'Throughfall', 'mm/d',  &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop2(6), dt, dw2)
             CALL histdef(hist2_id, 'drysoil_frac', 'Fraction of visibly dry soil', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(6), dt, dw2)
             CALL histdef(hist2_id, 'evapot', 'Potential evaporation', 'mm/d',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(6), dt, dw2)
             CALL histdef(hist2_id, 'evapot_corr', 'Potential evaporation', 'mm/d',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(6), dt, dw2)
             CALL histdef(hist2_id, 'snowmelt', 'snow melt', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(6), dt, dw2)

             !-
             !- SECHIBA_HISTLEVEL2 = 7
             !-
             CALL histdef(hist2_id, 'swnet', 'Net solar radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(7), dt, dw2)
             CALL histdef(hist2_id, 'swdown', 'Incident solar radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(7), dt, dw2)
             CALL histdef(hist2_id, 'lwdown', 'Absorbed downward longwave radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'lwnet', 'Net surface longwave radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             !-
             !- SECHIBA_HISTLEVEL2 = 8
             !-
             IF ( river_routing ) THEN
                CALL histdef(hist2_id, 'fastr', 'Fast flow reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt, dw2)
                CALL histdef(hist2_id, 'slowr', 'Slow flow reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt, dw2)
                CALL histdef(hist2_id, 'streamr', 'Stream flow reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt, dw2)
                CALL histdef(hist2_id, 'floodr', 'Floodplains reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt,dw2)
                CALL histdef(hist2_id, 'floodh', 'Floodplains height', 'mm', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt,dw2)
                CALL histdef(hist2_id, 'pondr', 'Ponds reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt,dw2)
                CALL histdef(hist2_id, 'lakevol', 'Volume in lake reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt, dw2)
                IF ( do_irrigation ) THEN
                   CALL histdef(hist2_id, 'irrigation', 'Net irrigation', 'mm/d', &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(8), dt, dw2)
                   CALL histdef(hist2_id, 'netirrig', 'Net irrigation requirement', 'mm/d', &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(8), dt, dw2)
                   CALL histdef(hist2_id, 'irrigmap', 'Map of irrigated areas', 'm^2', &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(8), dt, dw2)
                ENDIF
                CALL histdef(hist2_id, 'floodmap', 'Map of floodplains', 'm^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(8), dt,dw2)
                CALL histdef(hist2_id, 'swampmap', 'Map of swamps', 'm^2', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(8), dt,dw2)
             ENDIF
             !-
             !- SECHIBA_HISTLEVEL2 = 9
             !-
             CALL histdef(hist2_id, 'beta', 'Beta Function', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             CALL histdef(hist2_id, 'raero', 'Aerodynamic resistance', 's/m',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             CALL histdef(hist2_id, 'Wind', 'Wind speed', 'm/s',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             CALL histdef(hist2_id, 'qsatt' , 'Surface saturated humidity', 'g/g', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             CALL histdef(hist2_id, 'vbeta1', 'Beta for sublimation', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             CALL histdef(hist2_id, 'vbeta4', 'Beta for bare soil', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             IF ( ok_co2 ) THEN
                CALL histdef(hist2_id, 'gsmean', 'mean stomatal conductance', 'mol/m2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(9), dt, dw2)
             ENDIF
             CALL histdef(hist2_id, 'vbeta5', 'Beta for floodplains', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             IF (  hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'reinf_slope', 'Slope index for each grid box', '1', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(9),  dt,dw2)
                CALL histdef(hist2_id, 'soilindex', 'Soil index', '1', &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, once2(9),  dt,dw2)
             ENDIF
             CALL histdef(hist2_id, 'humrel', 'Soil moisture stress', '-',  &
                  & iim,jjm, hori_id2, nvm,1,nvm, vegax_id2, 32, avescatter2(9), dt, dw2)
             !-
             !- SECHIBA_HISTLEVEL2 = 10
             !-
             IF ( ok_co2 ) THEN
                CALL histdef(hist2_id, 'cimean', 'Stomatal CO2 concentation', 'mmole/m2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             ENDIF
             CALL histdef(hist2_id, 'vbeta3', 'Beta for Transpiration', 'mm/d', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             CALL histdef(hist2_id, 'rveget', 'Canopy resistance', 's/m', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             IF ( .NOT. hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'rsol', 'Soil resistance', 's/m',  &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt, dw2)
             ENDIF
             CALL histdef(hist2_id,'vbeta2','Beta for Interception loss','mm/d', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             CALL histdef(hist2_id, 'qsintmax', 'Maximum Interception capacity', 'Kg/m^2', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             
             IF ( ok_bvoc ) THEN
                CALL histdef(hist2_id, 'PAR', 'PAR', 'umol phot/m^2/s',  &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                IF ( ok_radcanopy ) THEN
                   CALL histdef(hist2_id, 'PARsun', 'Sunlit Leaf PAR', 'umol phot/m^2/s', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'PARsh', 'Shaded Leaf Area PAR', 'umol phot/m^2/s', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'laisun', 'Sunlit Leaf Area Index', '1', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'laish', 'Shaded Leaf Area Index', '1', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'Fdf', 'Fdf', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   IF ( ok_multilayer ) then 
                      CALL histdef(hist2_id, 'PARsuntab', 'Sunlit Leaf PAR', 'umol phot/m^2/s', &
                           & iim,jjm, hori_id2, nlai+1, 1, nlai+1, laiax_id2, 32, avescatter2(10), dt,dw2)
                      CALL histdef(hist2_id, 'PARshtab', 'Shaded Leaf Area PAR', 'umol phot/m^2/s', &
                           & iim,jjm, hori_id2, nlai+1, 1, nlai+1, laiax_id2, 32, avescatter2(10), dt,dw2)
                   ENDIF
                   CALL histdef(hist2_id, 'coszang', 'coszang', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'PARdf', 'PARdf', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'PARdr', 'PARdr', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'Trans', 'Trans', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                END IF
                
                CALL histdef(hist2_id, 'flx_fertil_no', 'flx_fertil_no', 'ngN/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'CRF', 'CRF', '1', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_co2_bbg_year', 'flx_co2_bbg_year', 'kgC/m^2/yr ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(10), dt,dw2)  
                CALL histdef(hist2_id, 'N_qt_WRICE_year', 'N_qt_WRICE_year', 'kgN/yr ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(10), dt,dw2)  
                CALL histdef(hist2_id, 'N_qt_OTHER_year', 'N_qt_OTHER_year', 'kgN/yr ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(10), dt,dw2)  
                CALL histdef(hist2_id, 'flx_iso', 'flx_iso', 'kgC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_mono', 'flx_mono', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_apinen', 'flx_apinen', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_bpinen', 'flx_bpinen', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_limonen', 'flx_limonen', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_myrcen', 'flx_myrcen', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_sabinen', 'flx_sabinen', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_camphen', 'flx_camphen', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_3caren', 'flx_3caren', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_tbocimen', 'flx_tbocimen', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_othermono', 'flx_othermono', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_sesquiter', 'flx_sesquiter', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_ORVOC', 'flx_ORVOC', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_MBO', 'flx_MBO', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_methanol', 'flx_methanol', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_acetone', 'flx_acetone', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(1), dt,dw2)   !! LDF TEST modificato a 1... !!
                CALL histdef(hist2_id, 'flx_acetal', 'flx_acetal', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(1), dt,dw2)   !! LDF TEST modificato a 1... !!
                CALL histdef(hist2_id, 'flx_formal', 'flx_formal', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_acetic', 'flx_acetic', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_formic', 'flx_formic', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_no_soil', 'flx_no_soil', 'ngN/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_no', 'flx_no', 'ngN/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
             ENDIF
         ELSE 
             !-
             !- This is the ALMA convention output now
             !-
             !- 
             IF ( GridType == "RegLonLat" ) THEN
#ifdef CPP_PARA
                CALL histbeg(histname2, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id,orch_domain_id)
#else
                CALL histbeg(histname2, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id)
#endif
                IF (printlev >= 2) WRITE(numout,*)  'HISTBEG2 --->',istp_old,date0,dt,dw2,hist2_id
             ELSE
#ifdef CPP_PARA
                CALL histbeg(histname2, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id,domain_id=orch_domain_id)
#else
                CALL histbeg(histname2, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id)
#endif 
             ENDIF
             !-
             CALL histvert(hist2_id, 'veget', 'Vegetation types', '1', &
                  &    nvm,   veg, vegax_id2)
             CALL histvert(hist2_id, 'solth', 'Soil levels',      'm', &
                  &    ngrnd, sol, solax_id2)
             CALL histvert(hist2_id, 'soiltyp', 'Soil types',      '1', &
                  &    nstm, soltyp, soltax_id2)
             CALL histvert(hist2_id, 'nobio', 'Other surface types',      '1', &
                  &    nnobio, nobiotyp, nobioax_id2)
             IF (  hydrol_cwrr ) THEN
                CALL histvert(hist2_id, 'solay', 'Hydrol soil levels',      'm', &
                     &    nslm, diaglev(1:nslm), solayax_id2)
             ENDIF
             !-
             !-  Vegetation
             !-
             CALL histdef(hist2_id, 'vegetfrac', 'Fraction of vegetation', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(3), dt, dw2)
             CALL histdef(hist2_id, 'maxvegetfrac', 'Maximum fraction of vegetation', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(3), dt, dw2)
             CALL histdef(hist2_id, 'nobiofrac', 'Fraction of other surface types', '1', &
                  & iim,jjm, hori_id2, nnobio, 1, nnobio, nobioax_id2, 32, avescatter2(3), dt, dw2)
             !- 
             !-  General energy balance
             !-
             CALL histdef(hist2_id, 'SWnet', 'Net shortwave radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'LWnet', 'Net longwave radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qh', 'Sensible heat flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qle', 'Latent heat flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qg', 'Ground heat flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qf', 'Energy of fusion', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'Qv', 'Energy of sublimation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'DelSurfHeat', 'Change in surface layer heat', 'J/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, sumscatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'DelColdCont', 'Change in snow surface layer cold content', 'J/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, sumscatter2(7), dt, dw2)
             !-
             !- General water balance
             !-
             CALL histdef(hist2_id, 'Snowf', 'Snowfall rate', 'kg/m^2/s',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Rainf', 'Rainfall rate', 'kg/m^2/s',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Evap', 'Total Evapotranspiration', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qs', 'Surface runoff', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qsb', 'Sub-surface runoff', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qsm', 'Snowmelt', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'DelSoilMoist', 'Change in soil moisture', 'kg/m^2',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, sumscatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'DelSurfStor', 'Change in Surface Water Storage','kg/m^2',&
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, sumscatter2(7), dt,dw2)
             CALL histdef(hist2_id, 'DelIntercept', 'Change in interception storage', 'kg/m^2',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, sumscatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'DelSWE', 'Change in interception storage Snow Water Equivalent', 'kg/m^2',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, sumscatter2(7), dt, dw2)
             !-
             !- Surface state
             !-
             CALL histdef(hist2_id, 'AvgSurfT', 'Average surface temperature', 'K', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'RadT', 'Surface radiative temperature', 'K', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'Albedo', 'Albedo', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'SurfStor', 'Surface Water Storage','kg/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'SWE', 'Snow Water Equivalent', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter2(1), dt,dw2)
             !!-
             !-  Sub-surface state
             !- 
             IF ( .NOT. hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'SoilMoist', '3D average layer soil moisture', 'kg/m^2',  &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(7), dt, dw2)
             ELSE
                CALL histdef(hist2_id, 'SoilMoist', '3D average layer soil moisture', 'kg/m^2',  &
                     & iim,jjm, hori_id2, nslm, 1, nslm, solayax_id2, 32, avescatter2(7), dt, dw2)
             ENDIF
             CALL histdef(hist2_id, 'SoilWet', 'Total soil wetness', '-',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'SoilTemp', 'Soil temperature profile', 'K', &
                  & iim,jjm, hori_id2, ngrnd, 1, ngrnd, solax_id2, 32, avescatter2(7), dt, dw2)
             !- 
             !-  Evaporation components
             !-
             CALL histdef(hist2_id, 'PotEvap', 'Potential evapotranspiration', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'ECanop', 'Interception evaporation', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, fluxop_scinsec2(5), dt, dw2)
             CALL histdef(hist2_id, 'TVeg', 'Transpiration', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, fluxop_scinsec2(5), dt, dw2)
             CALL histdef(hist2_id, 'ESoil', 'Bare soil evaporation', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(5), dt, dw2)
             CALL histdef(hist2_id, 'RootMoist','Root zone soil water', 'kg/m^2',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'SubSnow','Snow sublimation','kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'ACond', 'Aerodynamic conductance', 'm/s',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             !-
             !-
             !-  Cold Season Processes
             !-
             CALL histdef(hist2_id, 'SnowFrac', 'Snow cover fraction', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'SAlbedo', 'Snow albedo', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'SnowDepth', '3D snow depth', 'm', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             !-
             !- Hydrologic variables
             !-
             IF ( river_routing ) THEN
                !
                IF (do_floodplains) THEN
                   CALL histdef(hist2_id, 'EWater', 'Open water evaporation', 'kg/m^2/s', &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(5), dt, dw2)
                   CALL histdef(hist2_id, 'FloodFrac', 'Floodplain Fraction', '-', &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt,dw2)
                ENDIF
                !
                CALL histdef(hist2_id, 'IrrigationMap', 'Map of irrigated areas', 'm^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt, dw2)
                CALL histdef(hist2_id, 'FloodplainsMap', 'Map of flooded areas', 'm^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt,dw2)
                CALL histdef(hist2_id, 'SwampMap', 'Map of swamp areas', 'm^2', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt,dw2)
                CALL histdef(hist2_id, 'Dis', 'Simulated River Discharge', 'm^3/s', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt,dw2)
             ENDIF
             !-
             !-
             !-
             CALL histdef(hist2_id, 'humrel', 'Soil moisture stress', '-',  &
                  & iim,jjm, hori_id2, nvm,1,nvm, vegax_id2, 32, avescatter2(9), dt, dw2)
             !-
             !-  The carbon budget
             !-
             IF ( ok_co2 ) THEN
                CALL histdef(hist2_id, 'GPP', 'Net assimilation of carbon by the vegetation', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt, dw2)
             ENDIF
             IF ( ok_stomate ) THEN
                CALL histdef(hist2_id, 'NEE', 'Net Ecosystem Exchange', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt,dw2)
                CALL histdef(hist2_id, 'maint_resp', 'Maintenance respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt,dw2)
                CALL histdef(hist2_id, 'hetero_resp', 'Heterotrophic respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt,dw2)
                CALL histdef(hist2_id, 'growth_resp', 'Growth respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt, dw2)
                CALL histdef(hist2_id, 'npp', 'Net Primary Production', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt, dw2)
             ENDIF
             !
          ENDIF
          !-
          CALL histdef(hist2_id, 'LandPoints', 'Land Points', '1', &
               & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(10), dt, dw2)  
          CALL histdef(hist2_id, 'Areas', 'Mesh areas', 'm2', &
               & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt, dw2)
          CALL histdef(hist2_id, 'Contfrac', 'Continental fraction', '1', &
               & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt, dw2)
          !-
          ! Write the names of the pfts in the high frequency sechiba history files
          global_attribute="PFT_name"
          DO i=1,nvm
             WRITE(global_attribute(9:10),"(I2.2)") i
             CALL histglobal_attr(hist2_id, global_attribute, PFT_name(i))
          ENDDO
          !-
          CALL histend(hist2_id)
      ENDIF
  ENDIF

    !-
    !=====================================================================
    !- 3.2 STOMATE's history file
    !=====================================================================
    IF ( ok_stomate ) THEN
       !-
       ! STOMATE IS ACTIVATED
       !-
       !Config Key   = STOMATE_OUTPUT_FILE
       !Config Desc  = Name of file in which STOMATE's output is going to be written
       !Config If    = OK_STOMATE
       !Config Def   = stomate_history.nc
       !Config Help  = This file is going to be created by the model
       !Config         and will contain the output from the model.
       !Config         This file is a truly COADS compliant netCDF file.
       !Config         It will be generated by the hist software from
       !Config         the IOIPSL package.
       !Config Units = [FILE]
       !-
       stom_histname='stomate_history.nc'
       CALL getin_p('STOMATE_OUTPUT_FILE', stom_histname)       
       IF (printlev >= 2) WRITE(numout,*) 'STOMATE_OUTPUT_FILE', TRIM(stom_histname)
       !-
       !Config Key   = STOMATE_HIST_DT
       !Config Desc  = STOMATE history time step
       !Config If    = OK_STOMATE
       !Config Def   = 10.
       !Config Help  = Time step of the STOMATE history file
       !Config Units = [days]
       !-
       hist_days_stom = 10.
       CALL getin_p('STOMATE_HIST_DT', hist_days_stom)       

       IF ( hist_id < 0 ) THEN
          ! Deactivate all history files if sechiba_history file is deactivated
          hist_dt_stom=0
          IF (printlev >= 2) WRITE(numout,*) &
               'STOMATE history file will not be created because sechiba_history file is deactivated.'
       ELSE IF ( hist_days_stom == moins_un ) THEN
          hist_dt_stom = moins_un
          IF (printlev >= 2) WRITE(numout,*) 'output frequency for STOMATE history file (d): one month.'
       ELSE IF ( hist_days_stom == 0 ) THEN
          ! Deactivate this file
          hist_dt_stom=0
          IF (printlev >= 2) WRITE(numout,*) 'STOMATE history file will not be created'
       ELSE
          hist_dt_stom = NINT( hist_days_stom ) * one_day
          IF (printlev >= 2) WRITE(numout,*) 'output frequency for STOMATE history file (d): ', &
               hist_dt_stom/one_day
       ENDIF

       ! test consistency between STOMATE_HIST_DT and DT_STOMATE parameters
       dt_stomate_loc = one_day
       CALL getin_p('DT_STOMATE', dt_stomate_loc)
       IF ( hist_days_stom /= moins_un .AND. hist_dt_stom/=0) THEN
          IF (dt_stomate_loc > hist_dt_stom) THEN
             IF (printlev >= 2) WRITE(numout,*) "DT_STOMATE = ",dt_stomate_loc,"  , STOMATE_HIST_DT = ",hist_dt_stom
             CALL ipslerr_p (3,'ioipslctrl_history', &
                  &          'Problem with DT_STOMATE > STOMATE_HIST_DT','', &
                  &          '(must be less or equal)')
          ENDIF
       ENDIF
       !-
       !- Initialize stomate_history file
       IF ( hist_dt_stom == 0 ) THEN
          ! Case hist_dt_stom=0 : No creation of stomate_history.nc file
          ! Nothing will be done.
          hist_id_stom=-1
       ELSE
          ! Initialise stomate_history file
       IF (is_omp_root) THEN
          IF ( GridType == "RegLonLat" ) THEN
#ifdef CPP_PARA
             CALL histbeg(stom_histname, iim, lon_rect, jjm, lat_rect,  1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id_stom,domain_id=orch_domain_id)
#else
             CALL histbeg(stom_histname, iim, lon_rect, jjm, lat_rect,  1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id_stom)
#endif
          ELSE
#ifdef CPP_PARA
             CALL histbeg(stom_histname, iim, lon, jjm, lat,  1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id_stom,domain_id=orch_domain_id)
#else
             CALL histbeg(stom_histname, iim, lon, jjm, lat,  1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id_stom)
#endif
          ENDIF
          !- define PFT axis
          hist_PFTaxis = (/ ( REAL(i,r_std), i=1,nvm ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'PFT', 'Plant functional type', &
               & '1', nvm, hist_PFTaxis, hist_PFTaxis_id)
          ! deforestation
          !- define Pool_10 axis
          hist_pool_10axis = (/ ( REAL(i,r_std), i=1,10 ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'P10', 'Pool 10 years', &
               & '1', 10, hist_pool_10axis, hist_pool_10axis_id)
          
          !- define Pool_100 axis
          hist_pool_100axis = (/ ( REAL(i,r_std), i=1,100 ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'P100', 'Pool 100 years', &
               & '1', 100, hist_pool_100axis, hist_pool_100axis_id)
          
          !- define Pool_11 axis
          hist_pool_11axis = (/ ( REAL(i,r_std), i=1,11 ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'P11', 'Pool 10 years + 1', &
               & '1', 11, hist_pool_11axis, hist_pool_11axis_id)
          
          !- define Pool_101 axis
          hist_pool_101axis = (/ ( REAL(i,r_std), i=1,101 ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'P101', 'Pool 100 years + 1', &
               & '1', 101, hist_pool_101axis, hist_pool_101axis_id)
       ENDIF
       !- define STOMATE history file
       CALL ioipslctrl_histstom (hist_id_stom, nvm, iim, jjm, &
            & dt, hist_dt_stom, hori_id, hist_PFTaxis_id, &
            & hist_pool_10axis_id, hist_pool_100axis_id, &
            & hist_pool_11axis_id, hist_pool_101axis_id)
       
       !- Write the names of the pfts in the stomate history files 
       IF (is_omp_root) THEN
          global_attribute="PFT_name"
          DO i=1,nvm
             WRITE(global_attribute(9:10),"(I2.2)") i
             CALL histglobal_attr(hist_id_stom, global_attribute, PFT_name(i))
          ENDDO

       !- end definition
          CALL histend(hist_id_stom)
       ENDIF
    END IF ! IF ( hist_dt_stom == 0 )

       !-
       !-
       !-
       ! STOMATE IPCC OUTPUTS IS ACTIVATED
       !-
       !Config Key   = STOMATE_IPCC_OUTPUT_FILE
       !Config Desc  = Name of file in which STOMATE's output is going to be written
       !Config If    = OK_STOMATE
       !Config Def   = stomate_ipcc_history.nc
       !Config Help  = This file is going to be created by the model
       !Config         and will contain the output from the model.
       !Config         This file is a truly COADS compliant netCDF file.
       !Config         It will be generated by the hist software from
       !Config         the IOIPSL package.
       !Config Units = [FILE]
       !-
       stom_ipcc_histname='stomate_ipcc_history.nc'
       CALL getin_p('STOMATE_IPCC_OUTPUT_FILE', stom_ipcc_histname)       
       IF (printlev >= 2) WRITE(numout,*) 'STOMATE_IPCC_OUTPUT_FILE', TRIM(stom_ipcc_histname)
       !-
       !Config Key   = STOMATE_IPCC_HIST_DT
       !Config Desc  = STOMATE IPCC history time step
       !Config If    = OK_STOMATE
       !Config Def   = 0.
       !Config Help  = Time step of the STOMATE IPCC history file
       !Config Units = [days]
       !-
       hist_days_stom_ipcc = zero
       CALL getin_p('STOMATE_IPCC_HIST_DT', hist_days_stom_ipcc)       
       IF ( hist_days_stom_ipcc == moins_un ) THEN
          hist_dt_stom_ipcc = moins_un
          IF (printlev >= 2) WRITE(numout,*) 'output frequency for STOMATE IPCC history file (d): one month.'
       ELSE
          hist_dt_stom_ipcc = NINT( hist_days_stom_ipcc ) * one_day
          IF (printlev >= 2) WRITE(numout,*) 'output frequency for STOMATE IPCC history file (d): ', &
            hist_dt_stom_ipcc/one_day
       ENDIF
       
       IF ( hist_dt_stom_ipcc /= 0 .AND. hist_id < 0 ) THEN
          ! sechiba_history file is not created therefore STOMATE IPCC history file will be deactivated
          hist_dt_stom_ipcc=0
          hist_days_stom_ipcc=0
          IF (printlev >= 2) WRITE(numout,*) 'STOMATE IPCC history file is not created.'
       END IF

       ! test consistency between STOMATE_IPCC_HIST_DT and DT_STOMATE parameters
       dt_stomate_loc = one_day
       CALL getin_p('DT_STOMATE', dt_stomate_loc)
       IF ( hist_days_stom_ipcc > zero ) THEN
          IF (dt_stomate_loc > hist_dt_stom_ipcc) THEN
             IF (printlev >= 2) WRITE(numout,*) "DT_STOMATE = ",dt_stomate_loc,"  , STOMATE_IPCC_HIST_DT = ",hist_dt_stom_ipcc
             CALL ipslerr_p (3,'ioipslctrl_history', &
                  &          'Problem with DT_STOMATE > STOMATE_IPCC_HIST_DT','', &
                  &          '(must be less or equal)')
          ENDIF
       ENDIF

       !Config Key   = OK_HISTSYNC
       !Config Desc  = Syncronize and write IOIPSL output files at each time step
       !Config If    = 
       !Config Def   = FALSE
       !Config Help  = Setting this flag to true might affect run performance. Only use it for debug perpose.
       !Config Units = [FLAG]
       ok_histsync=.FALSE.
       CALL getin_p('OK_HISTSYNC', ok_histsync)       



       IF ( hist_dt_stom_ipcc == 0 ) THEN
          hist_id_stom_ipcc = -1
       ELSE
          !-
          !- initialize
          IF (is_omp_root) THEN
             IF ( GridType == "RegLonLat" ) THEN
#ifdef CPP_PARA
                CALL histbeg(stom_ipcc_histname, iim, lon_rect, jjm, lat_rect,  1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id, hist_id_stom_ipcc,domain_id=orch_domain_id)
#else
                CALL histbeg(stom_ipcc_histname, iim, lon_rect, jjm, lat_rect,  1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id, hist_id_stom_ipcc)
#endif
             ELSE
#ifdef CPP_PARA
                CALL histbeg(stom_ipcc_histname, iim, lon, jjm, lat,  1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id, hist_id_stom_ipcc,domain_id=orch_domain_id)
#else
                CALL histbeg(stom_ipcc_histname, iim, lon, jjm, lat,  1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id, hist_id_stom_ipcc)
#endif
             ENDIF
             !- declare this axis
             CALL histvert (hist_id_stom_IPCC, 'PFT', 'Plant functional type', &
                  & '1', nvm, hist_PFTaxis, hist_IPCC_PFTaxis_id)
             
             !- define STOMATE history file
             CALL ioipslctrl_histstomipcc (hist_id_stom_IPCC, nvm, iim, jjm, &
                  & dt, hist_dt_stom_ipcc, hori_id, hist_IPCC_PFTaxis_id)
             
             !- Write the names of the pfts in the stomate history files 
             global_attribute="PFT_name"
             DO i=1,nvm
                WRITE(global_attribute(9:10),"(I2.2)") i
                CALL histglobal_attr(hist_id_stom_IPCC, global_attribute, PFT_name(i))
             ENDDO

             !- end definition
             CALL histend(hist_id_stom_IPCC)
          ENDIF
      ENDIF
   ENDIF


  END SUBROUTINE ioipslctrl_history

!! ================================================================================================================================
!! SUBROUTINE    : ioipslctrl_histstom
!!
!>\BRIEF	 This subroutine initialize the IOIPSL stomate output file
!! 
!! DESCRIPTION   : This subroutine initialize the IOIPSL output file stomate_history.nc(default name).
!!                 This subroutine was previously named stom_define_history and where located in module intersurf.
!! RECENT CHANGE(S): None
!!
!! \n
!_ ================================================================================================================================
  SUBROUTINE ioipslctrl_histstom( &
       hist_id_stom, nvm, iim, jjm, dt, &
       hist_dt, hist_hori_id, hist_PFTaxis_id, &
       hist_pool_10axis_id, hist_pool_100axis_id, &
       hist_pool_11axis_id, hist_pool_101axis_id)
    ! deforestation axis added as arguments

    !---------------------------------------------------------------------
    !- Tell ioipsl which variables are to be written
    !- and on which grid they are defined
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !-
    !- Input
    !-
    !- File id
    INTEGER(i_std),INTENT(in) :: hist_id_stom
    !- number of PFTs
    INTEGER(i_std),INTENT(in) :: nvm
    !- Domain size
    INTEGER(i_std),INTENT(in) :: iim, jjm
    !- Time step of STOMATE (seconds)
    REAL(r_std),INTENT(in)    :: dt
    !- Time step of history file (s)
    REAL(r_std),INTENT(in)    :: hist_dt
    !- id horizontal grid
    INTEGER(i_std),INTENT(in) :: hist_hori_id
    !- id of PFT axis
    INTEGER(i_std),INTENT(in) :: hist_PFTaxis_id
    !- id of Deforestation axis
    INTEGER(i_std),INTENT(in) :: hist_pool_10axis_id,hist_pool_100axis_id
    INTEGER(i_std),INTENT(in) :: hist_pool_11axis_id,hist_pool_101axis_id
    !-
    !- 1 local
    !-
    !- maximum history level
    INTEGER(i_std), PARAMETER  :: max_hist_level = 10
    !- output level (between 0 and 10)
    !-  ( 0:nothing is written, 10:everything is written)
    INTEGER(i_std)             :: hist_level
    !- Character strings to define operations for histdef
    CHARACTER(LEN=40),DIMENSION(max_hist_level) :: ave

    !---------------------------------------------------------------------
    !=====================================================================
    !- 1 history level
    !=====================================================================
    !- 1.1 define history levelx
    !=====================================================================
    !Config Key   = STOMATE_HISTLEVEL
    !Config Desc  = STOMATE history output level (0..10)
    !Config If    = OK_STOMATE
    !Config Def   = 10
    !Config Help  = 0: nothing is written; 10: everything is written
    !Config Units = [-]
    !-
    hist_level = 10
    CALL getin_p('STOMATE_HISTLEVEL', hist_level)
    !-
    IF (printlev >= 2) WRITE(numout,*) 'STOMATE history level: ',hist_level
    IF ( (hist_level > max_hist_level).OR.(hist_level < 0) ) THEN
       STOP 'This history level is not allowed'
    ENDIF
    !=====================================================================
    !- 1.2 define operations according to output level
    !=====================================================================
    ave(1:hist_level) =  'ave(scatter(X))'
    ave(hist_level+1:max_hist_level) =  'never          '
    !=====================================================================
    !- 2 surface fields (2d)
    !- 3 PFT: 3rd dimension
    !=====================================================================


    ! structural litter above ground
    IF (is_omp_root) THEN
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTER_STR_AB       "), &
            &               TRIM("structural litter above ground                    "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! metabolic litter above ground                     
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTER_MET_AB       "), &
            &               TRIM("metabolic litter above ground                     "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! structural litter below ground               
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTER_STR_BE       "), &
            &               TRIM("structural litter below ground                    "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! metabolic litter below ground                
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTER_MET_BE       "), &
            &               TRIM("metabolic litter below ground                     "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! fraction of soil covered by dead leaves           
       CALL histdef (hist_id_stom, &
            &               TRIM("DEADLEAF_COVER      "), &
            &               TRIM("fraction of soil covered by dead leaves           "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)
       
       ! total soil and litter carbon
       CALL histdef (hist_id_stom, &
            &               TRIM("TOTAL_SOIL_CARB     "), &
            &               TRIM("total soil and litter carbon                      "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! active soil carbon in ground                 
       CALL histdef (hist_id_stom, &
            &               TRIM("CARBON_ACTIVE       "), &
            &               TRIM("active soil carbon in ground                      "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! slow soil carbon in ground                   
       CALL histdef (hist_id_stom, &
            &               TRIM("CARBON_SLOW         "), &
            &               TRIM("slow soil carbon in ground                        "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! passive soil carbon in ground                
       CALL histdef (hist_id_stom, &
            &               TRIM("CARBON_PASSIVE      "), &
            &               TRIM("passive soil carbon in ground                     "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! Long term 2 m temperature                           
       CALL histdef (hist_id_stom, &
            &               TRIM("T2M_LONGTERM        "), &
            &               TRIM("Longterm 2 m temperature                          "), &
            &               TRIM("K                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(9), dt, hist_dt)
       
       ! Monthly 2 m temperature                           
       CALL histdef (hist_id_stom, &
            &               TRIM("T2M_MONTH           "), &
            &               TRIM("Monthly 2 m temperature                           "), &
            &               TRIM("K                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(1), dt, hist_dt)
       
       ! Weekly 2 m temperature                            
       CALL histdef (hist_id_stom, &
            &               TRIM("T2M_WEEK            "), &
            &               TRIM("Weekly 2 m temperature                            "), &
            &               TRIM("K                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(1), dt, hist_dt)
       
       ! heterotr. resp. from ground                  
       CALL histdef (hist_id_stom, &
            &               TRIM("HET_RESP            "), &
            &               TRIM("heterotr. resp. from ground                       "), &
            &               TRIM("gC/m^2 tot/pft/day  "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)
       
       ! Fire fraction on ground
       CALL histdef (hist_id_stom, &
            &               TRIM("FIREFRAC            "), &
            &               TRIM("Fire fraction on ground                           "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)

       ! Fire index on ground                      
       CALL histdef (hist_id_stom, &
            &               TRIM("FIREINDEX           "), &
            &               TRIM("Fire index on ground                              "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
       
       ! Litter humidity                                   
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTERHUM           "), &
            &               TRIM("Litter humidity                                   "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)
       
       ! CO2 flux                                  
       CALL histdef (hist_id_stom, &
            &               TRIM("CO2FLUX             "), &
            &               TRIM("CO2 flux                                          "), &
            &               TRIM("gC/m^2/pft/mth      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)

       ! Output CO2 flux from fire                         
       CALL histdef (hist_id_stom, &
            &               TRIM("CO2_FIRE            "), &
            &               TRIM("Output CO2 flux from fire                         "), &
            &               TRIM("gC/day/m^2/pft      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)
       
       ! CO2 taken from atmosphere for initiate growth     
       CALL histdef (hist_id_stom, &
            &               TRIM("CO2_TAKEN           "), &
            &               TRIM("CO2 taken from atmosphere for initiate growth     "), &
            &               TRIM("gC/day/m^2/pft      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)

       IF (ok_dgvm) THEN
          ! total co2 flux (sum over 13 PFTs). when DGVM is activated, the previous
          ! SUM(CO2FLUX*veget_max) is wrong. We should look at this variable.
          CALL histdef (hist_id_stom, &
               &               TRIM("tCO2FLUX            "), &
               &               TRIM("total CO2flux of 13 PFTs (after adjustment)       "), &
               &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(1), dt, hist_dt)
          
          ! should be the same with tCO2FLUX
          CALL histdef (hist_id_stom, &
               &               TRIM("tCO2FLUX_OLD        "), &
               &               TRIM("total CO2flux of 13 PFTs(multiply by veget_max_old"), &
               &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(10), dt, hist_dt)
          
          CALL histdef (hist_id_stom, &
               &               TRIM("tGPP                 "), &
               &               TRIM("total GPP of 13 PFTs (after adjustment)           "), &
               &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(1), dt, hist_dt)
       
          CALL histdef (hist_id_stom, &
               &               TRIM("tRESP_GROWTH         "), &
               &               TRIM("total resp growth of 13 PFTs (after adjustment)   "), &
               &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(1), dt, hist_dt)
          
          CALL histdef (hist_id_stom, &
               &               TRIM("tRESP_MAINT          "), &
               &               TRIM("total resp maint  of 13 PFTs (after adjustment)   "), &
               &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(1), dt, hist_dt)
       
          CALL histdef (hist_id_stom, &
               &               TRIM("tRESP_HETERO         "), &
               &               TRIM("total resp hetero of 13 PFTs (after adjustment)   "), &
               &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(1), dt, hist_dt)
       
          CALL histdef (hist_id_stom, &
               &               TRIM("tCARBON              "), &
               &               TRIM("total carbon of 13 PFTs (after adjustment)        "), &
               &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(1), dt, hist_dt)
          
          CALL histdef (hist_id_stom, &
               &               TRIM("tBIOMASS             "), &
               &               TRIM("total biomass of 13 PFTs (after adjustment)       "), &
               &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(1), dt, hist_dt)
       
          CALL histdef (hist_id_stom, &
               &               TRIM("tLITTER              "), &
               &               TRIM("total litter of 13 PFTs (after adjustment)        "), &
               &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(1), dt, hist_dt)
       
          CALL histdef (hist_id_stom, &
               &               TRIM("tSOILC               "), &
               &               TRIM("total soil carbon of 13 PFTs (after adjustment)   "), &
               &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(1), dt, hist_dt)

          CALL histdef (hist_id_stom, &
               &               TRIM("tCO2_TAKEN           "), &
               &               TRIM("total co2_to_bm 13 PFTs (after adjustment)        "), &
               &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(10), dt, hist_dt)
          
          CALL histdef (hist_id_stom, &
               &               TRIM("tCO2_FIRE            "), &
               &               TRIM("total co2_fire 13 PFTs (after adjustment)         "), &
               &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
               &               1,1,1, -99,32, ave(10), dt, hist_dt)
       END IF
       

       CALL histdef (hist_id_stom, &
            &               TRIM("FPC_MAX             "), &
            &               TRIM("foliage projective cover                          "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("MAXFPC_LASTYEAR     "), &
            &               TRIM("foliage projective cover of last year             "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

       ! "seasonal" 2 m temperature                           
       CALL histdef (hist_id_stom, &
         &               TRIM("TSEASON             "), &
         &               TRIM("Seasonal 2 m temperature                             "), &
         &               TRIM("K                   "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(10), dt, hist_dt)

       ! how many days after onset                           
       CALL histdef (hist_id_stom, &
         &               TRIM("TMIN_SPRING_TIME    "), &
         &               TRIM("how many days after onset                            "), &
         &               TRIM("days                "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

       !                           
       CALL histdef (hist_id_stom, &
         &               TRIM("ONSET_DATE          "), &
         &               TRIM("onset date                                           "), &
         &               TRIM("day                 "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

       ! Leaf Area Index                                   
       CALL histdef (hist_id_stom, &
            &               TRIM("LAI                 "), &
            &               TRIM("Leaf Area Index                                   "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)
       
       ! Maximum vegetation fraction (LAI -> infinity)     
       CALL histdef (hist_id_stom, &
            &               TRIM("VEGET_COV_MAX       "), &
            &               TRIM("Maximum vegetation fraction (LAI -> infinity)     "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)
       
       ! Net primary productivity                          
       CALL histdef (hist_id_stom, &
            &               TRIM("NPP                 "), &
            &               TRIM("Net primary productivity                          "), &
            &               TRIM("gC/day/(m^2 tot)    "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)

       ! Gross primary productivity                        
       CALL histdef (hist_id_stom, &
            &               TRIM("GPP                 "), &
            &               TRIM("Gross primary productivity                        "), &
            &               TRIM("gC/day/m^2          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)

       ! Density of individuals                            
       CALL histdef (hist_id_stom, &
            &               TRIM("IND                 "), &
            &               TRIM("Density of individuals                            "), &
            &               TRIM("1/ m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)

       ! Adaptation to climate
       CALL histdef (hist_id_stom, &
            &               TRIM("ADAPTATION          "), &
            &               TRIM("Adaptation to climate (DGVM)                      "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)
    
       ! Probability from regenerative
       CALL histdef (hist_id_stom, &
            &               TRIM("REGENERATION        "), &
            &               TRIM("Probability from regenerative (DGVM)               "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)
       
       ! crown area of individuals (m**2)
       CALL histdef (hist_id_stom, &
            &               TRIM("CN_IND              "), &
            &               TRIM("crown area of individuals                         "), &
            &               TRIM("m^2                 "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)

       ! woodmass of individuals (gC)
       CALL histdef (hist_id_stom, &
            &               TRIM("WOODMASS_IND        "), &
            &               TRIM("Woodmass of individuals                           "), &
            &               TRIM("gC/pft              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)

       ! total living biomass
       CALL histdef (hist_id_stom, &
            &               TRIM("TOTAL_M             "), &
            &               TRIM("Total living biomass                              "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)
       
       ! Leaf mass                                         
       CALL histdef (hist_id_stom, &
            &               TRIM("LEAF_M              "), &
            &               TRIM("Leaf mass                                         "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! Sap mass above ground                             
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_M_AB            "), &
            &               TRIM("Sap mass above ground                             "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)

       ! Sap mass below ground                             
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_M_BE            "), &
            &               TRIM("Sap mass below ground                             "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! Heartwood mass above ground                       
       CALL histdef (hist_id_stom, &
            &               TRIM("HEART_M_AB          "), &
            &               TRIM("Heartwood mass above ground                       "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)

       ! Heartwood mass below ground                       
       CALL histdef (hist_id_stom, &
            &               TRIM("HEART_M_BE          "), &
            &               TRIM("Heartwood mass below ground                       "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)

       ! Root mass                                         
       CALL histdef (hist_id_stom, &
            &               TRIM("ROOT_M              "), &
            &               TRIM("Root mass                                         "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! Fruit mass                                        
       CALL histdef (hist_id_stom, &
            &               TRIM("FRUIT_M             "), &
            &               TRIM("Fruit mass                                        "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! Carbohydrate reserve mass                         
       CALL histdef (hist_id_stom, &
            &               TRIM("RESERVE_M           "), &
            &               TRIM("Carbohydrate reserve mass                         "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! total turnover rate
       CALL histdef (hist_id_stom, &
            &               TRIM("TOTAL_TURN          "), &
            &               TRIM("total turnover rate                               "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Leaf turnover                                     
       CALL histdef (hist_id_stom, &
            &               TRIM("LEAF_TURN           "), &
            &               TRIM("Leaf turnover                                     "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Sap turnover above                                
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_AB_TURN         "), &
            &               TRIM("Sap turnover above                                "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Root turnover                                     
       CALL histdef (hist_id_stom, &
            &               TRIM("ROOT_TURN           "), &
            &               TRIM("Root turnover                                     "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Fruit turnover                                    
       CALL histdef (hist_id_stom, &
            &               TRIM("FRUIT_TURN          "), &
            &               TRIM("Fruit turnover                                    "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! total conversion of biomass to litter 
       CALL histdef (hist_id_stom, &
            &               TRIM("TOTAL_BM_LITTER     "), &
            &               TRIM("total conversion of biomass to litter             "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Leaf death                                        
       CALL histdef (hist_id_stom, &
            &               TRIM("LEAF_BM_LITTER      "), &
            &               TRIM("Leaf death                                        "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)
       
       ! Sap death above ground                            
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_AB_BM_LITTER    "), &
            &               TRIM("Sap death above ground                            "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Sap death below ground                            
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_BE_BM_LITTER    "), &
            &               TRIM("Sap death below ground                            "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Heartwood death above ground                      
       CALL histdef (hist_id_stom, &
            &               TRIM("HEART_AB_BM_LITTER  "), &
            &               TRIM("Heartwood death above ground                      "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Heartwood death below ground                      
       CALL histdef (hist_id_stom, &
            &               TRIM("HEART_BE_BM_LITTER  "), &
            &               TRIM("Heartwood death below ground                      "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Root death                                        
       CALL histdef (hist_id_stom, &
            &               TRIM("ROOT_BM_LITTER      "), &
            &               TRIM("Root death                                        "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)
       
       ! Fruit death                                       
       CALL histdef (hist_id_stom, &
            &               TRIM("FRUIT_BM_LITTER     "), &
            &               TRIM("Fruit death                                       "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Carbohydrate reserve death                        
       CALL histdef (hist_id_stom, &
            &               TRIM("RESERVE_BM_LITTER   "), &
            &               TRIM("Carbohydrate reserve death                        "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Maintenance respiration                           
       CALL histdef (hist_id_stom, &
            &               TRIM("MAINT_RESP          "), &
            &               TRIM("Maintenance respiration                           "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)

       ! Growth respiration                                
       CALL histdef (hist_id_stom, &
            &               TRIM("GROWTH_RESP         "), &
            &               TRIM("Growth respiration                                "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! age                                               
       CALL histdef (hist_id_stom, &
            &               TRIM("AGE                 "), &
            &               TRIM("age                                               "), &
            &               TRIM("years               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(7), dt, hist_dt)
       
       ! height                                            
       CALL histdef (hist_id_stom, &
            &               TRIM("HEIGHT              "), &
            &               TRIM("height                                            "), &
            &               TRIM("m                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(7), dt, hist_dt)

       ! weekly moisture stress                            
       CALL histdef (hist_id_stom, &
            &               TRIM("MOISTRESS           "), &
            &               TRIM("weekly moisture stress                            "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)

       ! Maximum rate of carboxylation                     
       CALL histdef (hist_id_stom, &
            &               TRIM("VCMAX               "), &
            &               TRIM("Maximum rate of carboxylation                     "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! leaf age                                          
       CALL histdef (hist_id_stom, &
            &               TRIM("LEAF_AGE            "), &
            &               TRIM("leaf age                                          "), &
            &               TRIM("days                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)
       
       ! Fraction of trees that dies (gap)                 
       CALL histdef (hist_id_stom, &
            &               TRIM("MORTALITY           "), &
            &               TRIM("Fraction of trees that dies (gap)                 "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! Fraction of plants killed by fire                 
       CALL histdef (hist_id_stom, &
            &               TRIM("FIREDEATH           "), &
            &               TRIM("Fraction of plants killed by fire                 "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! Density of newly established saplings             
       CALL histdef (hist_id_stom, &
            &               TRIM("IND_ESTAB           "), &
            &               TRIM("Density of newly established saplings             "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! Establish tree
       CALL histdef (hist_id_stom, &
            &               TRIM("ESTABTREE           "), &
            &               TRIM("Rate of tree establishement                       "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(6), dt, hist_dt)

       ! Establish grass
       CALL histdef (hist_id_stom, &
            &               TRIM("ESTABGRASS          "), &
            &               TRIM("Rate of grass establishement                      "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(6), dt, hist_dt)

       ! Fraction of plants that dies (light competition)  
       CALL histdef (hist_id_stom, &
            &               TRIM("LIGHT_DEATH         "), &
            &               TRIM("Fraction of plants that dies (light competition)  "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! biomass allocated to leaves                       
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_LEAF       "), &
            &               TRIM("biomass allocated to leaves                       "), &
            &               TRIM("gC/m**2/pft/dt      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to sapwood above ground         
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_SAP_AB     "), &
            &               TRIM("biomass allocated to sapwood above ground         "), &
            &               TRIM("gC/m**2/pft/dt      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to sapwood below ground         
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_SAP_BE     "), &
            &               TRIM("biomass allocated to sapwood below ground         "), &
            &               TRIM("gC/m**2/pft/dt      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to roots                        
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_ROOT       "), &
            &               TRIM("biomass allocated to roots                        "), &
            &               TRIM("gC/m**2/pft/dt          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to fruits                       
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_FRUIT      "), &
            &               TRIM("biomass allocated to fruits                       "), &
            &               TRIM("gC/m**2/pft/dt          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to carbohydrate reserve         
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_RES        "), &
            &               TRIM("biomass allocated to carbohydrate reserve         "), &
            &               TRIM("gC/m**2/pft/dt          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! time constant of herbivore activity               
       CALL histdef (hist_id_stom, &
            &               TRIM("HERBIVORES          "), &
            &               TRIM("time constant of herbivore activity               "), &
            &               TRIM("days                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! turnover time for grass leaves                    
       CALL histdef (hist_id_stom, &
            &               TRIM("TURNOVER_TIME       "), &
            &               TRIM("turnover time for grass leaves                    "), &
            &               TRIM("days                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! 10 year wood product pool                         
       CALL histdef (hist_id_stom, &
            &               TRIM("PROD10              "), &
            &               TRIM("10 year wood product pool                         "), &
            &               TRIM("gC/m**2             "), iim,jjm, hist_hori_id, &
            &               11,1,11, hist_pool_11axis_id,32, ave(5), dt, hist_dt)
       
       ! annual flux for each 10 year wood product pool    
       CALL histdef (hist_id_stom, &
            &               TRIM("FLUX10              "), &
            &               TRIM("annual flux for each 10 year wood product pool    "), &
            &               TRIM("gC/m**2/yr          "), iim,jjm, hist_hori_id, &
            &               10,1,10, hist_pool_10axis_id,32, ave(5), dt, hist_dt)
       
       ! 100 year wood product pool                        
       CALL histdef (hist_id_stom, &
            &               TRIM("PROD100             "), &
            &               TRIM("100 year wood product pool                        "), &
            &               TRIM("gC/m**2             "), iim,jjm, hist_hori_id, &
            &               101,1,101, hist_pool_101axis_id,32, ave(5), dt, hist_dt)

       ! annual flux for each 100 year wood product pool   
       CALL histdef (hist_id_stom, &
            &               TRIM("FLUX100             "), &
            &               TRIM("annual flux for each 100 year wood product pool   "), &
            &               TRIM("gC/m**2/yr          "), iim,jjm, hist_hori_id, &
            &               100,1,100, hist_pool_100axis_id,32, ave(5), dt, hist_dt)

       ! annual release right after deforestation          
       CALL histdef (hist_id_stom, &
            &               TRIM("CONVFLUX            "), &
            &               TRIM("annual release right after deforestation          "), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

       ! annual release from all 10 year wood product pools 
       CALL histdef (hist_id_stom, &
            &               TRIM("CFLUX_PROD10        "), &
            &               TRIM("annual release from all 10 year wood product pools"), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

       ! annual release from all 100year wood product pools
       CALL histdef (hist_id_stom, &
            &               TRIM("CFLUX_PROD100       "), &
            &               TRIM("annual release from all 100year wood product pools"), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

       ! WOOD HARVEST
       ! 10 year wood product pool                         
       CALL histdef (hist_id_stom, &
            &               TRIM("PROD10_HARVEST      "), &
            &               TRIM("10 year wood product pool                         "), &
            &               TRIM("gC/m**2             "), iim,jjm, hist_hori_id, &
            &               11,1,11, hist_pool_11axis_id,32, ave(5), dt, hist_dt)
       
       ! annual flux for each 10 year wood product pool    
       CALL histdef (hist_id_stom, &
            &               TRIM("FLUX10_HARVEST      "), &
            &               TRIM("annual flux for each 10 year wood product pool    "), &
            &               TRIM("gC/m**2/yr          "), iim,jjm, hist_hori_id, &
            &               10,1,10, hist_pool_10axis_id,32, ave(5), dt, hist_dt)
       
       ! 100 year wood product pool                        
       CALL histdef (hist_id_stom, &
            &               TRIM("PROD100_HARVEST     "), &
            &               TRIM("100 year wood product pool                        "), &
            &               TRIM("gC/m**2             "), iim,jjm, hist_hori_id, &
            &               101,1,101, hist_pool_101axis_id,32, ave(5), dt, hist_dt)

       ! annual flux for each 100 year wood product pool   
       CALL histdef (hist_id_stom, &
            &               TRIM("FLUX100_HARVEST     "), &
            &               TRIM("annual flux for each 100 year wood product pool   "), &
            &               TRIM("gC/m**2/yr          "), iim,jjm, hist_hori_id, &
            &               100,1,100, hist_pool_100axis_id,32, ave(5), dt, hist_dt)

       ! annual release right after deforestation          
       CALL histdef (hist_id_stom, &
            &               TRIM("CONVFLUX_HARVEST      "), &
            &               TRIM("annual release right after deforestation          "), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

       ! annual release from all 10 year wood product pools 
       CALL histdef (hist_id_stom, &
            &               TRIM("CFLUX_PROD10_HARVEST   "), &
            &               TRIM("annual release from all 10 year wood product pools"), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

       ! annual release from all 100year wood product pools
       ! Note removed last letter T from HARVEST in the variable name to limit number of authorized charcters
       CALL histdef (hist_id_stom, &
            &               TRIM("CFLUX_PROD100_HARVES"), &
            &               TRIM("annual release from all 100year wood product pools"), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

       CALL histdef (hist_id_stom, &
            &               TRIM("WOOD_HARVEST  "), &
            &               TRIM("harvested wood biomass"), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

       CALL histdef (hist_id_stom, &
            &               TRIM("WOOD_HARVEST_PFT  "), &
            &               TRIM("harvested wood biomass per PFT"), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! agriculure product
       CALL histdef (hist_id_stom, &
            &               TRIM("HARVEST_ABOVE       "), &
            &               TRIM("annual release product after harvest              "), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)


       CALL histdef(hist_id_stom, 'RESOLUTION_X', 'E-W resolution', 'm', &
            & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
       CALL histdef(hist_id_stom, 'RESOLUTION_Y', 'N-S resolution', 'm', &
            & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
       CALL histdef(hist_id_stom, 'CONTFRAC', 'Continental fraction', '1', &
            & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
       CALL histdef(hist_id_stom, 'Areas', 'Mesh areas', 'm2', &
            & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
       
       !  Special outputs for phenology
       CALL histdef (hist_id_stom, &
            &               TRIM("WHEN_GROWTHINIT     "), &
            &               TRIM("Time elapsed from season beginning                "), &
            &               TRIM("d                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("PFTPRESENT          "), &
            &               TRIM("PFT exists                                        "), &
            &               TRIM("d                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("GDD_MIDWINTER       "), &
            &               TRIM("Growing degree days, since midwinter              "), &
            &               TRIM("degK                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)

       CALL histdef (hist_id_stom, &
            &               TRIM("GDD_M5_DORMANCE     "), &
            &               TRIM("Growing degree days, since dormance               "), &
            &               TRIM("degK                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("NCD_DORMANCE        "), &
            &               TRIM("Number of chilling days, since leaves were lost   "), &
            &               TRIM("d                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("ALLOW_INITPHENO     "), &
            &               TRIM("Allow to declare beginning of the growing season  "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("BEGIN_LEAVES        "), &
            &               TRIM("Signal to start putting leaves on                 "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
    ENDIF

  END SUBROUTINE ioipslctrl_histstom

!! ================================================================================================================================
!! SUBROUTINE    : ioipslctrl_histstomipcc
!!
!>\BRIEF	 This subroutine initialize the IOIPSL stomate second output file (ipcc file)
!! 
!! DESCRIPTION   : This subroutine initialize the IOIPSL stomate second output file named stomate_ipcc_history.nc(default name).
!!                 This subroutine was previously called stom_IPCC_define_history and located in module intersurf.
!!
!! RECENT CHANGE(S): None
!!
!! \n
!_ ================================================================================================================================
  SUBROUTINE ioipslctrl_histstomipcc( &
       hist_id_stom_IPCC, nvm, iim, jjm, dt, &
       hist_dt, hist_hori_id, hist_PFTaxis_id)
    ! deforestation axis added as arguments

    !---------------------------------------------------------------------
    !- Tell ioipsl which variables are to be written
    !- and on which grid they are defined
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !-
    !- Input
    !-
    !- File id
    INTEGER(i_std),INTENT(in) :: hist_id_stom_IPCC
    !- number of PFTs
    INTEGER(i_std),INTENT(in) :: nvm
    !- Domain size
    INTEGER(i_std),INTENT(in) :: iim, jjm
    !- Time step of STOMATE (seconds)
    REAL(r_std),INTENT(in)    :: dt
    !- Time step of history file (s)
    REAL(r_std),INTENT(in)    :: hist_dt
    !- id horizontal grid
    INTEGER(i_std),INTENT(in) :: hist_hori_id
    !- id of PFT axis
    INTEGER(i_std),INTENT(in) :: hist_PFTaxis_id
    !-
    !- 1 local
    !-
    !- Character strings to define operations for histdef
    CHARACTER(LEN=40),DIMENSION(max_hist_level) :: ave

    !=====================================================================
    !- 1 define operations
    !=====================================================================
    ave(1) =  'ave(scatter(X))'
    !=====================================================================
    !- 2 surface fields (2d)
    !=====================================================================
    ! Carbon in Vegetation
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cVeg"), &
         &               TRIM("Carbon in Vegetation"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Litter Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cLitter"), &
         &               TRIM("Carbon in Litter Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Soil Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cSoil"), &
         &               TRIM("Carbon in Soil Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Products of Land Use Change
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cProduct"), &
         &               TRIM("Carbon in Products of Land Use Change"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon Mass Variation
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cMassVariation"), &
         &               TRIM("Terrestrial Carbon Mass Variation"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Leaf Area Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("lai"), &
         &               TRIM("Leaf Area Fraction"), &
         &               TRIM("1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Gross Primary Production
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("gpp"), &
         &               TRIM("Gross Primary Production"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Autotrophic Respiration
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("ra"), &
         &               TRIM("Autotrophic Respiration"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Net Primary Production
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("npp"), &
         &               TRIM("Net Primary Production"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Heterotrophic Respiration
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("rh"), &
         &               TRIM("Heterotrophic Respiration"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Emission from Fire
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fFire"), &
         &               TRIM("CO2 Emission from Fire"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! CO2 Flux to Atmosphere from Crop Harvesting
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fHarvest"), &
         &               TRIM("CO2 Flux to Atmosphere from Crop Harvesting"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Flux to Atmosphere from Land Use Change
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fLuc"), &
         &               TRIM("CO2 Flux to Atmosphere from Land Use Change"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Flux to Atmosphere from Wood Harvest                                                                                
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fWoodharvest"), &
         &               TRIM("CO2 Flux to Atmosphere from Wood Harvest"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! Net Biospheric Production
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nbp"), &
         &               TRIM("Net Biospheric Production"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Total Carbon Flux from Vegetation to Litter
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fVegLitter"), &
         &               TRIM("Total Carbon Flux from Vegetation to Litter"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Total Carbon Flux from Litter to Soil
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fLitterSoil"), &
         &               TRIM("Total Carbon Flux from Litter to Soil"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! Carbon in Leaves
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cLeaf"), &
         &               TRIM("Carbon in Leaves"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Stem
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cStem"), &
         &               TRIM("Carbon in Stem"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Roots
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cRoot"), &
         &               TRIM("Carbon in Roots"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Other Living Compartments
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cMisc"), &
         &               TRIM("Carbon in Other Living Compartments"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! Carbon in Above-Ground Litter
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cLitterAbove"), &
         &               TRIM("Carbon in Above-Ground Litter"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Below-Ground Litter
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cLitterBelow"), &
         &               TRIM("Carbon in Below-Ground Litter"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Fast Soil Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cSoilFast"), &
         &               TRIM("Carbon in Fast Soil Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Medium Soil Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cSoilMedium"), &
         &               TRIM("Carbon in Medium Soil Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Slow Soil Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cSoilSlow"), &
         &               TRIM("Carbon in Slow Soil Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    !- 3 PFT: 3rd dimension
    ! Fractional Land Cover of PFT
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("landCoverFrac"), &
         &               TRIM("Fractional Land Cover of PFT"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)


    ! Total Primary Deciduous Tree Cover Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("treeFracPrimDec"), &
         &               TRIM("Total Primary Deciduous Tree Cover Fraction"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! Total Primary Evergreen Tree Cover Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("treeFracPrimEver"), &
         &               TRIM("Total Primary Evergreen Tree Cover Fraction"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! Total C3 PFT Cover Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("c3PftFrac"), &
         &               TRIM("Total C3 PFT Cover Fraction"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Total C4 PFT Cover Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("c4PftFrac"), &
         &               TRIM("Total C4 PFT Cover Fraction"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Growth Autotrophic Respiration
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("rGrowth"), &
         &               TRIM("Growth Autotrophic Respiration"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Maintenance Autotrophic Respiration
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("rMaint"), &
         &               TRIM("Maintenance Autotrophic Respiration"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Flux from Atmosphere due to NPP Allocation to Leaf
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nppLeaf"), &
         &               TRIM("CO2 Flux from Atmosphere due to NPP Allocation to Leaf"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Flux from Atmosphere due to NPP Allocation to Wood
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nppStem"), &
         &               TRIM("CO2 Flux from Atmosphere due to NPP Allocation to Stem"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Flux from Atmosphere due to NPP Allocation to Root
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nppRoot"), &
         &               TRIM("CO2 Flux from Atmosphere due to NPP Allocation to Root"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Net Carbon Mass Flux out of Atmophere due to Net Ecosystem Productivity on Land.
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nep"), &
         &               TRIM("Net Carbon Mass Flux out of Atmophere due to Net Ecosystem Productivity."), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    CALL histdef(hist_id_stom_IPCC, 'RESOLUTION_X', 'E-W resolution', 'm', &
         & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
    CALL histdef(hist_id_stom_IPCC, 'RESOLUTION_Y', 'N-S resolution', 'm', &
         & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
    CALL histdef(hist_id_stom_IPCC, 'CONTFRAC', 'Continental fraction', '1', &
         & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
    CALL histdef(hist_id_stom_IPCC, 'Areas', 'Mesh areas', 'm2', &
         & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)

  END SUBROUTINE ioipslctrl_histstomipcc

!! ================================================================================================================================
!! SUBROUTINE    : ioipslctrl_restini
!!
!>\BRIEF	 This subroutine initialize the restart files in ORCHDIEE. 
!! 
!! DESCRIPTION   : This subroutine initialize restart files in ORCHIDEE. IOIPSL is used for manipulating the restart files. 
!!                 This subroutine was previously called intsurf_restart and located in module intersurf.
!!
!! RECENT CHANGE(S): None
!!
!! \n
!_ ================================================================================================================================
  SUBROUTINE ioipslctrl_restini(istp, date0, dt, rest_id, rest_id_stom, itau_offset, date0_shifted)

    USE mod_orchidee_para
    !
    !  This subroutine initialized the restart file for the land-surface scheme
    !
    IMPLICIT NONE
    !
    INTEGER(i_std), INTENT(in)                  :: istp      !! Time step of the restart file
    REAL(r_std)                                 :: date0     !! The date at which itau = 0
    REAL(r_std)                                 :: dt        !! Time step
    INTEGER(i_std), INTENT(out)                 :: rest_id, rest_id_stom   !! ID of the restart file
    INTEGER(i_std), INTENT(out)                 :: itau_offset    !! Note the result is always itau_offset=0 as overwrite_time=TRUE
    REAL(r_std), INTENT(out)                    :: date0_shifted  !! Note the result is always date0_shifted=date0 as overwrite_time=TRUE


    !  LOCAL
    !
    REAL(r_std)                 :: dt_rest, date0_rest
    INTEGER(i_std)              :: itau_dep
    INTEGER(i_std),PARAMETER    :: llm=1
    REAL(r_std), DIMENSION(llm) :: lev
    LOGICAL, PARAMETER          :: overwrite_time=.TRUE. !! Always override the date from the restart files for SECHIBA and STOMATE. 
                                                         !! The date is taken from the gcm or from the driver restart file. 
    REAL(r_std)                 :: in_julian, rest_julian
    INTEGER(i_std)              :: yy, mm, dd
    REAL(r_std)                 :: ss
    !
    !Config Key   = SECHIBA_restart_in
    !Config Desc  = Name of restart to READ for initial conditions
    !Config If    = OK_SECHIBA 
    !Config Def   = NONE
    !Config Help  = This is the name of the file which will be opened
    !Config         to extract the initial values of all prognostic
    !Config         values of the model. This has to be a netCDF file.
    !Config         Not truly COADS compliant. NONE will mean that
    !Config         no restart file is to be expected.
    !Config Units = [FILE]
!-
    CALL getin_p('SECHIBA_restart_in',restname_in)
    IF (printlev >= 2) WRITE(numout,*) 'Restart file for sechiba: ', restname_in
    !-
    !Config Key   = SECHIBA_rest_out
    !Config Desc  = Name of restart files to be created by SECHIBA
    !Config If    = OK_SECHIBA
    !Config Def   = sechiba_rest_out.nc
    !Config Help  = This variable give the name for
    !Config         the restart files. The restart software within
    !Config         IOIPSL will add .nc if needed.
    !Config Units = [FILE]
    !
    CALL getin_p('SECHIBA_rest_out', restname_out)
 
    lev(:) = zero
    itau_dep = istp
    in_julian = itau2date(istp, date0, dt)
    date0_rest = date0
    dt_rest = dt
    !
    IF (is_root_prc) THEN
       CALL restini( restname_in, iim_g, jjm_g, lon_g, lat_g, llm, lev, &
            &  restname_out, itau_dep, date0_rest, dt_rest, rest_id, overwrite_time)
    ELSE
       rest_id=0
    ENDIF
    CALL bcast (itau_dep)
    CALL bcast (date0_rest)
    CALL bcast (dt_rest)
    !
    !  itau_dep of SECHIBA is phased with the GCM if needed
    !
    rest_julian = itau2date(itau_dep, date0_rest, dt_rest)

    ! Note by JG
    ! restini never modifies itau_dep and date0_rest when overwrite_time=TRUE. 
    ! This means that itau_dep=istp and date0_rest=date0 => rest_julian=in_julian. 
    ! The result of below IF will therfor always be itau_offset=0 and date0_shifted=date0
    IF ( ABS( in_julian - rest_julian) .GT. dt/one_day .AND. .NOT. OFF_LINE_MODE ) THEN
       WRITE(numout,*) 'The SECHIBA restart is not for the same timestep as the GCM,'
       WRITE(numout,*) 'the two are synchronized. The land-surface conditions can not impose'
       WRITE(numout,*) 'the chronology of the simulation.'
       WRITE(numout,*) 'Time step of the GCM :', istp, 'Julian day : ', in_julian
       CALL ju2ymds(in_julian, yy, mm, dd, ss)
       WRITE(numout,*) 'In other word (yy,mm,dd,ss) : ', yy, mm, dd, ss
       WRITE(numout,*) 'Time step of SECHIBA :', itau_dep, 'Julian day : ', rest_julian
       CALL ju2ymds(rest_julian, yy, mm, dd, ss)
       WRITE(numout,*) 'In other word (yy,mm,dd,ss) : ', yy, mm, dd, ss
       
       itau_offset = itau_dep - istp
       date0_shifted = date0 - itau_offset*dt/one_day
       
       WRITE(numout,*) 'The new starting date is :', date0_shifted
       CALL ju2ymds(date0_shifted, yy, mm, dd, ss)
       WRITE(numout,*) 'In other word (yy,mm,dd,ss) : ', yy, mm, dd, ss
    ELSE
       itau_offset = 0
       date0_shifted = date0
    ENDIF

    !=====================================================================
    !- 1.5 Restart file for STOMATE
    !=====================================================================
    IF ( ok_stomate ) THEN 
       !-
       ! STOMATE IS ACTIVATED
       !-
       !Config Key   = STOMATE_RESTART_FILEIN
       !Config Desc  = Name of restart to READ for initial conditions of STOMATE
       !Config If    = STOMATE_OK_STOMATE
       !Config Def   = NONE
       !Config Help  = This is the name of the file which will be opened
       !Config         to extract the initial values of all prognostic
       !Config         values of STOMATE.
       !Config Units = [FILE]
       !-
       CALL getin_p('STOMATE_RESTART_FILEIN',stom_restname_in)
       IF (printlev >= 2) WRITE(numout,*) 'STOMATE INPUT RESTART_FILE', stom_restname_in
       !-
       !Config Key   = STOMATE_RESTART_FILEOUT
       !Config Desc  = Name of restart files to be created by STOMATE
       !Config If    = STOMATE_OK_STOMATE
       !Config Def   = stomate_rest_out.nc
       !Config Help  = This is the name of the file which will be opened
       !Config         to write the final values of all prognostic values
       !Config         of STOMATE.
       !Config Units = [FILE]
       !-
       CALL getin_p('STOMATE_RESTART_FILEOUT', stom_restname_out)
       IF (printlev >= 2) WRITE(numout,*) 'STOMATE OUTPUT RESTART_FILE', stom_restname_out
       !-
       IF (is_root_prc) THEN
         CALL restini (stom_restname_in, iim_g, jjm_g, lon_g, lat_g, llm, lev, &
            &  stom_restname_out, itau_dep, date0_rest, dt_rest, rest_id_stom, overwrite_time)
       ELSE
         rest_id_stom=0
       ENDIF
       CALL bcast (itau_dep)
       CALL bcast (date0_rest)
       CALL bcast (dt_rest)
       !-
    ENDIF
  END SUBROUTINE ioipslctrl_restini

END MODULE ioipslctrl
