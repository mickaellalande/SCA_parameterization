!
! $Id: surf_seaice_mod.F90 3102 2017-12-03 20:27:42Z oboucher $
!
MODULE surf_seaice_mod

  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE surf_seaice( & 
       rlon, rlat, swnet, lwnet, alb1, fder, &
       itime, dtime, jour, knon, knindex, &
       lafin, &
       tsurf, p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, gustiness, pctsrf, &
       snow, qsurf, qsol, agesno, tsoil, &
       z0m, z0h, SFRWL, alb_dir_new, alb_dif_new, evap, fluxsens, fluxlat, &  
       tsurf_new, dflux_s, dflux_l, &
       flux_u1, flux_v1)

  USE dimphy
  USE surface_data
  USE ocean_forced_mod, ONLY : ocean_forced_ice
  USE ocean_cpl_mod, ONLY    : ocean_cpl_ice
  USE ocean_slab_mod, ONLY   : ocean_slab_ice
  USE indice_sol_mod

!
! This subroutine will make a call to ocean_XXX_ice according to the ocean mode (force, 
! slab or couple). The calculation of rugosity for the sea-ice surface is also done
! in here because it is the same calculation for the different modes of ocean.
!
    INCLUDE "dimsoil.h"
    INCLUDE "clesphys.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN)     :: knindex
    LOGICAL, INTENT(IN)                      :: lafin
    REAL, INTENT(IN)                         :: dtime
    REAL, DIMENSION(klon), INTENT(IN)        :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)        :: swnet  ! net shortwave radiation at surface  
    REAL, DIMENSION(klon), INTENT(IN)        :: lwnet  ! net longwave radiation at surface  
    REAL, DIMENSION(klon), INTENT(IN)        :: alb1   ! albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(IN)        :: fder
    REAL, DIMENSION(klon), INTENT(IN)        :: tsurf
    REAL, DIMENSION(klon), INTENT(IN)        :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)        :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)        :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)        :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)        :: ps
    REAL, DIMENSION(klon), INTENT(IN)        :: u1, v1, gustiness
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)  :: pctsrf

! In/Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)          :: snow, qsurf, qsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: agesno
    REAL, DIMENSION(klon, nsoilmx), INTENT(INOUT) :: tsoil

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: z0m, z0h
!albedo SB >>>
!    REAL, DIMENSION(klon), INTENT(OUT)       :: alb1_new  ! new albedo in visible SW interval
!    REAL, DIMENSION(klon), INTENT(OUT)       :: alb2_new  ! new albedo in near IR interval
    REAL, DIMENSION(6), INTENT(IN)    :: SFRWL
    REAL, DIMENSION(klon,nsw), INTENT(OUT)   :: alb_dir_new,alb_dif_new
!albedo SB <<<
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)       :: flux_u1, flux_v1

! Local arguments
!****************************************************************************************
    REAL, DIMENSION(klon)  :: radsol

!albedo SB >>>
    REAL, DIMENSION(klon) :: alb1_new,alb2_new
!albedo SB <<<
!
! End definitions
!****************************************************************************************


!****************************************************************************************
! Calculate total net radiance at surface
!
!****************************************************************************************
    radsol(:) = 0.0
    radsol(1:knon) = swnet(1:knon) + lwnet(1:knon)

!****************************************************************************************
! Switch according to type of ocean (couple, slab or forced)
!
!****************************************************************************************
    IF (type_ocean == 'couple') THEN
       
       CALL ocean_cpl_ice( &
            rlon, rlat, swnet, lwnet, alb1, & 
            fder, & 
            itime, dtime, knon, knindex, &
            lafin,&
            p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum,&
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, gustiness, pctsrf, &
            radsol, snow, qsurf, &
            alb1_new, alb2_new, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l)
       
    ELSE IF (type_ocean == 'slab'.AND.version_ocean=='sicINT') THEN
       CALL ocean_slab_ice( & 
          itime, dtime, jour, knon, knindex, &
          tsurf, p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum,&
          AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
          ps, u1, v1, gustiness, &
          radsol, snow, qsurf, qsol, agesno, &
          alb1_new, alb2_new, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
          tsurf_new, dflux_s, dflux_l, swnet)

      ELSE ! type_ocean=force or slab +sicOBS or sicNO
       CALL ocean_forced_ice( &
            itime, dtime, jour, knon, knindex, &
            tsurf, p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum,&
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, gustiness, &
            radsol, snow, qsol, agesno, tsoil, &
            qsurf, alb1_new, alb2_new, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l)

    END IF

!****************************************************************************************
! Calculate rugosity
!
!****************************************************************************************

    z0m=z0m_seaice
    z0h = z0h_seaice

!albedo SB >>>
     select case(NSW)
     case(2)
       alb_dir_new(1:knon,1)=alb1_new(1:knon)
       alb_dir_new(1:knon,2)=alb2_new(1:knon)
     case(4)
       alb_dir_new(1:knon,1)=alb1_new(1:knon)
       alb_dir_new(1:knon,2)=alb2_new(1:knon)
       alb_dir_new(1:knon,3)=alb2_new(1:knon)
       alb_dir_new(1:knon,4)=alb2_new(1:knon)
     case(6)
       alb_dir_new(1:knon,1)=alb1_new(1:knon)
       alb_dir_new(1:knon,2)=alb1_new(1:knon)
       alb_dir_new(1:knon,3)=alb1_new(1:knon)
       alb_dir_new(1:knon,4)=alb2_new(1:knon)
       alb_dir_new(1:knon,5)=alb2_new(1:knon)
       alb_dir_new(1:knon,6)=alb2_new(1:knon)
     end select
alb_dif_new=alb_dir_new
!albedo SB <<<




  END SUBROUTINE surf_seaice
!
!****************************************************************************************
!
END MODULE surf_seaice_mod

