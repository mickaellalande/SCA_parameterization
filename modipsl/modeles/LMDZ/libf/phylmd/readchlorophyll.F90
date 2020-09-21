!
! $Id$
!

subroutine readchlorophyll(debut)

    use netcdf95, only: nf95_close, nf95_gw_var, nf95_inq_dimid, & 
                        nf95_inq_varid, nf95_open
    use netcdf, only: nf90_get_var, nf90_noerr, nf90_nowrite

    USE phys_cal_mod, ONLY : mth_cur
    USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, klon_glo, &
                                 grid2dto1d_glo
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    USE mod_phys_lmdz_para, ONLY: scatter 
    USE phys_state_var_mod, ONLY: chl_con

    implicit none

    include "YOMCST.h"

! Variable input
    logical debut

! Variables locales
    integer n_lat   ! number of latitudes in the input data
    integer n_lon   ! number of longitudes in the input data
    integer n_lev   ! number of levels in the input data
    integer n_month ! number of months in the input data
    real, pointer:: latitude(:)
    real, pointer:: longitude(:)
    real, pointer:: time(:)
    integer i, k
    integer, save :: mth_pre
!$OMP THREADPRIVATE(mth_pre)

! Champs reconstitues
    real, allocatable:: chlorocon(:, :, :)
    real, allocatable:: chlorocon_mois(:, :)
    real, allocatable:: chlorocon_mois_glo(:)

! For NetCDF:
    integer ncid_in  ! IDs for input files
    integer varid, ncerr


!--------------------------------------------------------


!--only read file if beginning of run or start of new month
    IF (debut.OR.mth_cur.NE.mth_pre) THEN

    IF (is_mpi_root) THEN


    CALL nf95_open("chlorophyll.nc", nf90_nowrite, ncid_in)

    CALL nf95_inq_varid(ncid_in, "lon", varid)
    CALL nf95_gw_var(ncid_in, varid, longitude)
    n_lon = size(longitude)
!    print *, 'LON chlorophyll=', n_lon, longitude
    IF (n_lon.NE.nbp_lon) THEN
       print *,'Le nombre de lon n est pas egal a nbp_lon'
       STOP
    ENDIF


    CALL nf95_inq_varid(ncid_in, "lat", varid)
    CALL nf95_gw_var(ncid_in, varid, latitude)
    n_lat = size(latitude)
!    print *, 'LAT chlorophyll=', n_lat, latitude
    IF (n_lat.NE.nbp_lat) THEN 
       print *,'Le nombre de lat n est pas egal a jnbp_lat'
       STOP
    ENDIF

    CALL nf95_inq_varid(ncid_in, "time", varid)
    CALL nf95_gw_var(ncid_in, varid, time)
    n_month = size(time)
!    print *, 'TIME aerosol strato=', n_month, time
    IF (n_month.NE.12) THEN 
       print *,'Le nombre de month n est pas egal a 12'
       STOP
    ENDIF

    IF (.not.ALLOCATED(chlorocon))          ALLOCATE(chlorocon(n_lon, n_lat, n_month))
    IF (.not.ALLOCATED(chlorocon_mois))     ALLOCATE(chlorocon_mois(n_lon, n_lat))
    IF (.not.ALLOCATED(chlorocon_mois_glo)) ALLOCATE(chlorocon_mois_glo(klon_glo))

!--reading stratospheric AOD at 550 nm
    CALL nf95_inq_varid(ncid_in, "CHL", varid)
    ncerr = nf90_get_var(ncid_in, varid, chlorocon)
    print *,'code erreur readchlorophyll=', ncerr, varid

    CALL nf95_close(ncid_in)

!---select the correct month
    IF (mth_cur.LT.1.OR.mth_cur.GT.12) THEN
      print *,'probleme avec le mois dans readchlorophyll =', mth_cur
    ENDIF
    chlorocon_mois(:,:) = chlorocon(:,:,mth_cur)

!---reduce to a klon_glo grid 
    CALL grid2dTo1d_glo(chlorocon_mois,chlorocon_mois_glo)


    print*,"chrolophyll current month",mth_cur
    do i=1,klon_glo
!      if(isnan(chlorocon_mois_glo(i)))then ! isnan() is not in the Fortran standard...
!      Another way to check for NaN:
       if(chlorocon_mois_glo(i).ne.chlorocon_mois_glo(i)) then
         chlorocon_mois_glo(i)=0.
      endif
      !print*,"high chl con",i,chlorocon_mois_glo(i)
    enddo

!    DEALLOCATE(chlorocon)
!    DEALLOCATE(chlorocon_mois)
!    DEALLOCATE(chlorocon_mois_glo)
  
    ENDIF !--is_mpi_root

!--scatter on all proc
    CALL scatter(chlorocon_mois_glo,chl_con)

!--keep memory of previous month
    mth_pre=mth_cur

    ENDIF !--debut ou nouveau mois

end subroutine readchlorophyll
