!
! $Id: surface_data.F90 3324 2018-05-15 15:56:55Z musat $
!
MODULE surface_data

  REAL, PARAMETER        :: calice=1.0/(5.1444e+06*0.15)
  REAL, PARAMETER        :: calsno=1./(2.3867e+06*.15)
  
  LOGICAL, SAVE          :: ok_veget      ! true for use of vegetation model ORCHIDEE
  !$OMP THREADPRIVATE(ok_veget)

  CHARACTER(len=10), SAVE :: type_veget   ! orchidee/y/bucket/n/betaclim
  !$OMP THREADPRIVATE(type_veget)

  LOGICAL, SAVE          :: ok_snow       ! true for coupling to snow model SISVAT
  !$OMP THREADPRIVATE(ok_snow)

  CHARACTER(len=6), SAVE :: type_ocean    ! force/slab/couple
  !$OMP THREADPRIVATE(type_ocean)

  ! if type_ocean=couple : version_ocean=opa8 ou nemo
  ! if type_ocean=slab   : version_ocean=sicOBS or sicINT or sicNO
  CHARACTER(len=6), SAVE :: version_ocean 
  !$OMP THREADPRIVATE(version_ocean)

  ! Pas de temps couplage atm/oce (en secondes)
  REAL, SAVE             :: t_coupl
  !$OMP THREADPRIVATE(t_coupl)

END MODULE surface_data
