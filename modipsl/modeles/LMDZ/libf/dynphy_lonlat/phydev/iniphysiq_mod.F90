!
! $Id: iniphysiq.F 1403 2010-07-01 09:02:53Z fairhead $
!
MODULE iniphysiq_mod

CONTAINS

SUBROUTINE iniphysiq(iim,jjm,nlayer, &
                     nbp, communicator, &
                     punjours, pdayref,ptimestep, &
                     rlatu,rlatv,rlonu,rlonv,aire,cu,cv, &
                     prad,pg,pr,pcpp,iflag_phys)
  USE dimphy, ONLY: init_dimphy
  USE inigeomphy_mod, ONLY: inigeomphy
  USE mod_phys_lmdz_para, ONLY: klon_omp ! number of columns (on local omp grid)
  USE infotrac, ONLY: nqtot, type_trac
  USE infotrac_phy, ONLY: init_infotrac_phy
  USE inifis_mod, ONLY: inifis
  USE phyaqua_mod, ONLY: iniaqua
  USE nrtype, ONLY: pi
  IMPLICIT NONE
  !
  !=======================================================================
  !   Initialisation of the physical constants and some positional and 
  !   geometrical arrays for the physics
  !=======================================================================
 
 
  include "iniprint.h"

  REAL,INTENT(IN) :: prad ! radius of the planet (m)
  REAL,INTENT(IN) :: pg ! gravitational acceleration (m/s2)
  REAL,INTENT(IN) :: pr ! ! reduced gas constant R/mu
  REAL,INTENT(IN) :: pcpp ! specific heat Cp
  REAL,INTENT(IN) :: punjours ! length (in s) of a standard day
  INTEGER, INTENT (IN) :: nlayer ! number of atmospheric layers
  INTEGER, INTENT (IN) :: iim ! number of atmospheric coulumns along longitudes
  INTEGER, INTENT (IN) :: jjm  ! number of atompsheric columns along latitudes
  INTEGER, INTENT(IN) :: nbp ! number of physics columns for this MPI process
  INTEGER, INTENT(IN) :: communicator ! MPI communicator
  REAL, INTENT (IN) :: rlatu(jjm+1) ! latitudes of the physics grid
  REAL, INTENT (IN) :: rlatv(jjm) ! latitude boundaries of the physics grid
  REAL, INTENT (IN) :: rlonv(iim+1) ! longitudes of the physics grid
  REAL, INTENT (IN) :: rlonu(iim+1) ! longitude boundaries of the physics grid
  REAL, INTENT (IN) :: aire(iim+1,jjm+1) ! area of the dynamics grid (m2)
  REAL, INTENT (IN) :: cu((iim+1)*(jjm+1)) ! cu coeff. (u_covariant = cu * u)
  REAL, INTENT (IN) :: cv((iim+1)*jjm) ! cv coeff. (v_covariant = cv * v)
  INTEGER, INTENT (IN) :: pdayref ! reference day of for the simulation
  REAL,INTENT(IN) :: ptimestep !physics time step (s)
  INTEGER,INTENT(IN) :: iflag_phys ! type of physics to be called

  INTEGER :: ibegin,iend,offset
  INTEGER :: i,j,k
  CHARACTER (LEN=20) :: modname='iniphysiq'
  CHARACTER (LEN=80) :: abort_message


  ! --> initialize physics distribution, global fields and geometry
  ! (i.e. things in phy_common or dynphy_lonlat)
  CALL inigeomphy(iim,jjm,nlayer, &
               nbp, communicator, &
               rlatu,rlatv, &
               rlonu,rlonv, &
               aire,cu,cv)

  ! --> now initialize things specific to the phydev physics package

!$OMP PARALLEL 

  ! Initialize physical constants in physics:
  CALL inifis(prad,pg,pr,pcpp)
  
  ! Initialize dimphy module
  CALL Init_dimphy(klon_omp,nlayer)

  ! Initialize tracer names, numbers, etc. for physics
  CALL init_infotrac_phy(nqtot,type_trac)

  ! Additional initializations for aquaplanets 
  IF (iflag_phys>=100) THEN
    CALL iniaqua(klon_omp,iflag_phys)
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE iniphysiq

END MODULE iniphysiq_mod
