!
! $Id: time_phylmdz_mod.F90 2805 2017-03-01 16:50:11Z fairhead $
!
MODULE time_phylmdz_mod

    IMPLICIT NONE
    REAL,SAVE    :: pdtphys     ! physics time step (s)
!$OMP THREADPRIVATE(pdtphys)
    INTEGER,SAVE :: day_step_phy    ! number of physical steps per day
!$OMP THREADPRIVATE(day_step_phy)
    INTEGER,SAVE :: ndays       ! number of days to run
!$OMP THREADPRIVATE(ndays)
    INTEGER,SAVE :: annee_ref   ! reference year from the origin
!$OMP THREADPRIVATE(annee_ref)
    INTEGER,SAVE :: day_ref     ! reference year of the origin
!$OMP THREADPRIVATE(day_ref)
    INTEGER,SAVE :: day_ini     ! initial day of the run starting from 1st january of annee_ref
!$OMP THREADPRIVATE(day_ini)
    INTEGER,SAVE :: day_end     ! final day of the run starting from 1st january of annee_ref
!$OMP THREADPRIVATE(day_end)
    REAL,SAVE    :: start_time  ! starting time from the begining of the initial day
!$OMP THREADPRIVATE(start_time)
    INTEGER,SAVE :: raz_date
!$OMP THREADPRIVATE(raz_date)

    INTEGER,SAVE :: itau_phy     ! number of physiq iteration from origin
!$OMP THREADPRIVATE(itau_phy)
    INTEGER,SAVE :: itaufin_phy      ! final iteration (in itau_phy steps)
!$OMP THREADPRIVATE(itaufin_phy)
    REAL,SAVE    :: current_time ! current elapsed time in seconds from the begining of the run
!$OMP THREADPRIVATE(current_time)
    

CONTAINS

  SUBROUTINE init_time(annee_ref_, day_ref_, day_ini_, start_time_, &
                       ndays_, pdtphys_)
  USE ioipsl_getin_p_mod, ONLY : getin_p
  USE phys_cal_mod, ONLY: phys_cal_init
  IMPLICIT NONE
  INCLUDE 'YOMCST.h'
    INTEGER, INTENT(IN) :: annee_ref_  
    INTEGER, INTENT(IN) :: day_ref_    
    INTEGER, INTENT(IN) :: day_ini_    
    REAL,    INTENT(IN) :: start_time_ 
    INTEGER, INTENT(IN) :: ndays_      
    REAL,    INTENT(IN) :: pdtphys_    
    
    annee_ref    = annee_ref_
    day_ref      = day_ref_
    day_ini      = day_ini_
    start_time   = start_time_
    ndays        = ndays_
    pdtphys      = pdtphys_
    
    ! Initialize module variable not inherited from dynamics
    day_step_phy = NINT(rday/pdtphys)
    day_end  = day_ini + ndays
  
    raz_date = 0
    CALL getin_p('raz_date', raz_date)

    current_time=0.
    
    CALL phys_cal_init(annee_ref,day_ref)
    
  END SUBROUTINE init_time

  SUBROUTINE init_iteration(itau_phy_)
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: itau_phy_
    itau_phy=itau_phy_
    IF (raz_date==1) itau_phy=0
    
    itaufin_phy=itau_phy+NINT(ndays/pdtphys)
    
  END SUBROUTINE init_iteration

  SUBROUTINE update_time(pdtphys_)
  ! This subroutine updates the module saved variables.
  USE ioipsl, ONLY : ymds2ju
  USE phys_cal_mod, ONLY: phys_cal_update
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE
  INCLUDE 'YOMCST.h'
  REAL,INTENT(IN) :: pdtphys_
  REAL            :: julian_date
  INTEGER         :: cur_day
  REAL            :: cur_sec

    ! Check if the physics timestep has changed
    IF ( ABS( (pdtphys-pdtphys_) / ((pdtphys+pdtphys_)/2))> 10.*EPSILON(pdtphys_)) THEN
       WRITE(lunout,*) "WARNING ! Physics time step changes from a call to the next",pdtphys_,pdtphys
       WRITE(lunout,*) "Not sure the physics parametrizations can handle this..."
    ENDIF
    pdtphys=pdtphys_
    
    ! Update elapsed time since begining of run:
    current_time = current_time + pdtphys
    cur_day = int(current_time/rday)
    cur_sec = current_time - (cur_day * rday) 

    ! Compute corresponding Julian date and update calendar
    cur_day = cur_day + day_ini
    cur_sec = cur_sec + (start_time * rday)
    CALL ymds2ju(annee_ref,1, cur_day, cur_sec, julian_date)
    CALL phys_cal_update(julian_date)
    
  END SUBROUTINE update_time

END MODULE time_phylmdz_mod      

