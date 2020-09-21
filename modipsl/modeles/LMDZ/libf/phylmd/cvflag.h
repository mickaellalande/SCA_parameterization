!
! $Header$
!
      logical cvflag_grav
      logical cvflag_ice

      COMMON /cvflag/ cvflag_grav, cvflag_ice  
!$OMP THREADPRIVATE(/cvflag/)
