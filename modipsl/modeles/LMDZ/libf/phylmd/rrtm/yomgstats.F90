MODULE YOMGSTATS

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
! Module for timing statistics. Module is internal to the GSTATS package -
! routines GSTATS, SUSTATS and STATS_OUTPUT. The logical switches are
! re-initialized in SUMPINI

! LSTATS - TRUE for gathering timing statistics
! LSTATSCPU - TRUE for gathering CPU timing  statistics
! LSYNCSTATS - TRUE for syncronization (call to barrier) at the 
!              start of timing event
! LDETAILED_STATS - TRUE for more detail in output
! LSTATS_OMP - TRUE for gathering timing statistics on OpenMP regions
!                 1001-1999
! LSTATS_COMMS - TRUE for gathering detailed timing of Message passing
!                 501-1000
! NTRACE_STATS    - max number of entries in trace
! LTRACE_STATS    - True for trace of all calls to gstats
! LGSTATS_LABEL   - True after GSTATS-labels have been set
! JPMAXSTAT - max number of separate  timers in gstats
! JPOBCOUNT_BASE - first counter for obs types
! NCALLS - number of times a timer has been switched on
! TIMESUM - total time spent with timer on
! TIMESQSUM - sum of the squares of times
! TIMEMAX - max time of all calls
! TIMESUMB - sum of times between previous timer was invoked and this
!            timer was switched on ( to be used for finding out which parts
!            of the code that is not being timed)
! TIMELCALL - time when event was switched on or resumed
! TTCPUSUM - total cpu time
! TVCPUSUM - total vector cpu time
! THISTIME - total accumulated time for this call to timing event (necessary
!            to be able to suspend and resume timer and still have it counted
!            as one timing event)
! THISTCPU - as THISTIME but for CPU time
! THISVCPU - as THISTIME but for vector CPU time
! TTCPULCALL - as TIMELCALL but for CPU time
! TVCPULCALL - as TIMELCALL but for vector CPU time
! TIME_LAST_CALL - last time GSTATS was called
! TIME_START - used for recording parallel startup time


LOGICAL :: LSTATS = .TRUE.
LOGICAL :: LSTATS_OMP = .FALSE.
LOGICAL :: LSTATS_COMMS = .FALSE.
LOGICAL :: LSTATS_MEM = .FALSE.
LOGICAL :: LSTATS_ALLOC = .FALSE.
LOGICAL :: LSTATSCPU = .TRUE.
LOGICAL :: LSYNCSTATS = .FALSE.
LOGICAL :: LDETAILED_STATS = .TRUE.
LOGICAL :: LBARRIER_STATS = .FALSE.
LOGICAL :: LTRACE_STATS = .FALSE.
LOGICAL :: LGSTATS_LABEL = .FALSE.

INTEGER(KIND=JPIM),PARAMETER :: JPMAXSTAT=2500

INTEGER(KIND=JPIM),PARAMETER :: JPOBCOUNT_BASE=201
INTEGER(KIND=JPIM) :: NTRACE_STATS=0
INTEGER(KIND=JPIM) :: NCALLS(0:JPMAXSTAT)
INTEGER(KIND=JPIM) :: NCALLS_TOTAL=0
INTEGER(KIND=JPIM),ALLOCATABLE :: NCALL_TRACE(:)

REAL(KIND=JPRB) :: TIMESUM(0:JPMAXSTAT)
REAL(KIND=JPRB) :: TIMESQSUM(0:JPMAXSTAT)
REAL(KIND=JPRB) :: TIMEMAX(0:JPMAXSTAT)
REAL(KIND=JPRB) :: TIMESUMB(0:JPMAXSTAT)
REAL(KIND=JPRB) :: TIMELCALL(0:JPMAXSTAT)
REAL(KIND=JPRB) :: TTCPUSUM(0:JPMAXSTAT)
REAL(KIND=JPRB) :: TVCPUSUM(0:JPMAXSTAT)
REAL(KIND=JPRB) :: THISTIME(0:JPMAXSTAT)
REAL(KIND=JPRB) :: THISTCPU(0:JPMAXSTAT)
REAL(KIND=JPRB) :: THISVCPU(0:JPMAXSTAT)
REAL(KIND=JPRB) :: TTCPULCALL(0:JPMAXSTAT)
REAL(KIND=JPRB) :: TVCPULCALL(0:JPMAXSTAT)
REAL(KIND=JPRB) :: TIME_LAST_CALL

REAL(KIND=JPRB),ALLOCATABLE :: TIME_START(:)
REAL(KIND=JPRB),ALLOCATABLE :: TIME_TRACE(:)
INTEGER(KIND=JPIM),PARAMETER :: JPERR=0
INTEGER(KIND=JPIM),PARAMETER :: JPTAGSTAT=20555

CHARACTER*50 :: CCDESC(0:JPMAXSTAT) = ""
CHARACTER*3  :: CCTYPE(0:JPMAXSTAT) = ""

INTEGER(KIND=JPIM) :: NPROC_STATS = 1
INTEGER(KIND=JPIM) :: MYPROC_STATS = 1
INTEGER(KIND=JPIM),ALLOCATABLE :: NPRCIDS_STATS(:)

INTEGER(KIND=JPIM) :: NTMEM(0:JPMAXSTAT,5)
INTEGER(KIND=JPIM) :: NSTATS_MEM=0

INTEGER(KIND=JPIM) :: NPRNT_STATS=3

!$OMP THREADPRIVATE(ccdesc,cctype,lbarrier_stats,ldetailed_stats,lgstats_label,lstats,lstats_alloc)
!$OMP THREADPRIVATE(lstats_comms,lstats_mem,lstats_omp,lstatscpu,lsyncstats,ltrace_stats,myproc_stats)
!$OMP THREADPRIVATE(ncalls,ncalls_total,nprnt_stats,nproc_stats,nstats_mem,ntmem,ntrace_stats,thistcpu)
!$OMP THREADPRIVATE(thistime,thisvcpu,time_last_call,timelcall,timemax,timesqsum,timesum,timesumb)
!$OMP THREADPRIVATE(ttcpulcall,ttcpusum,tvcpulcall,tvcpusum)
!$OMP THREADPRIVATE(ncall_trace,nprcids_stats,time_start,time_trace)
END MODULE YOMGSTATS




