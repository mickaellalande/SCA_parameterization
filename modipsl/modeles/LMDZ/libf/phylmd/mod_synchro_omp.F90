MODULE mod_synchro_omp

    INTEGER,SAVE :: exit_omp  

CONTAINS

  SUBROUTINE Init_synchro_omp
  USE mod_phys_lmdz_para 
  IMPLICIT NONE

  CALL omp_barrier
!$OMP MASTER
  exit_omp=0
!$OMP END MASTER    
  CALL omp_barrier

  END SUBROUTINE Init_Synchro_omp
  
  SUBROUTINE Synchro_omp  
  USE mod_phys_lmdz_para
  IMPLICIT NONE
  LOGICAL :: out
  
    out=.FALSE.
!$OMP BARRIER
!$OMP BARRIER
!$OMP ATOMIC  
    exit_omp=exit_omp+1
!$OMP BARRIER
!$OMP BARRIER
    IF (exit_omp==omp_size) THEN
      out=.TRUE.
    ENDIF
    
    DO WHILE (.NOT. out)
!$OMP BARRIER
      IF (exit_omp==omp_size) out=.TRUE.
!$OMP BARRIER
    ENDDO

!$OMP BARRIER
!$OMP MASTER
    exit_omp=0
!$OMP END MASTER
!$OMP BARRIER

    IF (exit_omp/=0) THEN
      STOP 'synchro_omp'
    ENDIF

  END SUBROUTINE Synchro_omp


END MODULE mod_synchro_omp
