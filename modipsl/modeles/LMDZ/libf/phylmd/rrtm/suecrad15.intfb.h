INTERFACE
SUBROUTINE SUECRAD15 (KULOUT, KLEV, PETAH )
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KULOUT
REAL(KIND=JPRB) ,INTENT(IN) :: PETAH(KLEV+1)
END SUBROUTINE SUECRAD15
END INTERFACE
