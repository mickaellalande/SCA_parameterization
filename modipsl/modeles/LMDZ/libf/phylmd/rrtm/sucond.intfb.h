INTERFACE
SUBROUTINE SUCOND ( KULOUT , KLEV , PETA )
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KULOUT
REAL(KIND=JPRB) ,INTENT(IN) :: PETA(KLEV)
END SUBROUTINE SUCOND
END INTERFACE
