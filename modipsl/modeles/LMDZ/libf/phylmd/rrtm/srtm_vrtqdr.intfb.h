INTERFACE
SUBROUTINE SRTM_VRTQDR&
 & ( KLEV , KW,&
 & PREF , PREFD, PTRA , PTRAD,&
 & PDBT , PRDND, PRUP , PRUPD , PTDBT,&
 & PFD , PFU&
 & ) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE PARSRTM , ONLY : JPLAY, JPGPT
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KW
REAL(KIND=JPRB) ,INTENT(IN) :: PREF(JPLAY+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PREFD(JPLAY+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PTRA(JPLAY+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PTRAD(JPLAY+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PDBT(JPLAY+1)
REAL(KIND=JPRB) ,INTENT(OUT) :: PRDND(JPLAY+1)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PRUP(JPLAY+1)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PRUPD(JPLAY+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PTDBT(JPLAY+1)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFD(JPLAY+1,JPGPT)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFU(JPLAY+1,JPGPT)
END SUBROUTINE SRTM_VRTQDR
END INTERFACE
