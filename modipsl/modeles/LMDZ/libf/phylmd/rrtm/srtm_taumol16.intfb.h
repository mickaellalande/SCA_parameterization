INTERFACE
SUBROUTINE SRTM_TAUMOL16&
 & ( KLEV,&
 & P_FAC00 , P_FAC01 , P_FAC10 , P_FAC11,&
 & K_JP , K_JT , K_JT1 , P_ONEMINUS,&
 & P_COLH2O , P_COLCH4 , P_COLMOL,&
 & K_LAYTROP , P_SELFFAC , P_SELFFRAC, K_INDSELF , P_FORFAC , P_FORFRAC, K_INDFOR,&
 & P_SFLUXZEN, P_TAUG , P_TAUR&
 & ) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE PARSRTM , ONLY : JPLAY, JPG, NG16
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC00(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC01(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC10(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC11(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_JP(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_JT(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_JT1(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_ONEMINUS
REAL(KIND=JPRB) ,INTENT(IN) :: P_COLH2O(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_COLCH4(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_COLMOL(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_LAYTROP
REAL(KIND=JPRB) ,INTENT(IN) :: P_SELFFAC(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_SELFFRAC(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_INDSELF(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FORFAC(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FORFRAC(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_INDFOR(JPLAY)
REAL(KIND=JPRB) ,INTENT(OUT) :: P_SFLUXZEN(JPG)
REAL(KIND=JPRB) ,INTENT(OUT) :: P_TAUG(JPLAY,JPG)
REAL(KIND=JPRB) ,INTENT(OUT) :: P_TAUR(JPLAY,JPG)
END SUBROUTINE SRTM_TAUMOL16
END INTERFACE
