INTERFACE
SUBROUTINE RRTM_TAUMOL3 (KLEV,P_TAU,&
 & P_TAUAERL,P_FAC00,P_FAC01,P_FAC10,P_FAC11,P_FORFAC,K_JP,K_JT,K_JT1,P_ONEMINUS,&
 & P_COLH2O,P_COLCO2,P_COLN2O,K_LAYTROP,P_SELFFAC,P_SELFFRAC,K_INDSELF,PFRAC) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE PARRRTM , ONLY : JPLAY ,JPBAND ,JPGPT ,NG3 ,NGS2
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB) ,INTENT(OUT) :: P_TAU(JPGPT,JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_TAUAERL(JPLAY,JPBAND)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC00(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC01(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC10(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FAC11(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_FORFAC(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_JP(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_JT(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_JT1(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_ONEMINUS
REAL(KIND=JPRB) ,INTENT(IN) :: P_COLH2O(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_COLCO2(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_COLN2O(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_LAYTROP
REAL(KIND=JPRB) ,INTENT(IN) :: P_SELFFAC(JPLAY)
REAL(KIND=JPRB) ,INTENT(IN) :: P_SELFFRAC(JPLAY)
INTEGER(KIND=JPIM),INTENT(IN) :: K_INDSELF(JPLAY)
REAL(KIND=JPRB) ,INTENT(OUT) :: PFRAC(JPGPT,JPLAY)
END SUBROUTINE RRTM_TAUMOL3
END INTERFACE
