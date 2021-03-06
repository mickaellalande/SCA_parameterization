!***************************************************************************
SUBROUTINE RRTM_CMBGB12
!***************************************************************************

!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOERRTO12, ONLY : KAO     ,SELFREFO   ,FRACREFAO
USE YOERRTA12, ONLY : KA      ,SELFREF    ,FRACREFA
USE YOERRTRWT, ONLY : FREFA    ,FREFADF  ,RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JN, JP, JT

REAL(KIND=JPRB) :: Z_SUMF, Z_SUMK
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB12',0,ZHOOK_HANDLE)
DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(12)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(11)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+176)
        ENDDO

        KA(JN,JT,JP,IGC) = Z_SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(12)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(11)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+176)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(12)
    Z_SUMF = 0.0_JPRB
    DO IPR = 1, NGN(NGS(11)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMF = Z_SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = Z_SUMF
  ENDDO
ENDDO

DO JP = 1,9
  DO IGC = 1,NGC(12)

    FREFA(NGS(11)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,8
  DO IGC = 1,NGC(12)

    FREFADF(NGS(11)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB12',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB12
