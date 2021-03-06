!
! $Id: srtm_srtm_224gp.F90 2027 2014-04-29 13:38:53Z fairhead $
!
SUBROUTINE SRTM_SRTM_224GP &
 & ( KIDIA , KFDIA  , KLON  , KLEV  , KSW , KOVLP ,&
 &   PAER  , PALBD  , PALBP , PAPH  , PAP ,&
 &   PTS   , PTH    , PT    ,&
 &   PQ    , PCCO2  , POZN  , PRMU0 ,&
 &   PFRCL , PTAUC  , PASYC , POMGC ,&
 &   PALBT , PFSUX  , PFSUC &
 & )  

!-- interface to RRTM_SW
!     JJMorcrette 030225

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARSRTM  , ONLY : JPLAY
!USE YOERDI   , ONLY : RCH4   , RN2O   
USE YOERAD   , ONLY : NAER
USE YOESRTAER, ONLY : RSRTAUA, RSRPIZA, RSRASYA
USE YOMPHY3  , ONLY : RII0
USE YOMCST   , ONLY : RI0 



IMPLICIT NONE

#include "clesphys.h"

!-- Input arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM)               :: KLEV! UNDETERMINED INTENT 
INTEGER(KIND=JPIM)               :: KSW! UNDETERMINED INTENT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KOVLP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV)    ! top to bottom
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCO2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRCL(KLON,KLEV)     ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUC(KLON,KSW,KLEV) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PASYC(KLON,KSW,KLEV) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMGC(KLON,KSW,KLEV) ! bottom to top
REAL(KIND=JPRB)                  :: PALBT(KLON,KSW) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSUX(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSUC(KLON,2,KLEV+1) 
!INTEGER_M :: KMOL, KCLDATM, KNFLAG, KCEFLAG, KIQFLAG, KSTR  

!-- Output arguments

!-----------------------------------------------------------------------

!-- dummy integers

INTEGER(KIND=JPIM) :: ICLDATM, INFLAG, ICEFLAG, I_LIQFLAG, I_NMOL, I_NSTR

INTEGER(KIND=JPIM) :: IK, IMOL, J1, J2, JAE, JL, JK, JSW

!-- dummy reals

REAL(KIND=JPRB) :: Z_PZ(0:JPLAY)   , Z_TZ(0:JPLAY)   , Z_PAVEL(JPLAY)  , Z_TAVEL(JPLAY)
REAL(KIND=JPRB) :: Z_COLDRY(JPLAY) , Z_COLMOL(JPLAY) , Z_WKL(35,JPLAY)
REAL(KIND=JPRB) :: Z_CO2MULT(JPLAY), Z_COLCH4(JPLAY) , Z_COLCO2(JPLAY) , Z_COLH2O(JPLAY)
REAL(KIND=JPRB) :: Z_COLN2O(JPLAY) , Z_COLO2(JPLAY)  , Z_COLO3(JPLAY)
REAL(KIND=JPRB) :: Z_FORFAC(JPLAY) , Z_FORFRAC(JPLAY), Z_SELFFAC(JPLAY), Z_SELFFRAC(JPLAY)
REAL(KIND=JPRB) :: Z_FAC00(JPLAY)  , Z_FAC01(JPLAY)  , Z_FAC10(JPLAY)  , Z_FAC11(JPLAY)
REAL(KIND=JPRB) :: Z_TBOUND        , Z_ONEMINUS    , ZRMU0 , ZADJI0
REAL(KIND=JPRB) :: ZALBD(KSW)    , ZALBP(KSW)    , ZFRCL(JPLAY)
REAL(KIND=JPRB) :: ZTAUC(JPLAY,KSW), ZASYC(JPLAY,KSW), ZOMGC(JPLAY,KSW)
REAL(KIND=JPRB) :: ZTAUA(JPLAY,KSW), ZASYA(JPLAY,KSW), ZOMGA(JPLAY,KSW)

REAL(KIND=JPRB) :: ZBBCD(JPLAY+1), ZBBCU(JPLAY+1), ZBBFD(JPLAY+1), ZBBFU(JPLAY+1)
REAL(KIND=JPRB) :: ZUVCD(JPLAY+1), ZUVCU(JPLAY+1), ZUVFD(JPLAY+1), ZUVFU(JPLAY+1)
REAL(KIND=JPRB) :: ZVSCD(JPLAY+1), ZVSCU(JPLAY+1), ZVSFD(JPLAY+1), ZVSFU(JPLAY+1)
REAL(KIND=JPRB) :: ZNICD(JPLAY+1), ZNICU(JPLAY+1), ZNIFD(JPLAY+1), ZNIFU(JPLAY+1)

INTEGER(KIND=JPIM) :: I_LAYTROP, I_LAYSWTCH, I_LAYLOW
INTEGER(KIND=JPIM) :: INDFOR(JPLAY), INDSELF(JPLAY)
INTEGER(KIND=JPIM) :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

REAL(KIND=JPRB) :: Z_AMD                  ! Effective molecular weight of dry air (g/mol)
REAL(KIND=JPRB) :: Z_AMW                  ! Molecular weight of water vapor (g/mol)
REAL(KIND=JPRB) :: Z_AMCO2                ! Molecular weight of carbon dioxide (g/mol)
REAL(KIND=JPRB) :: Z_AMO                  ! Molecular weight of ozone (g/mol)
REAL(KIND=JPRB) :: Z_AMCH4                ! Molecular weight of methane (g/mol)
REAL(KIND=JPRB) :: Z_AMN2O                ! Molecular weight of nitrous oxide (g/mol)
REAL(KIND=JPRB) :: Z_AMC11                ! Molecular weight of CFC11 (g/mol) - CFCL3
REAL(KIND=JPRB) :: Z_AMC12                ! Molecular weight of CFC12 (g/mol) - CF2CL2
REAL(KIND=JPRB) :: Z_AVGDRO               ! Avogadro's number (molecules/mole)
REAL(KIND=JPRB) :: Z_GRAVIT               ! Gravitational acceleration (cm/sec2)
REAL(KIND=JPRB) :: Z_AMM

! Atomic weights for conversion from mass to volume mixing ratios; these
!  are the same values used in ECRT to assure accurate conversion to vmr
data Z_AMD   /  28.970_JPRB    /
data Z_AMW   /  18.0154_JPRB   /
data Z_AMCO2 /  44.011_JPRB    /
data Z_AMO   /  47.9982_JPRB   /
data Z_AMCH4 /  16.043_JPRB    /
data Z_AMN2O /  44.013_JPRB    /
data Z_AMC11 / 137.3686_JPRB   /
data Z_AMC12 / 120.9140_JPRB   /
data Z_AVGDRO/ 6.02214E23_JPRB /
data Z_GRAVIT/ 9.80665E02_JPRB /

REAL(KIND=JPRB) :: ZCLEAR, ZCLOUD, ZEPSEC, ZTOTCC

INTEGER(KIND=JPIM) :: IOVLP
REAL(KIND=JPRB) :: ZHOOK_HANDLE


#include "srtm_setcoef.intfb.h"
#include "srtm_spcvrt.intfb.h"


!-----------------------------------------------------------------------
!-- calculate information needed ny the radiative transfer routine 

IF (LHOOK) CALL DR_HOOK('SRTM_SRTM_224GP',0,ZHOOK_HANDLE)
ZEPSEC  = 1.E-06_JPRB
Z_ONEMINUS=1.0_JPRB -  ZEPSEC
ZADJI0 = RII0 / RI0
!-- overlap: 1=max-ran, 2=maximum, 3=random
IOVLP=3

!print *,'Entering srtm_srtm_224gp'

ICLDATM = 1
INFLAG    = 2
ICEFLAG    = 3
I_LIQFLAG = 1
I_NMOL    = 6
I_NSTR    = 2

DO JL = KIDIA, KFDIA
  ZRMU0=PRMU0(JL)
  IF (ZRMU0 > 0.0_JPRB) THEN

!- coefficients related to the cloud optical properties (original RRTM_SW)

!  print *,'just before SRTM_CLDPROP'

!  DO JK=1,KLEV 
!    CLDFRAC(JK) = PFRCL (JL,JK)
!    CLDDAT1(JK) = PSCLA1(JL,JK)
!    CLDDAT2(JK) = PSCLA2(JL,JK)
!    CLDDAT3(JK) = PSCLA3(JL,JK)
!    CLDDAT4(JK) = PSCLA4(JL,JK)
!    DO JMOM=0,16
!      CLDDATMOM(JMOM,JK)=PSCLMOM(JL,JMOM,JK)
!    ENDDO
!    print 9101,JK,CLDFRAC(JK),CLDDAT1(JK),CLDDAT2(JK),CLDDAT3(JK)&
!    &,CLDDAT4(JK),(CLDDATMOM(JMOM,JK),JMOM=0,NSTR)
    9101 format(1x,'srtm_srtm_224gp Cld :',I3,f7.4,7E12.5)
!  ENDDO

!  CALL SRTM_CLDPROP &
!    &( KLEV, ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NSTR &
!    &, CLDFRAC, CLDDAT1, CLDDAT2, CLDDAT3, CLDDAT4, CLDDATMOM &
!    &, TAUCLDORIG, TAUCLOUD, SSACLOUD, XMOM &
!    &)

!- coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients

    DO J1=1,35
      DO J2=1,KLEV
        Z_WKL(J1,J2)=0.0_JPRB 
      ENDDO
    ENDDO

    Z_TBOUND=PTS(JL)
    Z_PZ(0) = paph(JL,klev+1)/100._JPRB
    Z_TZ(0) = pth (JL,klev+1)

    ZCLEAR=1.0_JPRB
    ZCLOUD=0.0_JPRB
    ZTOTCC=0.0_JPRB
    DO JK = 1, KLEV
      Z_PAVEL(JK) = pap(JL,KLEV-JK+1) /100._JPRB
      Z_TAVEL(JK) = pt (JL,KLEV-JK+1)
      Z_PZ(JK)    = paph(JL,KLEV-JK+1)/100._JPRB
      Z_TZ(JK)    = pth (JL,KLEV-JK+1)
      Z_WKL(1,JK) = pq(JL,KLEV-JK+1)  *Z_AMD/Z_AMW
      Z_WKL(2,JK) = pcco2             *Z_AMD/Z_AMCO2
      Z_WKL(3,JK) = pozn(JL,KLEV-JK+1)*Z_AMD/Z_AMO
      Z_WKL(4,JK) = rn2o              *Z_AMD/Z_AMN2O
      Z_WKL(6,JK) = rch4              *Z_AMD/Z_AMCH4
      Z_AMM = (1-Z_WKL(1,JK))*Z_AMD + Z_WKL(1,JK)*Z_AMW
      Z_COLDRY(JK) = (Z_PZ(JK-1)-Z_PZ(JK))*1.E3_JPRB*Z_AVGDRO/(Z_GRAVIT*Z_AMM*(1+Z_WKL(1,JK)))
!    print 9200,JK,PAVEL(JK),TAVEL(JK),(WKL(JA,JK),JA=1,4),WKL(6,JK),COLDRY(JK)
      9200 format(1x,'SRTM ',I3,2F7.1,6E13.5)

      IF (KOVLP == 1) THEN
        ZCLEAR=ZCLEAR*(1.0_JPRB-MAX(PFRCL(JL,JK),ZCLOUD)) &
         & /(1.0_JPRB-MIN(ZCLOUD,1.0_JPRB-ZEPSEC))  
        ZCLOUD=PFRCL(JL,JK)
        ZTOTCC=1.0_JPRB-ZCLEAR
      ELSEIF (KOVLP == 2) THEN
        ZCLOUD=MAX(ZCLOUD,PFRCL(JL,JK))
        ZCLEAR=1.0_JPRB-ZCLOUD
        ZTOTCC=ZCLOUD
      ELSEIF (KOVLP == 3) THEN
        ZCLEAR=ZCLEAR*(1.0_JPRB-PFRCL(JL,JK))
        ZCLOUD=1.0_JPRB-ZCLEAR
        ZTOTCC=ZCLOUD
      ENDIF

    ENDDO

!  print *,'ZTOTCC ZCLEAR : ',ZTOTCC,' ',ZCLEAR

    DO IMOL=1,I_NMOL
      DO JK=1,KLEV
        Z_WKL(IMOL,JK)=Z_COLDRY(JK)* Z_WKL(IMOL,JK)
      ENDDO
    ENDDO

!    IF (ZTOTCC == 0.0_JPRB) THEN
!      DO JK=1,KLEV
!        ZFRCL(JK)=0.0_JPRB   
!      ENDDO
!    ELSE
!      DO JK=1,KLEV
!        ZFRCL(JK)=PFRCL(JL,JK)/ZTOTCC
!      ENDDO
!    ENDIF

!  print *,'just before SRTM_SETCOEF'

    ZFRCL(1:KLEV)=PFRCL(JL,1:KLEV)
    ZCLEAR=0._JPRB
    ZCLOUD=1._JPRB

    CALL SRTM_SETCOEF &
     & ( KLEV   , I_NMOL,&
     & Z_PAVEL  , Z_TAVEL   , Z_PZ     , Z_TZ     , Z_TBOUND,&
     & Z_COLDRY , Z_WKL,&
     & I_LAYTROP, I_LAYSWTCH, I_LAYLOW,&
     & Z_CO2MULT, Z_COLCH4  , Z_COLCO2 , Z_COLH2O , Z_COLMOL  , Z_COLN2O  , Z_COLO2 , Z_COLO3,&
     & Z_FORFAC , Z_FORFRAC , INDFOR , Z_SELFFAC, Z_SELFFRAC, INDSELF,&
     & Z_FAC00  , Z_FAC01   , Z_FAC10  , Z_FAC11,&
     & JP     , JT      , JT1     &
     & )  
  
!  print *,'just after SRTM_SETCOEF'

!- call the radiation transfer routine
  
    DO JSW=1,KSW
      ZALBD(JSW)=PALBD(JL,JSW)
      ZALBP(JSW)=PALBP(JL,JSW)
      DO JK=1,KLEV
        ZTAUC(JK,JSW) = PTAUC(JL,JSW,JK)
        ZASYC(JK,JSW) = PASYC(JL,JSW,JK)
        ZOMGC(JK,JSW) = POMGC(JL,JSW,JK)
!      print 9002,JSW,JK,ZFRCL(JK),ZTAUC(JK,JSW),ZASYC(JK,JSW),ZOMGC(JK,JSW)
        9002  format(1x,'srtm_224gp ClOPropECmodel ',2I3,f8.4,3E12.5)
      ENDDO
    ENDDO

!- mixing of aerosols
 
!  print *,'Aerosol optical properties computations'
!  DO JSW=1,KSW
!    print 9012,JSW,(JAE,RSRTAUA(JSW,JAE),RSRPIZA(JSW,JAE),RSRASYA(JSW,JAE),JAE=1,6)
    9012 format(I3,(/,I3,3E13.5))
!  ENDDO

!  DO JK=1,KLEV
!    print 9013,JK,(PAER(JL,JAE,JK),JAE=1,6)
    9013 format(1x,I3,6E12.5)
!  ENDDO

    IF (NAER == 0) THEN
      DO JSW=1,KSW
        DO JK=1,KLEV
          ZTAUA(JK,JSW)= 0.0_JPRB
          ZASYA(JK,JSW)= 0.0_JPRB
          ZOMGA(JK,JSW)= 1.0_JPRB
        ENDDO
      ENDDO
    ELSE
      DO JSW=1,KSW
        DO JK=1,KLEV
          IK=KLEV+1-JK
          ZTAUA(JK,JSW)=0.0_JPRB
          ZASYA(JK,JSW)=0.0_JPRB
          ZOMGA(JK,JSW)=0.0_JPRB
          DO JAE=1,6
            ZTAUA(JK,JSW)=ZTAUA(JK,JSW)+RSRTAUA(JSW,JAE)*PAER(JL,JAE,IK)
            ZOMGA(JK,JSW)=ZOMGA(JK,JSW)+RSRTAUA(JSW,JAE)*PAER(JL,JAE,IK) &
             & *RSRPIZA(JSW,JAE)  
            ZASYA(JK,JSW)=ZASYA(JK,JSW)+RSRTAUA(JSW,JAE)*PAER(JL,JAE,IK) &
             & *RSRPIZA(JSW,JAE)*RSRASYA(JSW,JAE)  
          ENDDO
          IF (ZOMGA(JK,JSW) /= 0.0_JPRB) THEN
            ZASYA(JK,JSW)=ZASYA(JK,JSW)/ZOMGA(JK,JSW)
          ENDIF
          IF (ZTAUA(JK,JSW) /= 0.0_JPRB) THEN
            ZOMGA(JK,JSW)=ZOMGA(JK,JSW)/ZTAUA(JK,JSW)
          ENDIF
!      print 9003,JSW,JK,ZTAUA(JK,JSW),ZOMGA(JK,JSW),ZASYA(JK,JSW)
9003  format(1x,'Aerosols ',2I3,3F10.4)
        ENDDO
      ENDDO
    ENDIF

    DO JK=1,KLEV+1
      ZBBCU(JK)=0.0_JPRB
      ZBBCD(JK)=0.0_JPRB
      ZBBFU(JK)=0.0_JPRB
      ZBBFD(JK)=0.0_JPRB
      ZUVCU(JK)=0.0_JPRB
      ZUVCD(JK)=0.0_JPRB
      ZUVFU(JK)=0.0_JPRB
      ZUVFD(JK)=0.0_JPRB
      ZVSCU(JK)=0.0_JPRB
      ZVSCD(JK)=0.0_JPRB
      ZVSFU(JK)=0.0_JPRB
      ZVSFD(JK)=0.0_JPRB
      ZNICU(JK)=0.0_JPRB
      ZNICD(JK)=0.0_JPRB
      ZNIFU(JK)=0.0_JPRB
      ZNIFD(JK)=0.0_JPRB
    ENDDO

!  print *,'just before calling STRM_SPCVRT for JL=',JL,' and ZRMU0=',ZRMU0

    CALL SRTM_SPCVRT &
     & ( KLEV   , I_NMOL    , KSW    , Z_ONEMINUS,&
     & Z_PAVEL  , Z_TAVEL   , Z_PZ     , Z_TZ     , Z_TBOUND  , ZALBD   , ZALBP,&
     & ZFRCL  , ZTAUC   , ZASYC  , ZOMGC  , ZTAUA   , ZASYA   , ZOMGA , ZRMU0,&
     & Z_COLDRY , Z_WKL,&
     & I_LAYTROP, I_LAYSWTCH, I_LAYLOW,&
     & Z_CO2MULT, Z_COLCH4  , Z_COLCO2 , Z_COLH2O , Z_COLMOL  , Z_COLN2O  , Z_COLO2 , Z_COLO3,&
     & Z_FORFAC , Z_FORFRAC , INDFOR , Z_SELFFAC, Z_SELFFRAC, INDSELF,&
     & Z_FAC00  , Z_FAC01   , Z_FAC10  , Z_FAC11,&
     & JP     , JT      , JT1,&
     & ZBBFD  , ZBBFU   , ZUVFD  , ZUVFU  , ZVSFD   , ZVSFU   , ZNIFD , ZNIFU,&
     & ZBBCD  , ZBBCU   , ZUVCD  , ZUVCU  , ZVSCD   , ZVSCU   , ZNICD , ZNICU &
     & )  

!  print *,'SRTM_SRTM_224GP before potential scaling'
!    IF (IOVLP == 3) THEN
!      DO JK=1,KLEV+1
!!      print 9004,JK,ZBBCU(JK),ZBBCD(JK),ZBBFU(JK),ZBBFD(JK)
        9004 format(1x,'Clear-sky and total fluxes U & D ',I3,4F10.3)
!        PFSUC(JL,1,JK)=ZBBCU(JK)
!        PFSUC(JL,2,JK)=ZBBCD(JK)
!        PFSUX(JL,1,JK)=ZBBFU(JK)
!        PFSUX(JL,2,JK)=ZBBFD(JK)
!      ENDDO
!    ELSE
!    print *,'SRTM_SRTM_224GP after potential scaling'
      DO JK=1,KLEV+1
        PFSUC(JL,1,JK)=ZADJI0 * ZBBCU(JK)
        PFSUC(JL,2,JK)=ZADJI0 * ZBBCD(JK)
        PFSUX(JL,1,JK)=ZADJI0 * ( (1.0_JPRB-ZCLEAR)*ZBBFU(JK)+ZCLEAR*ZBBCU(JK) )
        PFSUX(JL,2,JK)=ZADJI0 * ( (1.0_JPRB-ZCLEAR)*ZBBFD(JK)+ZCLEAR*ZBBCD(JK) )
      ENDDO
!    ENDIF

!  DO JK=1,KLEV+1
!    print 9005,JK,PFSUC(JL,1,JK),PFSUC(JL,2,JK),PFSUX(JL,1,JK),PFSUX(JL,2,JK)
    9005 format(1x,'Clear-sky and total fluxes U & D ',I3,4F10.3)
!  ENDDO
  
  ELSE
    DO JK=1,KLEV+1
      PFSUC(JL,1,JK)=0.0_JPRB
      PFSUC(JL,2,JK)=0.0_JPRB
      PFSUX(JL,1,JK)=0.0_JPRB
      PFSUX(JL,2,JK)=0.0_JPRB
    ENDDO
  ENDIF
ENDDO

!PRINT *,'OUT OF SRTM_224GP'

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_SRTM_224GP',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_SRTM_224GP

