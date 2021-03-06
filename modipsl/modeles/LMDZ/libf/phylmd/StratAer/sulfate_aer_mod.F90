MODULE sulfate_aer_mod

! microphysical routines based on UPMC aerosol model by Slimane Bekki
! adapted for stratospheric sulfate aerosol in LMDZ by Christoph Kleinschmitt 

CONTAINS

!********************************************************************
    SUBROUTINE STRACOMP(sh,t_seri,pplay)

!   AEROSOL H2SO4 WEIGHT FRACTION AS A FUNCTION OF PH2O AND TEMPERATURE
!   ----------------------------------------------------------------
!   INPUT: 
!   H2O: VMR of H2O
!   t_seri: temperature (K)
!   PMB: pressure (mb)
!   klon: number of latitude bands in the model domain
!   klev: number of altitude bands in the model domain
!   for IFS: perhaps add another dimension for longitude
!
!   OUTPUT: 
!   R2SO4: aerosol H2SO4 weight fraction (percent)
 
    USE dimphy, ONLY : klon,klev
    USE aerophys
    USE phys_local_var_mod, ONLY: R2SO4

    IMPLICIT NONE

    REAL,DIMENSION(klon,klev),INTENT(IN)          :: t_seri  ! Temperature
    REAL,DIMENSION(klon,klev),INTENT(IN)          :: pplay   ! pression pour le mileu de chaque couche (en Pa)
    REAL,DIMENSION(klon,klev),INTENT(IN)          :: sh      ! humidite specifique
     
    REAL PMB(klon,klev), H2O(klon,klev)
!
!   working variables
    INTEGER I,J,K
    REAL TP, PH2O, VAL, A, B
!     local variables to be saved on exit
    INTEGER INSTEP
    INTEGER, PARAMETER :: N=16, M=28
    DATA INSTEP/0/
    REAL F(N,M)
    REAL XC(N)
    REAL YC(M)
    REAL XC1, XC16, YC1, YC28
!
    SAVE INSTEP,F,XC,YC,XC1,XC16,YC1,YC28

! convert pplay (in Pa) to PMB (in mb)
    PMB(:,:)=pplay(:,:)/100.0

! convert specific humidity sh (in kg/kg) to VMR H2O
    H2O(:,:)=sh(:,:)*mAIRmol/mH2Omol

    IF(INSTEP.EQ.0) THEN
    
       INSTEP=1
       XC(1)=0.01
       XC(2)=0.1
       XC(3)=0.5
       XC(4)=1.0
       XC(5)=1.5
       XC(6)=2.0
       XC(7)=3.0
       XC(8)=5.0
       XC(9)=6.0
       XC(10)=8.0
       XC(11)=10.0
       XC(12)=12.0
       XC(13)=15.0
       XC(14)=20.0
       XC(15)=30.0
       XC(16)=100.0
!
       YC(1)=175.0
       DO I=2,28
         YC(I)=YC(I-1)+5.0
       ENDDO

!      CONVERSION mb IN 1.0E-4mB
       DO I=1,16
         XC(I)=XC(I)*1.0E-4
       ENDDO
!
       XC1=XC(1)+1.E-10
       XC16=XC(16)-1.E-8
       YC1=YC(1)+1.E-5
       YC28=YC(28)-1.E-5

       F(6,4)=43.45
       F(6,5)=53.96
       F(6,6)=60.62
       F(6,7)=65.57
       F(6,8)=69.42
       F(6,9)=72.56
       F(6,10)=75.17
       F(6,11)=77.38
       F(6,12)=79.3
       F(6,13)=80.99
       F(6,14)=82.5
       F(6,15)=83.92
       F(6,16)=85.32
       F(6,17)=86.79
       F(6,18)=88.32
!
!      ADD FACTOR  BECAUSE THE SLOP IS TOO IMPORTANT
!      NOT FOR THIS ONE BUT THE REST
!      LOG DOESN'T WORK
       A=(F(6,5)-F(6,4))/( (YC(5)-YC(4))*2.0)
       B=-A*YC(4) + F(6,4)
       F(6,1)=A*YC(1) + B
       F(6,2)=A*YC(2) + B
       F(6,3)=A*YC(3) + B
!
       F(7,4)=37.02
       F(7,5)=49.46
       F(7,6)=57.51
       F(7,7)=63.12
       F(7,8)=67.42
       F(7,9)=70.85
       F(7,10)=73.70
       F(7,11)=76.09
       F(7,12)=78.15
       F(7,13)=79.96
       F(7,14)=81.56
       F(7,15)=83.02
       F(7,16)=84.43
       F(7,17)=85.85
       F(7,18)=87.33
!
       A=(F(7,5)-F(7,4))/( (YC(5)-YC(4))*2.0)
       B=-A*YC(4) + F(7,4)
       F(7,1)=A*YC(1) + B
       F(7,2)=A*YC(2) + B
       F(7,3)=A*YC(3) + B
!
       F(8,4)=25.85
       F(8,5)=42.26
       F(8,6)=52.78
       F(8,7)=59.55
       F(8,8)=64.55
       F(8,9)=68.45
       F(8,10)=71.63
       F(8,11)=74.29
       F(8,12)=76.56
       F(8,13)=78.53
       F(8,14)=80.27
       F(8,15)=81.83
       F(8,16)=83.27
       F(8,17)=84.67
       F(8,18)=86.10
!
       A=(F(8,5)-F(8,4))/( (YC(5)-YC(4))*2.5 )
       B=-A*YC(4) + F(8,4)
       F(8,1)=A*YC(1) + B
       F(8,2)=A*YC(2) + B
       F(8,3)=A*YC(3) + B
!
       F(9,4)=15.38
       F(9,5)=39.35
       F(9,6)=50.73
       F(9,7)=58.11
       F(9,8)=63.41
       F(9,9)=67.52
       F(9,10)=70.83
       F(9,11)=73.6
       F(9,12)=75.95
       F(9,13)=77.98
       F(9,14)=79.77
       F(9,15)=81.38
       F(9,16)=82.84
       F(9,17)=84.25
       F(9,18)=85.66
!
       A=(F(9,5)-F(9,4))/( (YC(5)-YC(4))*7.0)
       B=-A*YC(4) + F(9,4)
       F(9,1)=A*YC(1) + B
       F(9,2)=A*YC(2) + B
       F(9,3)=A*YC(3) + B
!
       F(10,4)=0.0
       F(10,5)=34.02
       F(10,6)=46.93
       F(10,7)=55.61
       F(10,8)=61.47
       F(10,9)=65.94
       F(10,10)=69.49
       F(10,11)=72.44
       F(10,12)=74.93
       F(10,13)=77.08
       F(10,14)=78.96
       F(10,15)=80.63
       F(10,16)=82.15
       F(10,17)=83.57
       F(10,18)=84.97
!
       A=(F(10,6)-F(10,5))/( (YC(6)-YC(5))*1.5)
       B=-A*YC(5) + F(10,5)
       F(10,1)=A*YC(1) + B
       F(10,2)=A*YC(2) + B
       F(10,3)=A*YC(3) + B
       F(10,4)=A*YC(4) + B
!
       F(11,4)=0.0
       F(11,5)=29.02
       F(11,6)=43.69
       F(11,7)=53.44
       F(11,8)=59.83
       F(11,9)=64.62
       F(11,10)=68.39
       F(11,11)=71.48
       F(11,12)=74.10
       F(11,13)=76.33
       F(11,14)=78.29
       F(11,15)=80.02
       F(11,16)=81.58
       F(11,17)=83.03
       F(11,18)=84.44
!
       A=(F(11,6)-F(11,5))/( (YC(6)-YC(5))*2.5 )
       B=-A*YC(5) + F(11,5)
       F(11,1)=A*YC(1) + B
       F(11,2)=A*YC(2) + B
       F(11,3)=A*YC(3) + B
       F(11,4)=A*YC(4) + B
!
       F(12,4)=0.0
       F(12,5)=23.13
       F(12,6)=40.86
       F(12,7)=51.44
       F(12,8)=58.38
       F(12,9)=63.47
       F(12,10)=67.43
       F(12,11)=70.66
       F(12,12)=73.38
       F(12,13)=75.70
       F(12,14)=77.72
       F(12,15)=79.51
       F(12,16)=81.11
       F(12,17)=82.58
       F(12,18)=83.99
!
       A=(F(12,6)-F(12,5))/( (YC(6)-YC(5))*3.5 )
       B=-A*YC(5) + F(12,5)
       F(12,1)=A*YC(1) + B
       F(12,2)=A*YC(2) + B
       F(12,3)=A*YC(3) + B
       F(12,4)=A*YC(4) + B
!
       F(13,4)=0.0
       F(13,5)=0.0
       F(13,6)=36.89
       F(13,7)=48.63
       F(13,8)=56.46
       F(13,9)=61.96
       F(13,10)=66.19
       F(13,11)=69.6
       F(13,12)=72.45
       F(13,13)=74.89
       F(13,14)=76.99
       F(13,15)=78.85
       F(13,16)=80.50
       F(13,17)=82.02
       F(13,18)=83.44
!
       A=(F(13,7)-F(13,6))/( (YC(7)-YC(6))*2.0)
       B=-A*YC(6) + F(13,6)
       F(13,1)=A*YC(1) + B
       F(13,2)=A*YC(2) + B
       F(13,3)=A*YC(3) + B
       F(13,4)=A*YC(4) + B
       F(13,5)=A*YC(5) + B
!
       F(14,4)=0.0
       F(14,5)=0.0
       F(14,6)=30.82
       F(14,7)=44.49
       F(14,8)=53.69
       F(14,9)=59.83
       F(14,10)=64.47
       F(14,11)=68.15
       F(14,12)=71.19
       F(14,13)=73.77
       F(14,14)=76.0
       F(14,15)=77.95
       F(14,16)=79.69
       F(14,17)=81.26
       F(14,18)=82.72
!
       A=(F(14,7)-F(14,6))/( (YC(7)-YC(6))*2.5 )
       B=-A*YC(6) + F(14,6)
       F(14,1)=A*YC(1) + B
       F(14,2)=A*YC(2) + B
       F(14,3)=A*YC(3) + B
       F(14,4)=A*YC(4) + B
       F(14,5)=A*YC(5) + B
!
       F(15,4)=0.0
       F(15,5)=0.0
       F(15,6)=0.0
       F(15,7)=37.71
       F(15,8)=48.49
       F(15,9)=56.40
       F(15,10)=61.75
       F(15,11)=65.89
       F(15,12)=69.25
       F(15,13)=72.07
       F(15,14)=74.49
       F(15,15)=76.59
       F(15,16)=78.45
       F(15,17)=80.12
       F(15,18)=81.64
!
       A=(F(15,8)-F(15,7))/( (YC(8)-YC(7))*1.5)
       B=-A*YC(7) + F(15,7)
       F(15,1)=A*YC(1) + B
       F(15,2)=A*YC(2) + B
       F(15,3)=A*YC(3) + B
       F(15,4)=A*YC(4) + B
       F(15,5)=A*YC(5) + B
       F(15,6)=A*YC(6) + B

!      SUPPOSE THAT AT GIVEN  AND PH2O<2mB,
!      %H2SO4 = A *LOG(PH2O) +B
!      XC(1-5) :EXTENSION LEFT (LOW H2O)
       DO J=1,18
         A=(F(6,J)-F(7,J))/(LOG(XC(6))-LOG(XC(7)))
         B=-A*LOG(XC(6)) + F(6,J)
         DO K=1,5
           F(K,J)=A*LOG(XC(K)) + B
         ENDDO
       ENDDO

!      XC(16) :EXTENSION RIGHT (HIGH H2O)
       DO J=1,18
         A=(F(15,J)-F(14,J))/(XC(15)-XC(14))
         B=-A*XC(15) + F(15,J)
       F(16,J)=A*XC(16) + B
!       F(16,2)=1.0
       ENDDO

!      YC(16-25) :EXTENSION DOWN (HIGH T)
       DO I=1,16
         A=(F(I,18)-F(I,17))/(YC(18)-YC(17))
         B=-A*YC(18) + F(I,18)
         DO K=19,28
         F(I,K)=A*YC(K) + B
         ENDDO
       ENDDO

!      MANUAL CORRECTIONS
       DO J=1,10
       F(1,J)=94.0
       ENDDO

       DO J=1,6
       F(2,J)=77.0 +REAL(J)
       ENDDO

       DO J=1,7
       F(16,J)=9.0
       ENDDO

       DO I=1,16
       DO J=1,28
         IF (F(I,J).LT.9.0)  F(I,J)=30.0
         IF (F(I,J).GT.99.99) F(I,J)=99.99
       ENDDO
       ENDDO
      
    ENDIF

    DO I=1,klon
    DO J=1,klev
        TP=t_seri(I,J)
        IF (TP.LT.175.1) TP=175.1
!    Partial pressure of H2O (mb)
        PH2O =PMB(I,J)*H2O(I,J)
        IF (PH2O.LT.XC1) THEN 
          R2SO4(I,J)=99.99
!          PH2O=XC(1)+1.0E-10
        ELSE 
          IF (PH2O.GT.XC16) PH2O=XC16 
!         SIMPLE LINEAR INTERPOLATIONS
          CALL FIND(PH2O,TP,XC,YC,F,VAL,N,M)
          IF (PMB(I,J).GE.10.0.AND.VAL.LT.60.0) VAL=60.0
          R2SO4(I,J)=VAL
        ENDIF
    ENDDO
    ENDDO

    END SUBROUTINE

!****************************************************************
    SUBROUTINE STRAACT(ACTSO4)

!   H2SO4 ACTIVITY (GIAUQUE) AS A FUNCTION OF H2SO4 WP
!   ----------------------------------------
!   INPUT: 
!   H2SO4: VMR of H2SO4
!   klon: number of latitude bands in the model domain
!   klev: number of altitude bands in the model domain
!   for IFS: perhaps add another dimension for longitude
!
!   OUTPUT: 
!   ACTSO4: H2SO4 activity (percent)
 
    USE dimphy, ONLY : klon,klev
    USE phys_local_var_mod, ONLY: R2SO4

    IMPLICIT NONE
     
    REAL ACTSO4(klon,klev)
   
!   Working variables         
    INTEGER NN,I,J,JX,JX1
    REAL TC,TB,TA,XT
    PARAMETER (NN=109)
    REAL XC(NN),  X(NN)

!   H2SO4 activity 
    DATA X/ & 
     &   0.0,0.25,0.78,1.437,2.19,3.07,4.03,5.04,6.08 & 
     &  ,7.13,8.18,14.33,18.59,28.59,39.17,49.49 & 
     &  ,102.4,157.8,215.7,276.9,341.6,409.8,481.5,556.6 & 
     &  ,635.5,719.,808.,902.,1000.,1103.,1211.,1322.,1437.,1555. & 
     &  ,1677.,1800.,1926.,2054.,2183.,2312.,2442.,2572.,2701.,2829. & 
     &  ,2955.,3080.,3203.,3325.,3446.,3564.,3681.,3796.,3910.,4022. & 
     &  ,4134.,4351.,4564.,4771.,4974.,5171.,5364.,5551.,5732.,5908. & 
     &  ,6079.,6244.,6404.,6559.,6709.,6854.,6994.,7131.,7264.,7393. & 
     &  ,7520.,7821.,8105.,8373.,8627.,8867.,9093.,9308.,9511.,9703. & 
     &  ,9885.,10060.,10225.,10535.,10819.,11079.,11318.,11537. & 
     &  ,11740.,12097.,12407.,12676.,12915.,13126.,13564.,13910. & 
     &  ,14191.,14423.,14617.,14786.,10568.,15299.,15491.,15654. & 
     &  ,15811./
!   H2SO4 weight fraction (percent)
    DATA XC/ & 
     &   100.0,99.982,99.963,99.945,99.927,99.908,99.890,99.872 & 
     &  ,99.853,99.835,99.817,99.725,99.634,99.452,99.270 & 
     &  ,99.090,98.196,97.319,96.457,95.610,94.777,93.959,93.156 & 
     &  ,92.365,91.588,90.824,90.073,89.334,88.607,87.892,87.188 & 
     &  ,86.495,85.814,85.143,84.482,83.832,83.191,82.560,81.939 & 
     &  ,81.327,80.724,80.130,79.545,78.968,78.399,77.839,77.286 & 
     &  ,76.741,76.204,75.675,75.152,74.637,74.129,73.628,73.133 & 
     &  ,72.164,71.220,70.300,69.404,68.530,67.678,66.847,66.037 & 
     &  ,65.245,64.472,63.718,62.981,62.261,61.557,60.868,60.195 & 
     &  ,59.537,58.893,58.263,57.646,56.159,54.747,53.405,52.126 & 
     &  ,50.908,49.745,48.634,47.572,46.555,45.580,44.646,43.749 & 
     &  ,42.059,40.495,39.043,37.691,36.430,35.251,33.107,31.209 & 
     &  ,29.517,27.999,26.629,23.728,21.397,19.482,17.882,16.525 & 
     &  ,15.360,13.461,11.980,10.792,9.819,8.932/

    DO I=1,klon
    DO J=1,klev
!     HERE LINEAR INTERPOLATIONS
        XT=R2SO4(I,J)
        CALL POSACT(XT,XC,NN,JX)
        JX1=JX+1
        IF(JX.EQ.0) THEN
          ACTSO4(I,J)=0.0 
        ELSE IF(JX.GE.NN) THEN
          ACTSO4(I,J)=15811.0 
        ELSE 
          TC=XT            -XC(JX)
          TB=X(JX1)        -X(JX)
          TA=XC(JX1)       -XC(JX)
          TA=TB/TA
          ACTSO4(I,J)=X(JX)  + TA*TC
        ENDIF
    ENDDO
    ENDDO

    END SUBROUTINE

!****************************************************************
    SUBROUTINE DENH2SA(t_seri)

!   AERSOL DENSITY AS A FUNCTION OF H2SO4 WEIGHT PERCENT AND T
!   ---------------------------------------------
!   VERY ROUGH APPROXIMATION (SEE FOR WATER IN HANDBOOK
!   LINEAR 2% FOR 30 DEGREES with RESPECT TO WATER)
!   
!   INPUT: 
!   R2SO4: aerosol H2SO4 weight fraction (percent)
!   t_seri: temperature (K)
!   klon: number of latitude bands in the model domain
!   klev: number of altitude bands in the model domain
!   for IFS: perhaps add another dimension for longitude
!
!   OUTPUT: 
!   DENSO4: aerosol mass density (gr/cm3 = aerosol mass/aerosol volume)
!   
    USE dimphy, ONLY : klon,klev
    USE phys_local_var_mod, ONLY: R2SO4, DENSO4

    IMPLICIT NONE

    REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
        
    INTEGER I,J

!   Loop on model domain (2 dimension for UPMC model; 3 for IFS)
    DO I=1,klon
    DO J=1,klev
!     RO AT 20C
      DENSO4(I,J)=0.78681252E-5*R2SO4(I,J)*R2SO4(I,J)+ 0.82185978E-2*R2SO4(I,J)+0.97968381
      DENSO4(I,J)=DENSO4(I,J)* ( 1.0 - (t_seri(I,J)-293.0)*0.02/30.0 )
    ENDDO
    ENDDO

    END SUBROUTINE

!***********************************************************
    SUBROUTINE FIND(X,Y,XC,YC,F,VAL,N,M)
!
!   BI-LINEAR INTERPOLATION

!   INPUT: 
!   X: Partial pressure of H2O (mb)
!   Y: temperature (K)
!   XC: Table partial pressure of H2O (mb)
!   YC: Table temperature (K)
!   F: Table aerosol H2SO4 weight fraction=f(XC,YC) (percent)
!
!   OUTPUT: 
!   VAL: aerosol H2SO4 weight fraction (percent)
 
    IMPLICIT NONE
   
    INTEGER N,M
    REAL X,Y,XC(N),YC(M),F(N,M),VAL
!
!   working variables
    INTEGER  IERX,IERY,JX,JY,JXP1,JYP1
    REAL SXY,SX1Y,SX1Y1,SXY1,TA,TB,T,UA,UB,U

    IERX=0
    IERY=0
    CALL POSITION(XC,X,N,JX,IERX)
    CALL POSITION(YC,Y,M,JY,IERY)

    IF(JX.EQ.0.OR.IERY.EQ.1) THEN
       VAL=99.99
       RETURN
    ENDIF

    IF(JY.EQ.0.OR.IERX.EQ.1) THEN
       VAL=9.0
       RETURN
    ENDIF

    JXP1=JX+1
    JYP1=JY+1
    SXY=F(JX,  JY  )
    SX1Y=F(JXP1,JY  )
    SX1Y1=F(JXP1,JYP1)
    SXY1=F(JX,  JYP1)

!   x-slope.
    TA=X       -XC(JX)
    TB=XC(JXP1)-XC(JX)
    T=TA/TB

!   y-slope.
    UA=Y       -YC(JY)
    UB=YC(JYP1)-YC(JY)
    U=UA/UB

!   Use bilinear interpolation to determine function at point X,Y.
    VAL=(1.-T)*(1.-U)*SXY + T*(1.0-U)*SX1Y + T*U*SX1Y1 + (1.0-T)*U*SXY1

    IF(VAL.LT.9.0) VAL=9.0
    IF(VAL.GT.99.99) VAL=99.99
   
    RETURN
    END SUBROUTINE
!****************************************************************
       SUBROUTINE POSITION(XC,X,N,JX,IER)
 
       IMPLICIT NONE
   
       INTEGER N,JX,IER,I
       REAL X,XC(N)

       IER=0
       IF(X.LT.XC(1)) THEN
         JX=0
       ELSE
         DO 10 I=1,N
           IF (X.LT.XC(I)) GO TO 20
 10      CONTINUE
         IER=1
 20      JX=I-1
       ENDIF

       RETURN
       END SUBROUTINE
!********************************************************************
       SUBROUTINE POSACT(XT,X,N,JX)
   
!      POSITION OF XT IN THE ARRAY X
!    -----------------------------------------------
   
       IMPLICIT NONE
   
       INTEGER N
       REAL XT,X(N)
!      Working variables	  	   
       INTEGER JX,I
  
       IF(XT.GT.X(1)) THEN
         JX=0
       ELSE
         DO 10 I=1,N
           IF (XT.GT.X(I)) GO TO 20
 10      CONTINUE
 20      JX=I
       ENDIF
   
       RETURN
       END SUBROUTINE

END MODULE sulfate_aer_mod
