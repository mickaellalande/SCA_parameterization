!
! $Id: YOETHF.h 2799 2017-02-24 18:50:33Z jyg $
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!*    COMMON *YOETHF* DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
!
!     *R__ES*   *CONSTANTS USED FOR COMPUTATION OF SATURATION
!                MIXING RATIO OVER LIQUID WATER(*R_LES*) OR
!                ICE(*R_IES*).
!     *RVTMP2*  *RVTMP2=RCPV/RCPD-1.
!     *RHOH2O*  *DENSITY OF LIQUID WATER.   (RATM/100.)
!
      REAL R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES
      REAL RVTMP2, RHOH2O
      REAL R5ALVCP,R5ALSCP,RALVDCP,RALSDCP,RALFDCP,RTWAT,RTBER,RTBERCU
      REAL RTICE,RTICECU,RTWAT_RTICE_R,RTWAT_RTICECU_R,RKOOP1,RKOOP2
      LOGICAL OK_BAD_ECMWF_THERMO ! If TRUE, then variables set by rrtm/suphec.F90
                                  ! If FALSE, then variables set by suphel.F90
      COMMON /YOETHF/R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES,    &
     &               RVTMP2, RHOH2O,                                    &
     &               R5ALVCP,R5ALSCP,RALVDCP,RALSDCP,                   &
     &               RALFDCP,RTWAT,RTBER,RTBERCU,                       &
     &               RTICE,RTICECU,RTWAT_RTICE_R,RTWAT_RTICECU_R,RKOOP1,&
     &               RKOOP2,                                            &
     &               OK_BAD_ECMWF_THERMO

!$OMP THREADPRIVATE(/YOETHF/)
