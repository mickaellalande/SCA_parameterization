MODULE YOECLOP

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!     * YOECLOP* PARAMETERS FOR CLOUD OPTICAL PROPERTIES
!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: RYFWCA(4)
REAL(KIND=JPRB) :: RYFWCB(4)
REAL(KIND=JPRB) :: RYFWCC(4)
REAL(KIND=JPRB) :: RYFWCD(4)
REAL(KIND=JPRB) :: RYFWCE(4)
REAL(KIND=JPRB) :: RYFWCF(4)
REAL(KIND=JPRB) :: REBCUA(4)
REAL(KIND=JPRB) :: REBCUB(4)
REAL(KIND=JPRB) :: REBCUC(4)
REAL(KIND=JPRB) :: REBCUD(4)
REAL(KIND=JPRB) :: REBCUE(4)
REAL(KIND=JPRB) :: REBCUF(4)
REAL(KIND=JPRB) :: RASWCA(4)
REAL(KIND=JPRB) :: RASWCB(4)
REAL(KIND=JPRB) :: RASWCC(4)
REAL(KIND=JPRB) :: RASWCD(4)
REAL(KIND=JPRB) :: RASWCE(4)
REAL(KIND=JPRB) :: RASWCF(4)
REAL(KIND=JPRB) :: REBCUG
REAL(KIND=JPRB) :: REBCUH
REAL(KIND=JPRB) :: REBCUI(6)
REAL(KIND=JPRB) :: REBCUJ(6)
REAL(KIND=JPRB) :: REFFIA
REAL(KIND=JPRB) :: REFFIB
REAL(KIND=JPRB) :: RTIW
REAL(KIND=JPRB) :: RRIW

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
!*    FOUQUART (1987) WATER CLOUD OPTICAL PROPERTIES

! RYFWCA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! RYFWCB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
! RYFWCC :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RYFWCD :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RYFWCE :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RYFWCF :  REAL   : ASSYMETRY FACTOR

!*    SLINGO (1989) WATER CLOUD OPTICAL PROPERTIES

! RASWCA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! RASWCB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
! RASWCC :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RASWCD :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RASWCE :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RASWCF :  REAL   : ASSYMETRY FACTOR

!*    ICE CLOUD OPTICAL PROPERTIES DERIVED FROM EBERT-CURRY (1992)

! REBCUA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! REBCUB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
! REBCUC :  REAL   : 1-C3  IN SINGLE SCATTERING ALBEDO FORMULA
! REBCUD :  REAL   : C4 IN SINGLE SCATTERING ALBEDO FORMULA
! REBCUE :  REAL   : C5 IN ASSYMETRY FACTOR FORMULA
! REBCUF :  REAL   : C6 IN ASSYMETRY FACTOR FORMULA
! REBCUG :  REAL   : C7 IN MASS ABSORPTION COEFFICIENT FORMULA
! REBCUH :  REAL   : C8 IN MASS ABSORPTION COEFFICIENT FORMULA
! REBCUI :  REAL   : C7 IN LW MASS ABSORPTION COEFFICIENT FORMULA
! REBCUJ :  REAL   : C8 IN LW MASS ABSORPTION COEFFICIENT FORMULA

! REFFIA :  REAL   : C9  IN EFFECTIVE RADIUS FORMULA
! REFFIB :  REAL   : C10 IN EFFECTIVE RADIUS FORMULA

!*    TRANSITION BETWEEN LIQUID AND SOLID WATER

! RTIW   :  REAL   : TEMPERATURE THRESHOLD
! RRIW   :  REAL   : TRANSITION RANGE
!     -----------------------------------------------------------------


!$OMP THREADPRIVATE(raswca,raswcb,raswcc,raswcd,raswce,raswcf,rebcua)
!$OMP THREADPRIVATE(rebcub,rebcuc,rebcud,rebcue,rebcuf,rebcug,rebcuh)
!$OMP THREADPRIVATE(rebcui,rebcuj,reffia,reffib,rriw,rtiw,ryfwca)
!$OMP THREADPRIVATE(ryfwcb,ryfwcc,ryfwcd,ryfwce,ryfwcf)

END MODULE YOECLOP
