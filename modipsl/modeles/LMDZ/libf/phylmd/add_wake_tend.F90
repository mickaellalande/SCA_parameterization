SUBROUTINE add_wake_tend(zddeltat, zddeltaq, zds, zddens, zoccur, text, abortphy)
!===================================================================
! Ajoute les tendances liées aux diverses parametrisations physiques aux
! variables d'etat des poches froides.
!===================================================================
!======================================================================
! Declarations
!======================================================================

USE dimphy, ONLY: klon, klev
USE phys_state_var_mod, ONLY: wake_deltat, wake_deltaq, wake_s, wake_dens

USE print_control_mod, ONLY: prt_level
IMPLICIT none

! Arguments :
!------------
  REAL, DIMENSION(klon, klev),   INTENT (IN)         :: zddeltat, zddeltaq
  REAL, DIMENSION(klon),         INTENT (IN)         :: zds, zddens
  INTEGER, DIMENSION(klon),      INTENT (IN)         :: zoccur
  CHARACTER*(*),                 INTENT (IN)         :: text
  INTEGER,                       INTENT (IN)         :: abortphy

! Local :
!--------

INTEGER                                              :: i, l



     IF (prt_level >= 5) then
        write (*,*) "In add_wake_tend, after ",text
        call flush
     end if

     IF (abortphy==1) RETURN ! on n ajoute pas les tendance si le modele
                              ! a deja plante.

!======================================================================
!    Add tendencies to wake state variables
!======================================================================
         DO l = 1, klev
           DO i = 1, klon
             IF (zoccur(i) .GE. 1) THEN
               wake_deltat(i, l) = wake_deltat(i, l) + zddeltat(i,l)
               wake_deltaq(i, l) = wake_deltaq(i, l) + zddeltaq(i,l)
             ELSE
               wake_deltat(i, l) = 0.
               wake_deltaq(i, l) = 0.
             ENDIF   ! (zoccur(i) .GE. 1)
           END DO
         END DO
         DO i = 1, klon
           IF (zoccur(i) .GE. 1) THEN
             wake_s(i)    = wake_s(i)    + zds(i)
             wake_dens(i) = wake_dens(i) + zddens(i)
           ELSE
             wake_s(i)    = 0.
             wake_dens(i) = 0.
           ENDIF   ! (zoccur(i) .GE. 1)
         END DO

RETURN
END
