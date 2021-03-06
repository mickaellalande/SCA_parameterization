!
! $Id $
!
SUBROUTINE cltrac(dtime,coef,t,tr,flux,paprs,pplay,delp, &
                  d_tr,d_tr_dry,flux_tr_dry)                    !jyg

  USE dimphy
  IMPLICIT NONE
!======================================================================
! Auteur(s): O. Boucher (LOA/LMD) date: 19961127
!            inspire de clvent
! Objet: diffusion verticale de traceurs avec flux fixe a la surface
!        ou/et flux du type c-drag
!
! Arguments:
!-----------
! dtime.......input-R- intervalle du temps (en secondes)
! coef........input-R- le coefficient d'echange (m**2/s) l>1
! t...........input-R- temperature (K)
! tr..........input-R- la q. de traceurs
! flux........input-R- le flux de traceurs a la surface
! paprs.......input-R- pression a inter-couche (Pa)
! pplay.......input-R- pression au milieu de couche (Pa)
! delp........input-R- epaisseur de couche (Pa)
! cdrag.......input-R- cdrag pour le flux de surface (non active)
! tr0.........input-R- traceurs a la surface ou dans l'ocean (non active)
! d_tr........output-R- le changement de tr
! d_tr_dry....output-R- le changement de tr du au depot sec (1st layer)
! flux_tr_dry.output-R- depot sec
!!! flux_tr..output-R- flux de tr
!======================================================================
  include "YOMCST.h"
!
! Entree
! 
  REAL,INTENT(IN)                        :: dtime
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: coef
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: t, tr
  REAL,DIMENSION(klon),INTENT(IN)        :: flux !(at/s/m2)
  REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs 
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay, delp
!
! Sorties
!
  REAL ,DIMENSION(klon,klev),INTENT(OUT) :: d_tr
  REAL ,DIMENSION(klon),INTENT(OUT)       :: d_tr_dry          !jyg
  REAL ,DIMENSION(klon),INTENT(OUT)       :: flux_tr_dry       !jyg
!  REAL ,DIMENSION(klon,klev),INTENT(OUT) :: flux_tr
!
! Local
! 
  INTEGER                   :: i, k
  REAL,DIMENSION(klon)      :: cdrag, tr0
  REAL,DIMENSION(klon,klev) :: zx_ctr
  REAL,DIMENSION(klon,klev) :: zx_dtr
  REAL,DIMENSION(klon)      :: zx_buf
  REAL,DIMENSION(klon,klev) :: zx_coef
  REAL,DIMENSION(klon,klev) :: local_tr
  REAL,DIMENSION(klon)      :: zx_alf1,zx_alf2,zx_flux

!======================================================================

  DO k = 1, klev
     DO i = 1, klon
        local_tr(i,k) = tr(i,k)
     ENDDO
  ENDDO

!======================================================================

  DO i = 1, klon
     zx_alf1(i) = (paprs(i,1)-pplay(i,2))/(pplay(i,1)-pplay(i,2))
     zx_alf2(i) = 1.0 - zx_alf1(i)
     flux_tr_dry(i) = -flux(i)*dtime                              !jyg
     zx_flux(i) =  flux_tr_dry(i)*RG                              !jyg
!!     zx_flux(i) =  -flux(i)*dtime*RG                            !jyg
! Pour le moment le flux est prescrit cdrag et zx_coef(1) vaut 0
     cdrag(i) = 0.0 
     tr0(i) = 0.0
     zx_coef(i,1) = cdrag(i)*dtime*RG 
     zx_ctr(i,1)=0.
     zx_dtr(i,1)=0.
  ENDDO

!======================================================================

  DO k = 2, klev
     DO i = 1, klon
        zx_coef(i,k) = coef(i,k)*RG/(pplay(i,k-1)-pplay(i,k))   &
             *(paprs(i,k)*2/(t(i,k)+t(i,k-1))/RD)**2
        zx_coef(i,k) = zx_coef(i,k)*dtime*RG  
     ENDDO
  ENDDO

!======================================================================

  DO i = 1, klon
     zx_buf(i) = delp(i,1) + zx_coef(i,1)*zx_alf1(i) + zx_coef(i,2)
     !
     zx_ctr(i,2) = (local_tr(i,1)*delp(i,1)+                  &
          zx_coef(i,1)*tr0(i)-zx_flux(i))/zx_buf(i)
     !
     zx_dtr(i,2) = (zx_coef(i,2)-zx_alf2(i)*zx_coef(i,1)) /   & 
          zx_buf(i)
     d_tr_dry(i) = -zx_flux(i)/zx_buf(i)                          !jyg
  ENDDO

  DO k = 3, klev
     DO i = 1, klon
        zx_buf(i) = delp(i,k-1) + zx_coef(i,k)      &
             + zx_coef(i,k-1)*(1.-zx_dtr(i,k-1))
        zx_ctr(i,k) = (local_tr(i,k-1)*delp(i,k-1)  & 
             +zx_coef(i,k-1)*zx_ctr(i,k-1) )/zx_buf(i)
        zx_dtr(i,k) = zx_coef(i,k)/zx_buf(i)
     ENDDO
  ENDDO

  DO i = 1, klon
     local_tr(i,klev) = ( local_tr(i,klev)*delp(i,klev) &
          +zx_coef(i,klev)*zx_ctr(i,klev) )             &
          / ( delp(i,klev) + zx_coef(i,klev)            &
          -zx_coef(i,klev)*zx_dtr(i,klev) )
  ENDDO

  DO k = klev-1, 1, -1
     DO i = 1, klon
        local_tr(i,k) = zx_ctr(i,k+1) + zx_dtr(i,k+1)*local_tr(i,k+1)
     ENDDO
  ENDDO

!======================================================================
!== flux_tr est le flux de traceur (positif vers bas)
!      DO i = 1, klon
!         flux_tr(i,1) = zx_coef(i,1)/(RG*dtime)
!      ENDDO
!      DO k = 2, klev
!      DO i = 1, klon
!         flux_tr(i,k) = zx_coef(i,k)/(RG*dtime)
!     .               * (local_tr(i,k)-local_tr(i,k-1))
!      ENDDO
!      ENDDO
!======================================================================
  DO k = 1, klev
     DO i = 1, klon
        d_tr(i,k) = local_tr(i,k) - tr(i,k)
     ENDDO
  ENDDO
  
END SUBROUTINE cltrac
