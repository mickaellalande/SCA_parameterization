SUBROUTINE borne_var_surf(klon,nbsrf,                        &
         t1,q1,u1,v1,                                        &
         ftsol,pctsrf,                                       &
         t2m, q2m, u10m, v10m,                               &
         zt2m_cor, zq2m_cor, zu10m_cor, zv10m_cor)

IMPLICIT NONE

!==================================================================
! Declarations
!==================================================================

! arguments
INTEGER klon,nbsrf
REAL,DIMENSION(klon),INTENT(IN) :: t1, q1, u1, v1
REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: t2m, q2m, u10m, v10m
REAL,DIMENSION(klon,nbsrf),INTENT(IN) :: ftsol, pctsrf
REAL,DIMENSION (klon) :: zt2m_cor, zq2m_cor, zu10m_cor, zv10m_cor


! local
INTEGER i,nsrf
REAL,DIMENSION (klon,nbsrf) :: t2m_cor, q2m_cor, u10m_cor, v10m_cor

!==================================================================
! Correction of sub surface variables
!==================================================================

DO nsrf=1,nbsrf
   DO i=1,klon
      t2m_cor(i,nsrf)=MIN(t2m(i,nsrf),MAX(t1(i),ftsol(i,nsrf)))
      t2m_cor(i,nsrf)=MAX(t2m_cor(i,nsrf),MIN(t1(i),ftsol(i,nsrf)))
      q2m_cor(i,nsrf)=MAX(q2m(i,nsrf),0.)
      u10m_cor(i,nsrf)=SIGN(MIN(ABS(u1(i)),ABS(u10m(i,nsrf))),u1(i))
      v10m_cor(i,nsrf)=SIGN(MIN(ABS(v1(i)),ABS(v10m(i,nsrf))),v1(i))
   ENDDO
ENDDO

!==================================================================
! Agregation of sub surfaces
!==================================================================

zt2m_cor=0.
zq2m_cor=0.0001
zu10m_cor=0.
zv10m_cor=0.
DO nsrf = 1, nbsrf
   DO i = 1, klon
      zt2m_cor(i)  = zt2m_cor(i)  + t2m_cor(i,nsrf)  * pctsrf(i,nsrf)
      zq2m_cor(i)  = zq2m_cor(i)  + q2m_cor(i,nsrf)  * pctsrf(i,nsrf)
      zu10m_cor(i) = zu10m_cor(i) + u10m_cor(i,nsrf) * pctsrf(i,nsrf)
      zv10m_cor(i) = zv10m_cor(i) + v10m_cor(i,nsrf) * pctsrf(i,nsrf)
   ENDDO
ENDDO

RETURN
END


