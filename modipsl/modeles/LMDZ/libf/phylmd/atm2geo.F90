!
! $Id: atm2geo.F90 2429 2016-01-27 12:43:09Z fairhead $
!
SUBROUTINE atm2geo ( im, jm, pte, ptn, plon, plat, pxx, pyy, pzz )
  USE dimphy
  USE mod_phys_lmdz_para
  IMPLICIT NONE
  INCLUDE 'YOMCST.h'
!
! Change wind local atmospheric coordinates to geocentric
!
  INTEGER, INTENT (in)                 :: im, jm
  REAL, DIMENSION (im,jm), INTENT (in) :: pte, ptn
  REAL, DIMENSION (im,jm), INTENT (in) :: plon, plat
  REAL, DIMENSION (im,jm), INTENT(out) :: pxx, pyy, pzz
  
  REAL :: rad


  rad = rpi / 180.0E0
  
  pxx(:,:) = & 
       - pte(:,:) * SIN(rad * plon(:,:)) &
       - ptn(:,:) * SIN(rad * plat(:,:)) * COS(rad * plon(:,:))

  pyy(:,:) = &
       + pte(:,:) * COS(rad * plon(:,:)) &
       - ptn(:,:) * SIN(rad * plat(:,:)) * SIN(rad * plon(:,:))
  
  pzz(:,:) = &
       + ptn(:,:) * COS(rad * plat (:,:))
  
! Value at North Pole  
  IF (is_north_pole_dyn) THEN
     pxx(:, 1) = - pte (1, 1)
     pyy(:, 1) = - ptn (1, 1) 
     pzz(:, 1) = pzz(1,1)
  ENDIF

! Value at South Pole
  IF (is_south_pole_dyn) THEN
     pxx(:,jm) = pxx(1,jm)
     pyy(:,jm) = pyy(1,jm)
     pzz(:,jm) = pzz(1,jm)
  ENDIF
  
END SUBROUTINE atm2geo
