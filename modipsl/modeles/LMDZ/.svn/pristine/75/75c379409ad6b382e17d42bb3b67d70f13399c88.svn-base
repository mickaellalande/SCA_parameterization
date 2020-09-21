!
! $Id: $
!
SUBROUTINE phyredem (fichnom)

  USE geometry_mod, ONLY : longitude_deg, latitude_deg
  USE iostart, ONLY: open_restartphy, close_restartphy, put_var, put_field

  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: fichnom

  INTEGER,PARAMETER :: tab_cntrl_len=100
  REAL :: tab_cntrl(tab_cntrl_len)

  ! open file

  CALL open_restartphy(fichnom)

  ! tab_cntrl() contains run parameters

  tab_cntrl(:)=0.0
 

  CALL put_var("controle", "Control parameters", tab_cntrl)

  ! coordinates

  CALL put_field("longitude", "Longitudes on physics grid", longitude_deg)
     
  CALL put_field("latitude", "Latitudes on physics grid", latitude_deg)

  ! close file

  CALL close_restartphy
  !$OMP BARRIER

END SUBROUTINE phyredem
