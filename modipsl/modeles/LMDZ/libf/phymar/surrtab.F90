SUBROUTINE SURRTAB

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** AER'S RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!     -----------------------------------------------------------------

#include "tsmbkind.h"

USE YOERRTAB , ONLY : TRANS, BPADE

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: ITR

!     LOCAL REAL SCALARS
REAL_B :: ZTAU, ZTFN


BPADE=_ONE_/0.278_JPRB
TRANS(0)   =_ONE_
TRANS(5000)=_ZERO_
DO ITR=1,4970
  ZTFN=REAL(ITR)/5000._JPRB
  ZTAU=BPADE*ZTFN/(_ONE_-ZTFN)
  TRANS(ITR)=EXP(-ZTAU)
! write(6,*) ITR,TRANS(ITR)
ENDDO
DO ITR=4971,4999
  TRANS(ITR)=_ZERO_
! write(6,*) ITR,TRANS(ITR)
ENDDO

!     -----------------------------------------------------------------

RETURN
END SUBROUTINE SURRTAB
