
! $Id: cv3_cine.F90 2416 2015-12-24 11:58:33Z jyg $

SUBROUTINE cv3_cine(nloc, ncum, nd, icb, inb, pbase, plcl, p, ph, tv, tvp, &
    cina, cinb, plfc)

  ! **************************************************************
  ! *
  ! CV3_CINE                                                    *
  ! *
  ! *
  ! written by   :   Frederique Cheruy                          *
  ! vectorization:   Jean-Yves Grandpeix, 19/06/2003, 11.54.43  *
  ! modified by :                                               *
  ! **************************************************************

  IMPLICIT NONE

  include "YOMCST.h"
  include "cvthermo.h"
  include "cv3param.h"
  ! input:
  INTEGER ncum, nd, nloc
  INTEGER icb(nloc), inb(nloc)
  REAL pbase(nloc), plcl(nloc)
  REAL p(nloc, nd), ph(nloc, nd+1)
  REAL tv(nloc, nd), tvp(nloc, nd)

  ! output
  REAL cina(nloc), cinb(nloc), plfc(nloc)

  ! local variables
  INTEGER il, i, j, k
  INTEGER itop(nloc), ineg(nloc), ilow(nloc)
  INTEGER ifst(nloc), isublcl(nloc)
  LOGICAL lswitch(nloc), lswitch1(nloc), lswitch2(nloc), lswitch3(nloc)
  LOGICAL exist_lfc(nloc)
  REAL dpmax
  REAL deltap, dcin
  REAL buoylcl(nloc), tvplcl(nloc), tvlcl(nloc)
  REAL p0(nloc)
  REAL buoyz(nloc), buoy(nloc, nd)

  ! -------------------------------------------------------------
  ! Initialization
  ! -------------------------------------------------------------
  DO il = 1, ncum
    cina(il) = 0.
    cinb(il) = 0.
  END DO

  ! --------------------------------------------------------------
  ! Recompute buoyancies
  ! --------------------------------------------------------------
  DO k = 1, nd
    DO il = 1, ncum
      ! print*,'tvp tv=',tvp(il,k),tv(il,k)
      buoy(il, k) = tvp(il, k) - tv(il, k)
    END DO
  END DO
  ! ---------------------------------------------------------------

  ! calcul de la flottabilite a LCL (Buoylcl)
  ! ifst = first P-level above lcl
  ! isublcl = highest P-level below lcl.
  ! ---------------------------------------------------------------

  DO il = 1, ncum
    tvplcl(il) = tvp(il, 1)*(plcl(il)/p(il,1))**(2./7.) !For dry air, R/Cp=2/7
  END DO

  DO il = 1, ncum
    IF (plcl(il)>p(il,icb(il))) THEN
      ifst(il) = icb(il)
      isublcl(il) = icb(il) - 1
    ELSE
      ifst(il) = icb(il) + 1
      isublcl(il) = icb(il)
    END IF
  END DO

  DO il = 1, ncum
    tvlcl(il) = tv(il, ifst(il)-1) + (tv(il,ifst(il))-tv(il,ifst(il)-1))*( &
      plcl(il)-p(il,ifst(il)-1))/(p(il,ifst(il))-p(il,ifst(il)-1))
  END DO

  DO il = 1, ncum
    buoylcl(il) = tvplcl(il) - tvlcl(il)
  END DO

  ! ---------------------------------------------------------------
  ! premiere couche contenant un  niveau de flotabilite positive
  ! et premiere couche contenant un  niveau de flotabilite negative
  ! au dessus du niveau de condensation
  ! ---------------------------------------------------------------
  DO il = 1, ncum
    itop(il) = nl - 1
    ineg(il) = nl - 1
    exist_lfc(il) = .FALSE.
  END DO
  DO k = nl - 1, 1, -1
    DO il = 1, ncum
      IF (k>=ifst(il)) THEN
        IF (buoy(il,k)>0.) THEN
          itop(il) = k
          exist_lfc(il) = .TRUE.
        ELSE
          ineg(il) = k
        END IF
      END IF
    END DO
  END DO

  ! ---------------------------------------------------------------
  ! When there is no positive buoyancy level, set Plfc, Cina and Cinb
  ! to arbitrary extreme values.
  ! ---------------------------------------------------------------
  DO il = 1, ncum
    IF (.NOT. exist_lfc(il)) THEN
      plfc(il) = 1.111
      cinb(il) = -1111.
      cina(il) = -1112.
    END IF
  END DO


  ! ---------------------------------------------------------------
  ! -- Two cases : BUOYlcl >= 0 and BUOYlcl < 0.
  ! ---------------------------------------------------------------

  ! --------------------
  ! -- 1.0 BUOYlcl >=0.
  ! --------------------

  dpmax = 50.
  DO il = 1, ncum
    lswitch1(il) = buoylcl(il) >= 0. .AND. exist_lfc(il)
    lswitch(il) = lswitch1(il)
  END DO

  ! 1.1 No inhibition case
  ! ----------------------
  ! If buoyancy is positive at LCL and stays positive over a large enough
  ! pressure interval (=DPMAX), inhibition is set to zero,

  DO il = 1, ncum
    IF (lswitch(il)) THEN
      IF (p(il,ineg(il))<p(il,icb(il))-dpmax) THEN
        plfc(il) = plcl(il)
        cina(il) = 0.
        cinb(il) = 0.
      END IF
    END IF
  END DO

  ! 1.2 Upper inhibition only case
  ! ------------------------------
  DO il = 1, ncum
    lswitch2(il) = p(il, ineg(il)) >= p(il, icb(il)) - dpmax
    lswitch(il) = lswitch1(il) .AND. lswitch2(il)
  END DO

  ! 1.2.1 Recompute itop (=1st layer with positive buoyancy above ineg)
  ! -------------------------------------------------------------------

  DO il = 1, ncum
    IF (lswitch(il)) THEN
      itop(il) = nl - 1
    END IF
  END DO

  DO k = nl, 1, -1
    DO il = 1, ncum
      IF (lswitch(il)) THEN
        IF (k>=ineg(il) .AND. buoy(il,k)>0) THEN
          itop(il) = k
        END IF
      END IF
    END DO
  END DO

  ! If there is no layer with positive buoyancy above ineg, set Plfc, 
  ! Cina and Cinb to arbitrary extreme values.
  DO il = 1, ncum
    IF (lswitch(il) .AND. itop(il) == nl - 1) THEN
      plfc(il) = 1.121
      cinb(il) = -1121.
      cina(il) = -1122.
    END IF
  END DO

  DO il = 1, ncum
    lswitch3(il) = itop(il) < nl -1
    lswitch(il) = lswitch1(il) .AND. lswitch2(il) .AND. lswitch3(il)
  END DO

  DO il = 1, ncum
    IF (lswitch(il)) THEN
      cinb(il) = 0.

      ! 1.2.2  Calcul de la pression du niveau de flot. nulle juste au-dessus
      ! de LCL
      ! ---------------------------------------------------------------------------
      IF (ineg(il)>isublcl(il)+1) THEN
        ! In order to get P0, one may interpolate linearly buoyancies
        ! between P(ineg) and P(ineg-1).
        p0(il) = (buoy(il,ineg(il))*p(il,ineg(il)-1)-buoy(il,ineg(il)-1)*p(il,ineg(il)))/ &
          (buoy(il,ineg(il))-buoy(il,ineg(il)-1))
      ELSE
        ! In order to get P0, one has to interpolate between P(ineg) and
        ! Plcl.
        p0(il) = (buoy(il,ineg(il))*plcl(il)-buoylcl(il)*p(il,ineg(il)))/ &
          (buoy(il,ineg(il))-buoylcl(il))
      END IF
    END IF
  END DO

  ! 1.2.3 Computation of PLFC
  ! -------------------------
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      plfc(il) = (buoy(il,itop(il))*p(il,itop(il)-1)-buoy(il,itop( &
        il)-1)*p(il,itop(il)))/(buoy(il,itop(il))-buoy(il,itop(il)-1))
    END IF
  END DO

  ! 1.2.4 Computation of CINA
  ! -------------------------

  ! Upper part of CINA : integral from P(itop-1) to Plfc
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      deltap = p(il, itop(il)-1) - plfc(il)
      dcin = rd*buoy(il, itop(il)-1)*deltap/(p(il,itop(il)-1)+plfc(il))
      cina(il) = min(0., dcin)
    END IF
  END DO

  ! Middle part of CINA : integral from P(ineg) to P(itop-1)
  DO k = 1, nl
    DO il = 1, ncum
      IF (lswitch(il)) THEN
        IF (k>=ineg(il) .AND. k<=itop(il)-2) THEN
          deltap = p(il, k) - p(il, k+1)
          dcin = 0.5*rd*(buoy(il,k)+buoy(il,k+1))*deltap/ph(il, k+1)
          cina(il) = cina(il) + min(0., dcin)
        END IF
      END IF
    END DO
  END DO

  ! Lower part of CINA : integral from P0 to P(ineg)
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      deltap = p0(il) - p(il, ineg(il))
      dcin = rd*buoy(il, ineg(il))*deltap/(p(il,ineg(il))+p0(il))
      cina(il) = cina(il) + min(0., dcin)
    END IF
  END DO


  ! ------------------
  ! -- 2.0 BUOYlcl <0.
  ! ------------------

  DO il = 1, ncum
    lswitch1(il) = buoylcl(il) < 0. .AND. exist_lfc(il)
    lswitch(il) = lswitch1(il)
  END DO

  ! 2.0.1 Premiere  couche ou la flotabilite est negative au dessus du sol
  ! ----------------------------------------------------
  ! au cas ou elle existe  sinon ilow=1 (nk apres)
  ! on suppose que la parcelle part de la premiere couche

  DO il = 1, ncum
    IF (lswitch(il)) THEN
      ilow(il) = 1
    END IF
  END DO

  DO k = nl, 1, -1
    DO il = 1, ncum
      IF (lswitch(il) .AND. k<=icb(il)-1) THEN
        IF (buoy(il,k)<0.) THEN
          ilow(il) = k
        END IF
      END IF
    END DO
  END DO

  ! 2.0.2  Calcul de la pression du niveau de flot. nulle sous le nuage
  ! ----------------------------------------------------
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      IF (ilow(il)>1) THEN
        p0(il) = (buoy(il,ilow(il))*p(il,ilow(il)-1)-buoy(il,ilow( &
          il)-1)*p(il,ilow(il)))/(buoy(il,ilow(il))-buoy(il,ilow(il)-1))
        buoyz(il) = 0.
      ELSE
        p0(il) = p(il, 1)
        buoyz(il) = buoy(il, 1)
      END IF
    END IF
  END DO

  ! 2.1. Computation of CINB
  ! -----------------------

  DO il = 1, ncum
    lswitch2(il) = (isublcl(il)==1 .AND. ilow(il)==1) .OR. &
      (isublcl(il)==ilow(il)-1)
    lswitch(il) = lswitch1(il) .AND. lswitch2(il)
  END DO
  ! c      IF (    (isublcl .EQ. 1 .AND. ilow .EQ. 1)
  ! c     $    .OR.(isublcl .EQ. ilow-1)) THEN

  ! 2.1.1 First case : Plcl just above P0
  ! -------------------------------------
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      deltap = p0(il) - plcl(il)
      dcin = rd*(buoyz(il)+buoylcl(il))*deltap/(p0(il)+plcl(il))
      cinb(il) = min(0., dcin)
    END IF
  END DO

  DO il = 1, ncum
    lswitch(il) = lswitch1(il) .AND. .NOT. lswitch2(il)
  END DO
  ! c      ELSE

  ! 2.1.2 Second case : there is at least one P-level between P0 and Plcl
  ! ---------------------------------------------------------------------

  ! Lower part of CINB : integral from P0 to P(ilow)
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      deltap = p0(il) - p(il, ilow(il))
      dcin = rd*(buoyz(il)+buoy(il,ilow(il)))*deltap/(p0(il)+p(il,ilow(il)))
      cinb(il) = min(0., dcin)
    END IF
  END DO


  ! Middle part of CINB : integral from P(ilow) to P(isublcl)
  ! c      DO k = ilow,isublcl-1
  DO k = 1, nl
    DO il = 1, ncum
      IF (lswitch(il) .AND. k>=ilow(il) .AND. k<=isublcl(il)-1) THEN
        deltap = p(il, k) - p(il, k+1)
        dcin = 0.5*rd*(buoy(il,k)+buoy(il,k+1))*deltap/ph(il, k+1)
        cinb(il) = cinb(il) + min(0., dcin)
      END IF
    END DO
  END DO

  ! Upper part of CINB : integral from P(isublcl) to Plcl
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      deltap = p(il, isublcl(il)) - plcl(il)
      dcin = rd*(buoy(il,isublcl(il))+buoylcl(il))*deltap/ &
        (p(il,isublcl(il))+plcl(il))
      cinb(il) = cinb(il) + min(0., dcin)
    END IF
  END DO


  ! c      ENDIF

  ! 2.2 Computation of CINA
  ! ---------------------

  DO il = 1, ncum
    lswitch2(il) = plcl(il) > p(il, itop(il)-1)
    lswitch(il) = lswitch1(il) .AND. lswitch2(il)
  END DO

  ! 2.2.1 FIrst case : Plcl > P(itop-1)
  ! ---------------------------------
  ! In order to get Plfc, one may interpolate linearly buoyancies
  ! between P(itop) and P(itop-1).
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      plfc(il) = (buoy(il,itop(il))*p(il,itop(il)-1)-buoy(il,itop( &
        il)-1)*p(il,itop(il)))/(buoy(il,itop(il))-buoy(il,itop(il)-1))
    END IF
  END DO

  ! Upper part of CINA : integral from P(itop-1) to Plfc
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      deltap = p(il, itop(il)-1) - plfc(il)
      dcin = rd*buoy(il, itop(il)-1)*deltap/(p(il,itop(il)-1)+plfc(il))
      cina(il) = min(0., dcin)
    END IF
  END DO

  ! Middle part of CINA : integral from P(icb+1) to P(itop-1)
  DO k = 1, nl
    DO il = 1, ncum
      IF (lswitch(il) .AND. k>=icb(il)+1 .AND. k<=itop(il)-2) THEN
        deltap = p(il, k) - p(il, k+1)
        dcin = 0.5*rd*(buoy(il,k)+buoy(il,k+1))*deltap/ph(il, k+1)
        cina(il) = cina(il) + min(0., dcin)
      END IF
    END DO
  END DO

  ! Lower part of CINA : integral from Plcl to P(icb+1)
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      IF (plcl(il)>p(il,icb(il))) THEN
        IF (icb(il)<itop(il)-1) THEN
          deltap = p(il, icb(il)) - p(il, icb(il)+1)
          dcin = 0.5*rd*(buoy(il,icb(il))+buoy(il,icb(il)+1))*deltap/ &
            ph(il, icb(il)+1)
          cina(il) = cina(il) + min(0., dcin)
        END IF

        deltap = plcl(il) - p(il, icb(il))
        dcin = rd*(buoylcl(il)+buoy(il,icb(il)))*deltap/ &
          (plcl(il)+p(il,icb(il)))
        cina(il) = cina(il) + min(0., dcin)
      ELSE
        deltap = plcl(il) - p(il, icb(il)+1)
        dcin = rd*(buoylcl(il)+buoy(il,icb(il)+1))*deltap/ &
          (plcl(il)+p(il,icb(il)+1))
        cina(il) = cina(il) + min(0., dcin)
      END IF
    END IF
  END DO

  DO il = 1, ncum
    lswitch(il) = lswitch1(il) .AND. .NOT. lswitch2(il)
  END DO
  ! c      ELSE

  ! 2.2.2 Second case : Plcl lies between P(itop-1) and P(itop);
  ! ----------------------------------------------------------
  ! In order to get Plfc, one has to interpolate between P(itop) and Plcl.
  DO il = 1, ncum
    IF (lswitch(il)) THEN
      plfc(il) = (buoy(il,itop(il))*plcl(il)-buoylcl(il)*p(il,itop(il)))/ &
        (buoy(il,itop(il))-buoylcl(il))
    END IF
  END DO

  DO il = 1, ncum
    IF (lswitch(il)) THEN
      deltap = plcl(il) - plfc(il)
      dcin = rd*buoylcl(il)*deltap/(plcl(il)+plfc(il))
      cina(il) = min(0., dcin)
    END IF
  END DO
  ! c      ENDIF



  RETURN
END SUBROUTINE cv3_cine
