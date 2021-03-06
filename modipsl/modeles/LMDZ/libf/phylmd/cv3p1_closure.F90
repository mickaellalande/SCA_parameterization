
! $Id: cv3p1_closure.F90 2826 2017-03-16 13:31:15Z fhourdin $

SUBROUTINE cv3p1_closure(nloc, ncum, nd, icb, inb, pbase, plcl, p, ph, tv, &
    tvp, buoy, supmax, ok_inhib, ale, alp, omega,sig, w0, ptop2, cape, cin, m, &
    iflag, coef, plim1, plim2, asupmax, supmax0, asupmaxmin, cbmf, plfc, &
    wbeff)


  ! **************************************************************
  ! *
  ! CV3P1_CLOSURE                                               *
  ! Ale & Alp Closure of Convect3              *
  ! *
  ! written by   :   Kerry Emanuel                              *
  ! vectorization:   S. Bony                                    *
  ! modified by :    Jean-Yves Grandpeix, 18/06/2003, 19.32.10  *
  ! Julie Frohwirth,     14/10/2005  17.44.22  *
  ! **************************************************************

  USE print_control_mod, ONLY: prt_level, lunout
  IMPLICIT NONE

  include "cvthermo.h"
  include "cv3param.h"
  include "YOMCST2.h"
  include "YOMCST.h"
  include "conema3.h"

  ! input:
  INTEGER, INTENT (IN)                               :: ncum, nd, nloc
  INTEGER, DIMENSION (nloc), INTENT (IN)             :: icb, inb
  REAL, DIMENSION (nloc), INTENT (IN)                :: pbase, plcl
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: p
  REAL, DIMENSION (nloc, nd+1), INTENT (IN)          :: ph
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: tv, tvp, buoy
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: supmax
  LOGICAL, INTENT (IN)                               :: ok_inhib ! enable convection inhibition by dryness
  REAL, DIMENSION (nloc), INTENT (IN)                :: ale, alp
  REAL, DIMENSION (nloc, nd), INTENT (IN)            :: omega

  ! input/output:
  REAL, DIMENSION (nloc, nd), INTENT (INOUT)         :: sig, w0
  REAL, DIMENSION (nloc), INTENT (INOUT)             :: ptop2

  ! output:
  REAL, DIMENSION (nloc), INTENT (OUT)               :: cape, cin
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: m
  REAL, DIMENSION (nloc), INTENT (OUT)               :: plim1, plim2
  REAL, DIMENSION (nloc, nd), INTENT (OUT)           :: asupmax
  REAL, DIMENSION (nloc), INTENT (OUT)               :: supmax0
  REAL, DIMENSION (nloc), INTENT (OUT)               :: asupmaxmin
  REAL, DIMENSION (nloc), INTENT (OUT)               :: cbmf, plfc
  REAL, DIMENSION (nloc), INTENT (OUT)               :: wbeff
  INTEGER, DIMENSION (nloc), INTENT (OUT)            :: iflag

  ! local variables:
  INTEGER il, i, j, k, icbmax, i0(nloc), klfc(nloc)
  REAL deltap, fac, w, amu
  REAL rhodp, dz
  REAL pbmxup
  REAL dtmin(nloc, nd), sigold(nloc, nd)
  REAL coefmix(nloc, nd)
  REAL pzero(nloc), ptop2old(nloc)
  REAL cina(nloc), cinb(nloc)
  INTEGER ibeg(nloc)
  INTEGER nsupmax(nloc)
  REAL supcrit, temp(nloc, nd)
  REAL p1(nloc), pmin(nloc)
  REAL asupmax0(nloc)
  LOGICAL ok(nloc)
  REAL siglim(nloc, nd), wlim(nloc, nd), mlim(nloc, nd)
  REAL wb2(nloc)
  REAL cbmflim(nloc), cbmf1(nloc), cbmfmax(nloc)
  REAL cbmflast(nloc)
  REAL coef(nloc)
  REAL xp(nloc), xq(nloc), xr(nloc), discr(nloc), b3(nloc), b4(nloc)
  REAL theta(nloc), bb(nloc)
  REAL term1, term2, term3
  REAL alp2(nloc) ! Alp with offset
  !CR: variables for new erosion of adiabiatic ascent
  REAL mad(nloc, nd), me(nloc, nd), betalim(nloc, nd), beta_coef(nloc, nd)
  REAL med(nloc, nd), md(nloc,nd)
!jyg<
! coef_peel is now in the common cv3_param 
!!  REAL coef_peel
!!  PARAMETER (coef_peel=0.25)
!>jyg

  REAL sigmax
  PARAMETER (sigmax=0.1)

  CHARACTER (LEN=20) :: modname = 'cv3p1_closure'
  CHARACTER (LEN=80) :: abort_message

  ! print *,' -> cv3p1_closure, Ale ',ale(1)


  ! -------------------------------------------------------
  ! -- Initialization
  ! -------------------------------------------------------


  DO il = 1, ncum
    alp2(il) = max(alp(il), 1.E-5)
    ! IM
    alp2(il) = max(alp(il), 1.E-12)
  END DO

  pbmxup = 50. ! PBMXUP+PBCRIT = cloud depth above which mixed updraughts
  ! exist (if any)

  IF (prt_level>=20) PRINT *, 'cv3p1_param nloc ncum nd icb inb nl', nloc, &
    ncum, nd, icb(nloc), inb(nloc), nl
  DO k = 1, nd          !jyg: initialization up to nd
    DO il = 1, ncum
      m(il, k) = 0.0
    END DO
  END DO

!CR: initializations for erosion of adiabatic ascent 
  DO k = 1,nd           !jyg: initialization up to nd
    DO il = 1, ncum
        mad(il,k)=0.
        me(il,k)=0.
        betalim(il,k)=1.
        wlim(il,k)=0.
    ENDDO
  ENDDO

  ! -------------------------------------------------------
  ! -- Reset sig(i) and w0(i) for i>inb and i<icb
  ! -------------------------------------------------------

  ! update sig and w0 above LNB:

  DO k = 1, nl - 1
    DO il = 1, ncum
      IF ((inb(il)<(nl-1)) .AND. (k>=(inb(il)+1))) THEN
        sig(il, k) = beta*sig(il, k) + 2.*alpha*buoy(il, inb(il))*abs(buoy(il &
          ,inb(il)))
        sig(il, k) = amax1(sig(il,k), 0.0)
        w0(il, k) = beta*w0(il, k)
      END IF
    END DO
  END DO

  ! if(prt.level.GE.20) print*,'cv3p1_param apres 100'
  ! compute icbmax:

  icbmax = 2
  DO il = 1, ncum
    icbmax = max(icbmax, icb(il))
  END DO
  ! if(prt.level.GE.20) print*,'cv3p1_param apres 200'

  ! update sig and w0 below cloud base:

  DO k = 1, icbmax
    DO il = 1, ncum
      IF (k<=icb(il)) THEN
        sig(il, k) = beta*sig(il, k) - 2.*alpha*buoy(il, icb(il))*buoy(il, &
          icb(il))
        sig(il, k) = amax1(sig(il,k), 0.0)
        w0(il, k) = beta*w0(il, k)
      END IF
    END DO
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 300'
  ! -------------------------------------------------------------
  ! -- Reset fractional areas of updrafts and w0 at initial time
  ! -- and after 10 time steps of no convection
  ! -------------------------------------------------------------

  DO k = 1, nl - 1
    DO il = 1, ncum
      IF (sig(il,nd)<1.5 .OR. sig(il,nd)>12.0) THEN
        sig(il, k) = 0.0
        w0(il, k) = 0.0
      END IF
    END DO
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 400'

  ! -------------------------------------------------------------
  ! jyg1
  ! --  Calculate adiabatic ascent top pressure (ptop)
  ! -------------------------------------------------------------


  ! c 1. Start at first level where precipitations form
  DO il = 1, ncum
    pzero(il) = plcl(il) - pbcrit
  END DO

  ! c 2. Add offset
  DO il = 1, ncum
    pzero(il) = pzero(il) - pbmxup
  END DO
  DO il = 1, ncum
    ptop2old(il) = ptop2(il)
  END DO

  DO il = 1, ncum
    ! CR:c est quoi ce 300??
    p1(il) = pzero(il) - 300.
  END DO

  ! compute asupmax=abs(supmax) up to lnm+1

  DO il = 1, ncum
    ok(il) = .TRUE.
    nsupmax(il) = inb(il)
  END DO

  DO i = 1, nl
    DO il = 1, ncum
      IF (i>icb(il) .AND. i<=inb(il)) THEN
        IF (p(il,i)<=pzero(il) .AND. supmax(il,i)<0 .AND. ok(il)) THEN
          nsupmax(il) = i
          ok(il) = .FALSE.
        END IF ! end IF (P(i) ...  )
      END IF ! end IF (icb+1 le i le inb)
    END DO
  END DO

  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 2.'
  DO i = 1, nl
    DO il = 1, ncum
      asupmax(il, i) = abs(supmax(il,i))
    END DO
  END DO


  DO il = 1, ncum
    asupmaxmin(il) = 10.
    pmin(il) = 100.
    ! IM ??
    asupmax0(il) = 0.
  END DO

  ! c 3.  Compute in which level is Pzero

  ! IM bug      i0 = 18
  DO il = 1, ncum
    i0(il) = nl
  END DO

  DO i = 1, nl
    DO il = 1, ncum
      IF (i>icb(il) .AND. i<=inb(il)) THEN
        IF (p(il,i)<=pzero(il) .AND. p(il,i)>=p1(il)) THEN
          IF (pzero(il)>p(il,i) .AND. pzero(il)<p(il,i-1)) THEN
            i0(il) = i
          END IF
        END IF
      END IF
    END DO
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 3.'

  ! c 4.  Compute asupmax at Pzero

  DO i = 1, nl
    DO il = 1, ncum
      IF (i>icb(il) .AND. i<=inb(il)) THEN
        IF (p(il,i)<=pzero(il) .AND. p(il,i)>=p1(il)) THEN
          asupmax0(il) = ((pzero(il)-p(il,i0(il)-1))*asupmax(il,i0(il))-( &
            pzero(il)-p(il,i0(il)))*asupmax(il,i0(il)-1))/(p(il,i0(il))-p(il, &
            i0(il)-1))
        END IF
      END IF
    END DO
  END DO


  DO i = 1, nl
    DO il = 1, ncum
      IF (p(il,i)==pzero(il)) THEN
        asupmax(i, il) = asupmax0(il)
      END IF
    END DO
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 4.'

  ! c 5. Compute asupmaxmin, minimum of asupmax

  DO i = 1, nl
    DO il = 1, ncum
      IF (i>icb(il) .AND. i<=inb(il)) THEN
        IF (p(il,i)<=pzero(il) .AND. p(il,i)>=p1(il)) THEN
          IF (asupmax(il,i)<asupmaxmin(il)) THEN
            asupmaxmin(il) = asupmax(il, i)
            pmin(il) = p(il, i)
          END IF
        END IF
      END IF
    END DO
  END DO

  DO il = 1, ncum
    ! IM
    IF (prt_level>=20) THEN
      PRINT *, 'cv3p1_closure il asupmax0 asupmaxmin', il, asupmax0(il), &
        asupmaxmin(il), pzero(il), pmin(il)
    END IF
    IF (asupmax0(il)<asupmaxmin(il)) THEN
      asupmaxmin(il) = asupmax0(il)
      pmin(il) = pzero(il)
    END IF
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 5.'


  ! Compute Supmax at Pzero

  DO i = 1, nl
    DO il = 1, ncum
      IF (i>icb(il) .AND. i<=inb(il)) THEN
        IF (p(il,i)<=pzero(il)) THEN
          supmax0(il) = ((p(il,i)-pzero(il))*asupmax(il,i-1)-(p(il, &
            i-1)-pzero(il))*asupmax(il,i))/(p(il,i)-p(il,i-1))
          GO TO 425
        END IF ! end IF (P(i) ... )
      END IF ! end IF (icb+1 le i le inb)
    END DO
  END DO

425 CONTINUE
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 425.'

  ! c 6. Calculate ptop2

  DO il = 1, ncum
    IF (asupmaxmin(il)<supcrit1) THEN
      ptop2(il) = pmin(il)
    END IF

    IF (asupmaxmin(il)>supcrit1 .AND. asupmaxmin(il)<supcrit2) THEN
      ptop2(il) = ptop2old(il)
    END IF

    IF (asupmaxmin(il)>supcrit2) THEN
      ptop2(il) = ph(il, inb(il))
    END IF
  END DO

  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 6.'

  ! c 7. Compute multiplying factor for adiabatic updraught mass flux


  IF (ok_inhib) THEN

    DO i = 1, nl
      DO il = 1, ncum
        IF (i<=nl) THEN
          coefmix(il, i) = (min(ptop2(il),ph(il,i))-ph(il,i))/(ph(il,i+1)-ph( &
            il,i))
          coefmix(il, i) = min(coefmix(il,i), 1.)
        END IF
      END DO
    END DO


  ELSE ! when inhibition is not taken into account, coefmix=1



    DO i = 1, nl
      DO il = 1, ncum
        IF (i<=nl) THEN
          coefmix(il, i) = 1.
        END IF
      END DO
    END DO

  END IF ! ok_inhib
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 7.'
  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------


  ! jyg2

  ! ==========================================================================


  ! -------------------------------------------------------------
  ! -- Calculate convective inhibition (CIN)
  ! -------------------------------------------------------------

  ! do i=1,nloc
  ! print*,'avant cine p',pbase(i),plcl(i)
  ! enddo
  ! do j=1,nd
  ! do i=1,nloc
  ! print*,'avant cine t',tv(i),tvp(i)
  ! enddo
  ! enddo
  CALL cv3_cine(nloc, ncum, nd, icb, inb, pbase, plcl, p, ph, tv, tvp, cina, &
    cinb, plfc)

  DO il = 1, ncum
    cin(il) = cina(il) + cinb(il)
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres cv3_cine'
  ! -------------------------------------------------------------
  ! --Update buoyancies to account for Ale
  ! -------------------------------------------------------------

  CALL cv3_buoy(nloc, ncum, nd, icb, inb, pbase, plcl, p, ph, ale, cin, tv, &
    tvp, buoy)
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres cv3_buoy'

  ! -------------------------------------------------------------
  ! -- Calculate convective available potential energy (cape),
  ! -- vertical velocity (w), fractional area covered by
  ! -- undilute updraft (sig), and updraft mass flux (m)
  ! -------------------------------------------------------------

  DO il = 1, ncum
    cape(il) = 0.0
  END DO

  ! compute dtmin (minimum buoyancy between ICB and given level k):

  DO k = 1, nl
    DO il = 1, ncum
      dtmin(il, k) = 100.0
    END DO
  END DO

  DO k = 1, nl
    DO j = minorig, nl
      DO il = 1, ncum
        IF ((k>=(icb(il)+1)) .AND. (k<=inb(il)) .AND. (j>=icb(il)) .AND. (j<= &
            (k-1))) THEN
          dtmin(il, k) = amin1(dtmin(il,k), buoy(il,j))
        END IF
      END DO
    END DO
  END DO

  ! the interval on which cape is computed starts at pbase :

  DO k = 1, nl
    DO il = 1, ncum

      IF ((k>=(icb(il)+1)) .AND. (k<=inb(il))) THEN
        IF (iflag_mix_adiab.eq.1) THEN
!CR:computation of cape from LCL: keep flag or to modify in all cases?
        deltap = min(plcl(il), ph(il,k-1)) - min(plcl(il), ph(il,k))
        ELSE
        deltap = min(pbase(il), ph(il,k-1)) - min(pbase(il), ph(il,k))
        ENDIF
        cape(il) = cape(il) + rrd*buoy(il, k-1)*deltap/p(il, k-1)
        cape(il) = amax1(0.0, cape(il))
        sigold(il, k) = sig(il, k)


        ! jyg       Coefficient coefmix limits convection to levels where a
        ! sufficient
        ! fraction of mixed draughts are ascending.
        siglim(il, k) = coefmix(il, k)*alpha1*dtmin(il, k)*abs(dtmin(il,k))
        siglim(il, k) = amax1(siglim(il,k), 0.0)
        siglim(il, k) = amin1(siglim(il,k), 0.01)
        ! c         fac=AMIN1(((dtcrit-dtmin(il,k))/dtcrit),1.0)
        fac = 1.
        wlim(il, k) = fac*sqrt(cape(il))
        amu = siglim(il, k)*wlim(il, k)
        rhodp = 0.007*p(il, k)*(ph(il,k)-ph(il,k+1))/tv(il, k)
        mlim(il, k) = amu*rhodp
        ! print*, 'siglim ', k,siglim(1,k)
      END IF

    END DO
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 600'

  DO il = 1, ncum
    ! IM beg
    IF (prt_level>=20) THEN
      PRINT *, 'cv3p1_closure il icb mlim ph ph+1 ph+2', il, icb(il), &
        mlim(il, icb(il)+1), ph(il, icb(il)), ph(il, icb(il)+1), &
        ph(il, icb(il)+2)
    END IF

    IF (icb(il)+1<=inb(il)) THEN
      ! IM end
      mlim(il, icb(il)) = 0.5*mlim(il, icb(il)+1)*(ph(il,icb(il))-ph(il,icb( &
        il)+1))/(ph(il,icb(il)+1)-ph(il,icb(il)+2))
      ! IM beg
    END IF !(icb(il.le.inb(il))) then
    ! IM end
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres 700'

  ! jyg1
  ! ------------------------------------------------------------------------
  ! c     Correct mass fluxes so that power used to overcome CIN does not
  ! c     exceed Power Available for Lifting (PAL).
  ! ------------------------------------------------------------------------

  DO il = 1, ncum
    cbmflim(il) = 0.
    cbmf(il) = 0.
  END DO

  ! c 1. Compute cloud base mass flux of elementary system (Cbmf0=Cbmflim)

  DO k = 1, nl
    DO il = 1, ncum
      ! old       IF (k .ge. icb(il) .and. k .le. inb(il)) THEN
      ! IM        IF (k .ge. icb(il)+1 .and. k .le. inb(il)) THEN
      IF (k>=icb(il) .AND. k<=inb(il) & !cor jyg
          .AND. icb(il)+1<=inb(il)) THEN !cor jyg
        cbmflim(il) = cbmflim(il) + mlim(il, k)
      END IF
    END DO
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres cbmflim'

  ! 1.5 Compute cloud base mass flux given by Alp closure (Cbmf1), maximum
  !     allowed mass flux (Cbmfmax) and final target mass flux (Cbmf)
  !     Cbmf is set to zero if Cbmflim (the mass flux of elementary cloud)
  !     is exceedingly small.

  DO il = 1, ncum
    wb2(il) = sqrt(2.*max(ale(il)+cin(il),0.))
  END DO

  DO il = 1, ncum
    IF (plfc(il)<100.) THEN
      ! This is an irealistic value for plfc => no calculation of wbeff
      wbeff(il) = 100.1
    ELSE
      ! Calculate wbeff
      IF (flag_wb==0) THEN
        wbeff(il) = wbmax
      ELSE IF (flag_wb==1) THEN
        wbeff(il) = wbmax/(1.+500./(ph(il,1)-plfc(il)))
      ELSE IF (flag_wb==2) THEN
        wbeff(il) = wbmax*(0.01*(ph(il,1)-plfc(il)))**2
      ELSE ! Option provisoire ou le iflag_wb/10 est considere comme une vitesse
        wbeff(il) = flag_wb*0.01+wbmax/(1.+500./(ph(il,1)-plfc(il)))
      END IF
    END IF
  END DO

!CR:Compute k at plfc
  DO il=1,ncum
           klfc(il)=nl
  ENDDO
  DO k=1,nl
     DO il=1,ncum
        if ((plfc(il).lt.ph(il,k)).and.(plfc(il).ge.ph(il,k+1))) then
           klfc(il)=k
        endif
     ENDDO
  ENDDO
!RC

  DO il = 1, ncum
    ! jyg    Modification du coef de wb*wb pour conformite avec papier Wake
    ! c       cbmf1(il) = alp2(il)/(0.5*wb*wb-Cin(il))
    cbmf1(il) = alp2(il)/(2.*wbeff(il)*wbeff(il)-cin(il))
!CR: Add large-scale component to the mass-flux
!encore connu sous le nom "Experience du tube de dentifrice"
    if ((coef_clos_ls.gt.0.).and.(plfc(il).gt.0.)) then 
       cbmf1(il) = cbmf1(il) - coef_clos_ls*min(0.,1./RG*omega(il,klfc(il)))
    endif
!RC
    IF (cbmf1(il)==0 .AND. alp2(il)/=0.) THEN
      WRITE (lunout, *) 'cv3p1_closure cbmf1=0 and alp NE 0 il alp2 alp cin ' &
        , il, alp2(il), alp(il), cin(il)
      abort_message = ''
      CALL abort_physic(modname, abort_message, 1)
    END IF
    cbmfmax(il) = sigmax*wb2(il)*100.*p(il, icb(il))/(rrd*tv(il,icb(il)))
  END DO

  DO il = 1, ncum
    IF (cbmflim(il)>1.E-6) THEN
      ! ATTENTION TEST CR
      ! if (cbmfmax(il).lt.1.e-12) then
      cbmf(il) = min(cbmf1(il), cbmfmax(il))
      ! else
      ! cbmf(il) = cbmf1(il)
      ! endif
      ! print*,'cbmf',cbmf1(il),cbmfmax(il)
    END IF
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres cbmflim_testCR'

  ! c 2. Compute coefficient and apply correction

  DO il = 1, ncum
    coef(il) = (cbmf(il)+1.E-10)/(cbmflim(il)+1.E-10)
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres coef_plantePLUS'

  DO k = 1, nl
    DO il = 1, ncum
      IF (k>=icb(il)+1 .AND. k<=inb(il)) THEN
        amu = beta*sig(il, k)*w0(il, k) + (1.-beta)*coef(il)*siglim(il, k)* &
          wlim(il, k)
        w0(il, k) = wlim(il, k)
        w0(il, k) = max(w0(il,k), 1.E-10)
        sig(il, k) = amu/w0(il, k)
        sig(il, k) = min(sig(il,k), 1.)
        ! c         amu = 0.5*(SIG(il,k)+sigold(il,k))*W0(il,k)
        m(il, k) = amu*0.007*p(il, k)*(ph(il,k)-ph(il,k+1))/tv(il, k)
      END IF
    END DO
  END DO
  ! jyg2
  DO il = 1, ncum
    w0(il, icb(il)) = 0.5*w0(il, icb(il)+1)
    m(il, icb(il)) = 0.5*m(il, icb(il)+1)*(ph(il,icb(il))-ph(il,icb(il)+1))/ &
      (ph(il,icb(il)+1)-ph(il,icb(il)+2))
    sig(il, icb(il)) = sig(il, icb(il)+1)
    sig(il, icb(il)-1) = sig(il, icb(il))
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres w0_sig_M'

!CR: new erosion of adiabatic ascent: modification of m
!computation of the sum of ascending fluxes 
  IF (iflag_mix_adiab.eq.1) THEN

!Verification sum(me)=sum(m)
  DO k = 1,nd                         !jyg: initialization up to nd
    DO il = 1, ncum
       md(il,k)=0.
       med(il,k)=0.
    ENDDO
  ENDDO

  DO k = nl,1,-1
    DO il = 1, ncum
           md(il,k)=md(il,k+1)+m(il,k+1)
    ENDDO
  ENDDO

  DO k = nl,1,-1
    DO il = 1, ncum
        IF ((k>=(icb(il))) .AND. (k<=inb(il))) THEN
           mad(il,k)=mad(il,k+1)+m(il,k+1)
        ENDIF
!        print*,"mad",il,k,mad(il,k)
    ENDDO
  ENDDO

!CR: erosion of each adiabatic ascent during its ascent

!Computation of erosion coefficient beta_coef
  DO k = 1, nl
    DO il = 1, ncum
       IF ((k>=(icb(il)+1)) .AND. (k<=inb(il)) .AND. (mlim(il,k).gt.0.)) THEN     
!          print*,"beta_coef",il,k,icb(il),inb(il),buoy(il,k),tv(il,k),wlim(il,k),wlim(il,k+1)
          beta_coef(il,k)=RG*coef_peel*buoy(il,k)/tv(il,k)/((wlim(il,k)+wlim(il,k+1))/2.)**2
       ELSE
          beta_coef(il,k)=0.
       ENDIF
    ENDDO
  ENDDO

!  print*,"apres beta_coef"

  DO k = 1, nl
    DO il = 1, ncum

      IF ((k>=(icb(il)+1)) .AND. (k<=inb(il))) THEN

!        print*,"dz",il,k,tv(il, k-1)
        dz = (ph(il,k-1)-ph(il,k))/(p(il, k-1)/(rrd*tv(il, k-1))*RG)
        betalim(il,k)=betalim(il,k-1)*exp(-1.*beta_coef(il,k-1)*dz)
!        betalim(il,k)=betalim(il,k-1)*exp(-RG*coef_peel*buoy(il,k-1)/tv(il,k-1)/5.**2*dz)
!        print*,"me",il,k,mlim(il,k),buoy(il,k),wlim(il,k),mad(il,k)
        dz = (ph(il,k)-ph(il,k+1))/(p(il, k)/(rrd*tv(il, k))*RG)
!        me(il,k)=betalim(il,k)*(m(il,k)+RG*coef_peel*buoy(il,k)/tv(il,k)/((wlim(il,k)+wlim(il,k+1))/2.)**2*dz*mad(il,k))
        me(il,k)=betalim(il,k)*(m(il,k)+beta_coef(il,k)*dz*mad(il,k))
!        print*,"B/w2",il,k,RG*coef_peel*buoy(il,k)/tv(il,k)/((wlim(il,k)+wlim(il,k+1))/2.)**2*dz    
      
      END IF
        
!Modification of m
      m(il,k)=me(il,k) 
    END DO
  END DO
 
!  DO il = 1, ncum
!     dz = (ph(il,icb(il))-ph(il,icb(il)+1))/(p(il, icb(il))/(rrd*tv(il, icb(il)))*RG)
!     m(il,icb(il))=m(il,icb(il))+RG*coef_peel*buoy(il,icb(il))/tv(il,icb(il)) &
!                  /((wlim(il,icb(il))+wlim(il,icb(il)+1))/2.)**2*dz*mad(il,icb(il))
!     print*,"wlim(icb)",icb(il),wlim(il,icb(il)),m(il,icb(il))
!  ENDDO

!Verification sum(me)=sum(m)
  DO k = nl,1,-1
    DO il = 1, ncum
           med(il,k)=med(il,k+1)+m(il,k+1)
!           print*,"somme(me),somme(m)",il,k,icb(il),med(il,k),md(il,k),me(il,k),m(il,k),wlim(il,k)
    ENDDO
  ENDDO


  ENDIF !(iflag_mix_adiab)
!RC



  ! c 3. Compute final cloud base mass flux and set iflag to 3 if
  ! c    cloud base mass flux is exceedingly small and is decreasing (i.e. if
  ! c    the final mass flux (cbmflast) is greater than the target mass flux
  ! c    (cbmf)).

  DO il = 1, ncum
    cbmflast(il) = 0.
  END DO

  DO k = 1, nl
    DO il = 1, ncum
      IF (k>=icb(il) .AND. k<=inb(il)) THEN
          !IMpropo??      IF ((k.ge.(icb(il)+1)).and.(k.le.inb(il))) THEN
        cbmflast(il) = cbmflast(il) + m(il, k)
      END IF
    END DO
  END DO

  DO il = 1, ncum
    IF (cbmflast(il)<1.E-6 .AND. cbmflast(il)>=cbmf(il)) THEN
      iflag(il) = 3
    END IF
  END DO

  DO k = 1, nl
    DO il = 1, ncum
      IF (iflag(il)>=3) THEN
        m(il, k) = 0.
        sig(il, k) = 0.
        w0(il, k) = 0.
      END IF
    END DO
  END DO
  IF (prt_level>=20) PRINT *, 'cv3p1_param apres iflag'

  ! c 4. Introduce a correcting factor for coef, in order to obtain an
  ! effective
  ! c    sigdz larger in the present case (using cv3p1_closure) than in the
  ! old
  ! c    closure (using cv3_closure).
  IF (1==0) THEN
    DO il = 1, ncum
      ! c      coef(il) = 2.*coef(il)
      coef(il) = 5.*coef(il)
    END DO
    ! version CVS du ..2008
  ELSE
    IF (iflag_cvl_sigd==0) THEN
      ! test pour verifier qu on fait la meme chose qu avant: sid constant
      coef(1:ncum) = 1.
    ELSE
      coef(1:ncum) = min(2.*coef(1:ncum), 5.)
      coef(1:ncum) = max(2.*coef(1:ncum), 0.2)
    END IF
  END IF

  IF (prt_level>=20) PRINT *, 'cv3p1_param FIN'
  RETURN
END SUBROUTINE cv3p1_closure


