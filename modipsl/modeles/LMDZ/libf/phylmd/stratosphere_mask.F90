!
! $Id: stratosphere_mask.F90 3123 2017-12-13 13:16:38Z oboucher $
!
SUBROUTINE stratosphere_mask(missing_val, pphis, t_seri, pplay, xlat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! determination of tropopause height and temperature from gridded temperature data
!
! reference: Reichler, T., M. Dameris, and R. Sausen (GRL, 10.1029/2003GL018240, 2003)
! modified: 6/28/06 tjr
! adapted to LMDZ by C. Kleinschmitt (2016-02-15)
! committed to LMDz by O. Boucher (2016) with a mistake 
! mistake corrected by O. Boucher (2017-12-11)
!
! input:  temp(nlon,nlat,nlev)  3D-temperature field 
!         ps(nlon,nlat)         2D-surface pressure field
!         zs(nlon,nlat)         2D-surface height
!         nlon                  grid points in x
!         nlat                  grid points in y
!         pfull(nlon,nlat,nlev) full pressure levels in Pa
!         plimu                 upper limit for tropopause pressure
!         pliml                 lower limit for tropopause pressure 
!         gamma                 tropopause criterion, e.g. -0.002 K/m
!
! output: p_tropopause(klon)    tropopause pressure in Pa with missing values
!         t_tropopause(klon)    tropopause temperature in K with missing values
!         z_tropopause(klon)    tropopause height in m with missing values
!         stratomask            stratospheric mask withtout missing values
!         ifil                  # of undetermined values
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE dimphy
USE phys_local_var_mod, ONLY: stratomask
USE phys_local_var_mod, ONLY: p_tropopause, z_tropopause, t_tropopause
USE print_control_mod, ONLY: lunout, prt_level

IMPLICIT NONE

INCLUDE "YOMCST.h"

REAL, INTENT(IN)                       :: missing_val ! missing value, also XIOS
REAL,DIMENSION(klon),INTENT(IN)        :: pphis   ! Geopotentiel de surface
REAL,DIMENSION(klon,klev),INTENT(IN)   :: t_seri  ! Temperature
REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay   ! pression pour le mileu de chaque couche (en Pa)
REAL,DIMENSION(klon),INTENT(IN)        :: xlat    ! latitudes pour chaque point 

REAL, PARAMETER                        :: plimu=45000.
REAL, PARAMETER                        :: pliml=7500.
REAL, PARAMETER                        :: gamma=-0.002
LOGICAL, PARAMETER                     :: dofill=.true.
REAL,DIMENSION(klon)                   :: tp
REAL,DIMENSION(klev)                   :: t, p
INTEGER                                :: i, k, ifil
REAL                                   :: ptrp, ttrp, ztrp, psrf, zsrf, pi

pi     = 4.*ATAN(1.)

!--computing tropopause
DO i=1,klon
  DO k=1,klev 
    t(k)=t_seri(i,klev+1-k) 
    p(k)=pplay(i,klev+1-k)
  ENDDO
  psrf=pplay(i,1)
  zsrf=pphis(i)/RG           !--altitude de la surface
  call twmo(missing_val, klev, t, p, psrf, zsrf, plimu, pliml, gamma, ptrp, ttrp, ztrp)
  tp(i)=ptrp
  p_tropopause(i)=ptrp
  z_tropopause(i)=ztrp
  t_tropopause(i)=ttrp
ENDDO

!--filling holes in tp but not in p_tropopause
IF (dofill) THEN
  ifil=0
  DO i=1,klon
  IF (ABS(tp(i)/missing_val-1.0).LT.0.01) THEN
    !set missing values to very simple profile (neighbour averaging too expensive in LMDZ)
    tp(i)=50000.-20000.*cos(xlat(i)/360.*2.*pi)
    ifil=ifil+1
  ENDIF
  ENDDO
!
ENDIF
!
DO i=1, klon
DO k=1, klev
  IF (pplay(i,k).LT.tp(i)) THEN
    stratomask(i,k)=1.0
  ELSE
    stratomask(i,k)=0.0
  ENDIF
ENDDO
ENDDO

IF (ifil.GT.0 .AND. prt_level >5) THEN
  write(lunout,*)'Tropopause: number of undetermined values =', ifil
ENDIF

RETURN
END SUBROUTINE stratosphere_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! twmo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine twmo(missing_val, level, t, p, ps, zs, plimu, pliml, gamma, ptrp, ttrp, ztrp)

! reference: Reichler, T., M. Dameris, and R. Sausen (GRL, 10.1029/2003GL018240, 2003)

implicit none

include "YOMCST.h"

integer,intent(in)              :: level
real,intent(in)                 :: missing_val
real,intent(in),dimension(level):: t, p
real,intent(in)                 :: plimu, pliml, gamma, ps, zs
real,intent(out)                :: ptrp, ttrp, ztrp

real,parameter                  :: deltaz = 2000.0

real     :: faktor
real     :: pmk, pm, a, b, tm, dtdp, dtdz, dlnp, tdlnp
real     :: ag, bg, ptph, ttph, a0, b0
real     :: pm0, tm0, pmk0, dtdz0
real     :: p2km, asum, aquer
real     :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
integer  :: icount, jj, j

ptrp=missing_val
ttrp=missing_val
ztrp=missing_val

faktor = -RG/RD

do j=level,2,-1

   ! dt/dz
   pmk= 0.5 * (p(j-1)**rkappa+p(j)**rkappa)        ! p**k at half level
   pm = pmk**(1./rkappa)                    ! p at half level
   a = (t(j-1)-t(j))/(p(j-1)**rkappa-p(j)**rkappa)    ! dT/dp^k
   b = t(j)-(a*p(j)**rkappa)
   tm = a * pmk + b                    ! T at half level         
   dtdp = a * rkappa * pm**(rkappa-1.) ! dT/dp at half level
   dtdz = faktor*dtdp*pm/tm            ! dT/dz at half level

   ! dt/dz valid?
   if (j.eq.level)     go to 999                 ! no, start level, initialize first
   if (dtdz.le.gamma)  go to 999                 ! no, dt/dz < -2 K/km
   if (dtdz0.gt.gamma) go to 999                 ! no, dt/dz below > -2 K/km
   if (pm.gt.plimu)    go to 999                 ! no, pm too low

   ! dtdz is valid, calculate tropopause pressure ptph
   ag = (dtdz-dtdz0) / (pmk-pmk0)     
   bg = dtdz0 - (ag * pmk0)          
   ptph = exp(log((gamma-bg)/ag)/rkappa)

   ! calculate temperature at this ptph assuming linear gamma
   ! interpolation
   ttph = tm0
   ttph = ttph - (bg * log(pm0)  + ag * (pm0**rkappa) /rkappa) / faktor*t(j)
   ttph = ttph + (bg * log(ptph) + ag * (ptph**rkappa)/rkappa) / faktor*t(j)

   if (ptph.lt.pliml) go to 999 
   if (ptph.gt.plimu) go to 999     

   ! 2nd test: dtdz above 2 km must not exceed gamma
   p2km = ptph + deltaz*(pm/tm)*faktor      ! p at ptph + 2km
   asum = 0.0                               ! dtdz above
   icount = 0                               ! number of levels above

   ! test until apm < p2km
   do jj=j,2,-1

       pmk2 = .5 * (p(jj-1)**rkappa+p(jj)**rkappa)    ! p mean ^kappa
       pm2 = pmk2**(1/rkappa)                      ! p mean
       if(pm2.gt.ptph) go to 110                ! doesn't happen
       if(pm2.lt.p2km) go to 888                ! ptropo is valid

       a2 = (t(jj-1)-t(jj))                     ! a
       a2 = a2/(p(jj-1)**rkappa-p(jj)**rkappa)
       b2 = t(jj)-(a2*p(jj)**rkappa)               ! b
       tm2 = a2 * pmk2 + b2                     ! T mean
       dtdp2 = a2 * rkappa * (pm2**(rkappa-1))        ! dt/dp
       dtdz2 = faktor*dtdp2*pm2/tm2
       asum = asum+dtdz2
       icount = icount+1
       aquer = asum/float(icount)               ! dt/dz mean
   
       ! discard ptropo ?
        if (aquer.le.gamma) go to 999           ! dt/dz above < gamma

110 continue 
    enddo                   ! test next level

888 continue                    ! ptph is valid
    ptrp = ptph
    ttrp = ttph

! now calculate height of tropopause by integrating hypsometric equation
! from ps to ptrp
! linearly interpolate in p (results are identical to p^kappa)

jj = LEVEL                                       ! bottom
do while ((P(jj).gt.PS) .or. (T(jj).lt.100))     ! T must be valid too
  jj=jj-1
enddo

DLNP = log(PS/P(jj))                             ! from surface pressure
TM = T(jj)                                       ! take TM of lowest level (better: extrapolate)
TDLNP = TM*DLNP

do while ( (JJ.ge.2) .and. (PTRP.lt.P(jj-1)) ) 
  DLNP = log(P(jj)/P(jj-1))
  TM = 0.5 * (T(jj) + T(jj-1))
  TDLNP = TDLNP + TM*DLNP
  JJ=JJ-1
enddo

DLNP = log(P(jj)/PTRP)                           ! up to tropopause pressure
TM = 0.5 * (T(jj) + TTRP)                        ! use TTRP to get TM of this level
TDLNP = TDLNP + TM*DLNP

ZTRP = ZS + TDLNP*RD/RG

!!if (ZTRP .lt. 0) then
!!  print*,'ZTRP=',ZTRP
!!  print*,'PS=',PS
!!  print*,'P=',P
!!  print*,'T=',T
!!  print*,'ZS=',ZS
!!  stop
!!endif

return

999 continue                    ! continue search at next higher level 
    tm0 = tm
    pm0 = pm
    pmk0 = pmk
    dtdz0  = dtdz
    a0 = a
    b0 = b

enddo 

! no tropopouse found
return
end subroutine twmo
