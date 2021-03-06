MODULE SUSTAONL_MOD
CONTAINS
SUBROUTINE SUSTAONL(KMEDIAP,KRESTM)

!**** *SUSTAONL * - Routine to initialize parallel environment

!     Purpose.
!     --------
!           Initialize D%NSTA and D%NONL.
!           Calculation of distribution of grid points to processors :
!           Splitting of grid in B direction

!**   Interface.
!     ----------
!        *CALL* *SUSTAONL *

!        Explicit arguments : KMEDIAP - mean number of grid points per PE
!        -------------------- KRESTM  - number of PEs with one extra point

!        Implicit arguments :
!        --------------------


!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-10-01
!        Modified 98-08-10 by K. YESSAD: removal of LRPOLE option.
!          - removal of LRPOLE in YOMCT0.
!          - removal of code under LRPOLE.
!        Modified 98-12-04 C. Fischer: merge with SUESTAONL (Aladin)
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
!USE MPL_MODULE      ! MPL 4.12.08

USE TPM_GEN
USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_DISTR

USE SET2PE_MOD
USE ABORT_TRANS_MOD
USE EQ_REGIONS_MOD

IMPLICIT NONE


!     DUMMY 
INTEGER(KIND=JPIM),INTENT(IN) :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(IN) :: KRESTM

!     LOCAL

INTEGER(KIND=JPIM) :: IXPTLAT(R%NDGL), ILSTPTLAT(R%NDGL)
INTEGER(KIND=JPIM) :: ICHK(R%NDLON,R%NDGL), ICOMBUF(R%NDGL*N_REGIONS_EW*2)
INTEGER(KIND=JPIM) :: I1, I2, IBUFLEN, IDGLG, IDWIDE,&
             &IGL, IGL1, IGL2, IGLOFF, IGPTA, IGPTOT, &
             &IGPTPRSETS, IGPTS, IGPTSP, ILEN, ILRECV, &
             &ILSEND, INPLAT, INXLAT, IPART,  IPOS, &
             &IPROCB, IPTSRE, IRECV, IPE, &
             &IREST, ISEND, ITAG, JA, JB, JGL, JL, JNPTSRE

LOGICAL :: LLABORT, LLALLAT
LOGICAL :: LLP1,LLP2

REAL(KIND=JPRB) ::  ZLAT, ZLAT1
REAL(KIND=JPRB) :: ZDIVID(R%NDGL),ZXPTLAT(R%NDGL)

!      -----------------------------------------------------------------

LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1

IDWIDE  = R%NDGL/2
IBUFLEN = R%NDGL*N_REGIONS_EW*2
IDGLG   = R%NDGL

I1 = MAX(   1,D%NFRSTLAT(MY_REGION_NS)-D%NFRSTLOFF)
I2 = MIN(IDGLG,D%NLSTLAT (MY_REGION_NS)-D%NFRSTLOFF)

ILEN = D%NLSTLAT(MY_REGION_NS) - D%NFRSTLAT(MY_REGION_NS)+1

IGPTPRSETS = SUM(G%NLOEN(1:D%NFRSTLAT(MY_REGION_NS)-1))

IGPTOT = SUM(G%NLOEN(1:R%NDGL))

IF (D%LSPLIT) THEN
  IF( LEQ_REGIONS )THEN
    IPE=0
    IGPTA=0
    DO JA=1,MY_REGION_NS-1
      DO JB=1,N_REGIONS(JA)
        IPE=IPE+1
        IF( IPE <= KRESTM .OR. KRESTM  ==  0)THEN
          IGPTA  = IGPTA + KMEDIAP
        ELSE
          IGPTA  = IGPTA + (KMEDIAP-1)
        ENDIF
      ENDDO
    ENDDO
    IGPTS=0
    DO JB=1,N_REGIONS(MY_REGION_NS)
      IPE=IPE+1
      IF( IPE <= KRESTM .OR. KRESTM  ==  0 )THEN
        IGPTS = IGPTS + KMEDIAP
      ELSE
        IGPTS = IGPTS + (KMEDIAP-1)
      ENDIF
    ENDDO
  ELSE
    IF (MY_REGION_NS <= KRESTM.OR.KRESTM == 0) THEN
      IGPTS = KMEDIAP
      IGPTA = KMEDIAP*(MY_REGION_NS-1)
    ELSE
      IGPTS = KMEDIAP-1
      IGPTA = KMEDIAP*KRESTM+IGPTS*(MY_REGION_NS-1-KRESTM)
    ENDIF
  ENDIF
ELSE
  IGPTA = IGPTPRSETS
  IGPTS = SUM(G%NLOEN(D%NFRSTLAT(MY_REGION_NS):D%NLSTLAT(MY_REGION_NS)))
ENDIF

IGPTSP = IGPTS/N_REGIONS(MY_REGION_NS)
IREST = IGPTS-N_REGIONS(MY_REGION_NS)*IGPTSP
IXPTLAT(1) = IGPTA-IGPTPRSETS+1
ZXPTLAT(1) = REAL(IXPTLAT(1))
ILSTPTLAT(1) = G%NLOEN(D%NFRSTLAT(MY_REGION_NS))
INPLAT = G%NLOEN(D%NFRSTLAT(MY_REGION_NS))-IXPTLAT(1)+1
DO JGL=2,ILEN
  IXPTLAT(JGL) = 1
  ZXPTLAT(JGL) = 1.0_JPRB
  ILSTPTLAT(JGL) =  G%NLOEN(D%NFRSTLAT(MY_REGION_NS)+JGL-1)
  INPLAT = INPLAT+G%NLOEN(D%NFRSTLAT(MY_REGION_NS)+JGL-1)
ENDDO
ILSTPTLAT(ILEN) = G%NLOEN(D%NLSTLAT(MY_REGION_NS))-INPLAT+IGPTS

DO JB=1,N_REGIONS_EW
  DO JGL=1,R%NDGL+N_REGIONS_NS-1
    D%NSTA(JGL,JB) = 0
    D%NONL(JGL,JB) = 0
  ENDDO
ENDDO


!  grid point decomposition
!  ---------------------------------------
LLALLAT = (N_REGIONS_NS == 1)
DO JGL=1,ILEN
  ZDIVID(JGL)=REAL(G%NLOEN(D%NFRSTLAT(MY_REGION_NS)+JGL-1),JPRB)
ENDDO
DO JB=1,N_REGIONS(MY_REGION_NS)

  IF (JB <= IREST) THEN
    IPTSRE = IGPTSP+1
  ELSE
    IPTSRE = IGPTSP
  ENDIF

  IPART=0
  DO JNPTSRE=1,IPTSRE
    ZLAT  = 1._JPRB
    ZLAT1 = 1._JPRB
    IF (MY_REGION_NS <= D%NAPSETS .AND.(IPART /= 2.OR.LLALLAT)) THEN
!cdir novector
      DO JGL=1,ILEN
        IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
          ZLAT1  = (ZXPTLAT(JGL)-1.0_JPRB)/ZDIVID(JGL)
          ZLAT   = MIN(ZLAT1,ZLAT)
          INXLAT = JGL
          IPART  = 1
          EXIT
        ENDIF
      ENDDO
    ELSEIF (MY_REGION_NS > N_REGIONS_NS-D%NAPSETS.AND.(IPART /= 1.OR.LLALLAT)) THEN
!cdir novector
      DO JGL=1,ILEN
        IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
          ZLAT1  = (ZXPTLAT(JGL)-1.0_JPRB)/ZDIVID(JGL)
          ZLAT   = MIN(ZLAT1,ZLAT)
          INXLAT = JGL
          IPART  = 2
          EXIT
        ENDIF
      ENDDO
    ELSE
!cdir novector
      DO JGL=1,ILEN
        IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
          ZLAT1 = (ZXPTLAT(JGL)-1.0_JPRB)/ZDIVID(JGL)
          IF (ZLAT1 < ZLAT) THEN
            ZLAT   = ZLAT1
            INXLAT = JGL
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    IF (INXLAT >= I1 .AND. INXLAT <= I2) THEN
      IF (D%NSTA(D%NPTRFLOFF+INXLAT,JB) == 0) THEN
        D%NSTA(D%NPTRFLOFF+INXLAT,JB) = IXPTLAT(INXLAT)
      ENDIF
      D%NONL(D%NPTRFLOFF+INXLAT,JB) = D%NONL(D%NPTRFLOFF+INXLAT,JB)+1
    ENDIF
    IXPTLAT(INXLAT) = IXPTLAT(INXLAT)+1
    ZXPTLAT(INXLAT) = REAL(IXPTLAT(INXLAT),JPRB)
  ENDDO
ENDDO


! Exchange local partitioning info to produce global view
!

IF( NPROC > 1 )THEN

  IF( LEQ_REGIONS )THEN

    ITAG = MTAGPART
    IPOS = 0
    DO JGL=1,D%NLSTLAT(MY_REGION_NS)-D%NFRSTLAT(MY_REGION_NS)+1
      IPOS = IPOS+1
      ICOMBUF(IPOS) = D%NSTA(D%NPTRFLOFF+JGL,MY_REGION_EW)
      IPOS = IPOS+1
      ICOMBUF(IPOS) = D%NONL(D%NPTRFLOFF+JGL,MY_REGION_EW)
    ENDDO
    IF( IPOS > IBUFLEN )THEN
      CALL ABORT_TRANS(' SUSTAONL: SEND BUFFER TOO SMALL FOR GLOBAL INFO')
    ENDIF
    ILSEND = IPOS

    DO JA=1,N_REGIONS_NS
      DO JB=1,N_REGIONS(JA)
        CALL SET2PE(ISEND,JA,JB,0,0)
        IF(ISEND /= MYPROC) THEN
!         CALL MPL_SEND(ICOMBUF(1:ILSEND),KDEST=NPRCIDS(ISEND),KTAG=ITAG, &
!          &   CDSTRING='SUSTAONL:')
!         MPL 4.12.08
          CALL ABOR1(' SUSTAONL: JUSTE APRES MPL_SEND')
        ENDIF
      ENDDO
    ENDDO

    DO JA=1,N_REGIONS_NS
      IGL1 = D%NFRSTLAT(JA)
      IGL2 = D%NLSTLAT(JA)
      DO JB=1,N_REGIONS(JA)
        CALL SET2PE(IRECV,JA,JB,0,0)
        IF(IRECV /= MYPROC) THEN
          ILEN = (D%NLSTLAT(JA)-D%NFRSTLAT(JA)+1)*2
!         CALL MPL_RECV(ICOMBUF(1:ILEN),KSOURCE=NPRCIDS(IRECV),KTAG=ITAG, &
!          & KOUNT=ILRECV,CDSTRING='SUSTAONL:')
!         MPL 4.12.08
          CALL ABOR1(' SUSTAONL: JUSTE APRES MPL_RCV')
          IPOS = 0
          DO JGL=IGL1,IGL2
            IGL = D%NPTRFRSTLAT(JA)+JGL-IGL1
            IPOS = IPOS+1
            D%NSTA(IGL,JB) = ICOMBUF(IPOS)
            IPOS = IPOS+1
            D%NONL(IGL,JB) = ICOMBUF(IPOS)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

  ELSE

    ITAG = MTAGPART
    IPOS = 0
    DO JB=1,N_REGIONS(MY_REGION_NS)
      DO JGL=1,D%NLSTLAT(MY_REGION_NS)-D%NFRSTLAT(MY_REGION_NS)+1
        IPOS = IPOS+1
        ICOMBUF(IPOS) = D%NSTA(D%NPTRFLOFF+JGL,JB)
        IPOS = IPOS+1
        ICOMBUF(IPOS) = D%NONL(D%NPTRFLOFF+JGL,JB)
      ENDDO
    ENDDO
    IF( IPOS > IBUFLEN )THEN
      CALL ABORT_TRANS(' SUSTAONL: SEND BUFFER TOO SMALL FOR GLOBAL INFO')
    ENDIF
    ILSEND = IPOS
    DO JA=1,N_REGIONS_NS
      CALL SET2PE(ISEND,JA,MY_REGION_EW,0,0)
      IF(ISEND /= MYPROC) THEN
!       CALL MPL_SEND(ICOMBUF(1:ILSEND),KDEST=NPRCIDS(ISEND),KTAG=ITAG, &
!        &   CDSTRING='SUSTAONL:')
!         MPL 4.12.08
          CALL ABOR1(' SUSTAONL: JUSTE APRES MPL_SEND')
      ENDIF
    ENDDO

    DO JA=1,N_REGIONS_NS
      CALL SET2PE(IRECV,JA,MY_REGION_EW,0,0)
      IF(IRECV /= MYPROC) THEN
        ILEN = (D%NLSTLAT(JA)-D%NFRSTLAT(JA)+1)*N_REGIONS(JA)*2
!       CALL MPL_RECV(ICOMBUF(1:ILEN),KSOURCE=NPRCIDS(IRECV),KTAG=ITAG, &
!        & KOUNT=ILRECV,CDSTRING='SUSTAONL:')
!         MPL 4.12.08
          CALL ABOR1(' SUSTAONL: JUSTE APRES MPL_RCV')
        IGL1 = D%NFRSTLAT(JA)
        IGL2 = D%NLSTLAT(JA)
        IPOS = 0
        DO JB=1,N_REGIONS(JA)
          DO JGL=IGL1,IGL2
            IGL = D%NPTRFRSTLAT(JA)+JGL-IGL1
            IPOS = IPOS+1
            D%NSTA(IGL,JB) = ICOMBUF(IPOS)
            IPOS = IPOS+1
            D%NONL(IGL,JB) = ICOMBUF(IPOS)
          ENDDO
        ENDDO
      ENDIF
    ENDDO

  ENDIF

ENDIF

! Confirm consistency of global partitioning, specifically testing for
! multiple assignments of same grid point and unassigned grid points

LLABORT = .FALSE.
DO JGL=1,R%NDGL
  DO JL=1,G%NLOEN(JGL)
    ICHK(JL,JGL) = 1
  ENDDO
ENDDO
DO JA=1,N_REGIONS_NS
  IGLOFF = D%NPTRFRSTLAT(JA)
  DO JB=1,N_REGIONS(JA)
    IGL1 = D%NFRSTLAT(JA)
    IGL2 = D%NLSTLAT(JA)
    DO JGL=IGL1,IGL2
      IGL = IGLOFF+JGL-IGL1
      DO JL=D%NSTA(IGL,JB),D%NSTA(IGL,JB)+D%NONL(IGL,JB)-1
        IF( ICHK(JL,JGL) /= 1 )THEN
          WRITE(NOUT,'(" SUSTAONL : seta=",i4," setb=",i4,&
           &" row=",I4," sta=",I4," INVALID GRID POINT")')&
           &JA,JB,JGL,JL
          WRITE(0,'(" SUSTAONL : seta=",i4," setb=",i4,&
           &" ROW=",I4," sta=",I4," INVALID GRID POINT")')&
           &JA,JB,JGL,JL
          LLABORT = .TRUE.
        ENDIF
        ICHK(JL,JGL) = 2
      ENDDO
    ENDDO
  ENDDO
ENDDO
DO JGL=1,R%NDGL
  DO JL=1,G%NLOEN(JGL)
    IF( ICHK(JL,JGL) /= 2 )THEN
      WRITE(NOUT,'(" SUSTAONL : row=",i4," sta=",i4,&
       &" GRID POINT NOT ASSIGNED")') JGL,JL
      LLABORT = .TRUE.
    ENDIF
  ENDDO
ENDDO
IF( LLABORT )THEN
  WRITE(NOUT,'(" SUSTAONL : inconsistent partitioning")')
  CALL ABORT_TRANS(' SUSTAONL: inconsistent partitioning')
ENDIF


IF (LLP1) THEN
  WRITE(UNIT=NOUT,FMT='('' OUTPUT FROM ROUTINE SUSTAONL '')')
  WRITE(UNIT=NOUT,FMT='('' '')')
  WRITE(UNIT=NOUT,FMT='('' PARTITIONING INFORMATION '')')
  WRITE(UNIT=NOUT,FMT='('' '')')
  IPROCB = MIN(32,N_REGIONS_EW)
  WRITE(UNIT=NOUT,FMT='(17X," SETB=",32(1X,I3))') (JB,JB=1,IPROCB)
  DO JA=1,N_REGIONS_NS
    IPROCB = MIN(32,N_REGIONS(JA))
    WRITE(UNIT=NOUT,FMT='('' '')')
    IGLOFF = D%NPTRFRSTLAT(JA)
    IGL1 = D%NFRSTLAT(JA)
    IGL2 = D%NLSTLAT(JA)
    DO JGL=IGL1,IGL2
      IGL=IGLOFF+JGL-IGL1
      WRITE(UNIT=NOUT,FMT='(" SETA=",I3," LAT=",I3," NSTA=",&
       &32(1X,I3))') JA,JGL,(D%NSTA(IGL,JB),JB=1,IPROCB)
      WRITE(UNIT=NOUT,FMT='(" SETA=",I3," LAT=",I3," D%NONL=",&
       &32(1X,I3))') JA,JGL,(D%NONL(IGL,JB),JB=1,IPROCB)
      WRITE(UNIT=NOUT,FMT='('' '')')
    ENDDO
    WRITE(UNIT=NOUT,FMT='('' '')')
  ENDDO
  WRITE(UNIT=NOUT,FMT='('' '')')
  WRITE(UNIT=NOUT,FMT='('' '')')
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE SUSTAONL
END MODULE SUSTAONL_MOD

