! $Id: $

PROGRAM gcm

#ifdef CPP_IOIPSL
  USE IOIPSL
#endif

  USE mod_const_mpi, ONLY: init_const_mpi
  USE parallel_lmdz
  USE infotrac
!#ifdef CPP_PHYS
!  USE mod_interface_dyn_phys, ONLY: init_interface_dyn_phys
!#endif
  USE mod_hallo
  USE Bands
  USE filtreg_mod
  USE control_mod

#ifdef CPP_PHYS
  USE iniphysiq_mod, ONLY: iniphysiq
#endif
  USE comconst_mod, ONLY: cpp, daysec, dtphys, dtvr, g, r, rad
  USE logic_mod ! all of it, because of copyin clause when calling leapfrog
  USE temps_mod, ONLY: calend,start_time,annee_ref,day_ref, &
                       itau_dyn,itau_phy,day_ini,jD_ref,jH_ref,day_end, &
                       dt,hour_ini,itaufin

  IMPLICIT NONE

  !      ......   Version  du 10/01/98    ..........

  !             avec  coordonnees  verticales hybrides 
  !   avec nouveaux operat. dissipation * ( gradiv2,divgrad2,nxgraro2 )

  !=======================================================================
  !
  !   Auteur:  P. Le Van /L. Fairhead/F.Hourdin
  !   -------
  !
  !   Objet:
  !   ------
  !
  !   GCM LMD nouvelle grille
  !
  !=======================================================================
  !
  !  ... Dans inigeom , nouveaux calculs pour les elongations  cu , cv
  !      et possibilite d'appeler une fonction f(y)  a derivee tangente
  !      hyperbolique a la  place de la fonction a derivee sinusoidale.
  !  ... Possibilite de choisir le schema pour l'advection de
  !        q  , en modifiant iadv dans traceur.def  (MAF,10/02) .
  !
  !      Pour Van-Leer + Vapeur d'eau saturee, iadv(1)=4. (F.Codron,10/99)
  !      Pour Van-Leer iadv=10
  !
  !-----------------------------------------------------------------------
  !   Declarations:
  !   -------------
  include "dimensions.h"
  include "paramet.h"
  include "comdissnew.h"
  include "comgeom.h"
  include "description.h"
  include "iniprint.h"
  include "tracstoke.h"


  REAL zdtvr

  !   variables dynamiques
  REAL,ALLOCATABLE,SAVE  :: vcov(:,:),ucov(:,:) ! vents covariants
  REAL,ALLOCATABLE,SAVE  :: teta(:,:)     ! temperature potentielle 
  REAL, ALLOCATABLE,SAVE :: q(:,:,:)      ! champs advectes
  REAL,ALLOCATABLE,SAVE  :: ps(:)         ! pression  au sol
  !      REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
  REAL,ALLOCATABLE,SAVE  :: masse(:,:)    ! masse d'air
  REAL,ALLOCATABLE,SAVE  :: phis(:)       ! geopotentiel au sol
  !      REAL phi(ip1jmp1,llm)                  ! geopotentiel
  !      REAL w(ip1jmp1,llm)                    ! vitesse verticale

  ! variables dynamiques intermediaire pour le transport

  !   variables pour le fichier histoire
  REAL dtav      ! intervalle de temps elementaire

  REAL time_0

  LOGICAL lafin

  real time_step, t_wrt, t_ops

  !+jld variables test conservation energie
  !      REAL ecin(ip1jmp1,llm),ecin0(ip1jmp1,llm)
  !     Tendance de la temp. potentiel d (theta)/ d t due a la 
  !     tansformation d'energie cinetique en energie thermique
  !     cree par la dissipation
  !      REAL dhecdt(ip1jmp1,llm)
  !      REAL vcont(ip1jm,llm),ucont(ip1jmp1,llm)
  !      REAL      d_h_vcol, d_qt, d_qw, d_ql, d_ec
  !      CHARACTER (len=15) :: ztit
  !-jld 


  character (len=80) :: dynhist_file, dynhistave_file
  character (len=20) :: modname
  character (len=80) :: abort_message
  ! locales pour gestion du temps
  INTEGER :: an, mois, jour
  REAL :: heure


  !-----------------------------------------------------------------------
  !   Initialisations:
  !   ----------------

  abort_message = 'last timestep reached'
  modname = 'gcm'
  descript = 'Run GCM LMDZ'
  lafin    = .FALSE.
  dynhist_file = 'dyn_hist'
  dynhistave_file = 'dyn_hist_ave'



  !----------------------------------------------------------------------
  !  lecture des fichiers gcm.def ou run.def
  !  ---------------------------------------
  !
  CALL conf_gcm( 99, .TRUE. )
  if (mod(iphysiq, iperiod) /= 0) call abort_gcm("conf_gcm", &
       "iphysiq must be a multiple of iperiod", 1)
  !
  !
  !------------------------------------
  !   Initialisation partie parallele
  !------------------------------------
  CALL init_const_mpi

  call init_parallel
  call Read_Distrib

!#ifdef CPP_PHYS
!  CALL Init_Phys_lmdz(iim,jjp1,llm,mpi_size,distrib_phys)
  !#endif
  !      CALL set_bands
  !#ifdef CPP_PHYS
!  CALL Init_interface_dyn_phys
!#endif
  CALL barrier

  CALL set_bands
  if (mpi_rank==0) call WriteBands
  call Set_Distrib(distrib_caldyn)

  !$OMP PARALLEL
  call Init_Mod_hallo
  !$OMP END PARALLEL

  !#ifdef CPP_PHYS
  !c$OMP PARALLEL
  !      call InitComgeomphy ! now done in iniphysiq
  !c$OMP END PARALLEL 
  !#endif

  !-----------------------------------------------------------------------
  !   Choix du calendrier
  !   -------------------

  !      calend = 'earth_365d'

#ifdef CPP_IOIPSL
  if (calend == 'earth_360d') then
     call ioconf_calendar('360d')
     write(lunout,*)'CALENDRIER CHOISI: Terrestre a 360 jours/an'
  else if (calend == 'earth_365d') then
     call ioconf_calendar('noleap')
     write(lunout,*)'CALENDRIER CHOISI: Terrestre a 365 jours/an'
  else if (calend == 'gregorian') then
     call ioconf_calendar('gregorian')
     write(lunout,*)'CALENDRIER CHOISI: Terrestre bissextile'
  else
     abort_message = 'Mauvais choix de calendrier'
     call abort_gcm(modname,abort_message,1)
  endif
#endif


  !-----------------------------------------------------------------------
  !   Initialisation des traceurs
  !   ---------------------------
  !  Choix du nombre de traceurs et du schema pour l'advection
  !  dans fichier traceur.def, par default ou via INCA
  call infotrac_init

  ! Allocation de la tableau q : champs advectes   
  ALLOCATE(ucov(ijb_u:ije_u,llm))
  ALLOCATE(vcov(ijb_v:ije_v,llm))
  ALLOCATE(teta(ijb_u:ije_u,llm))
  ALLOCATE(masse(ijb_u:ije_u,llm))
  ALLOCATE(ps(ijb_u:ije_u))
  ALLOCATE(phis(ijb_u:ije_u))
  ALLOCATE(q(ijb_u:ije_u,llm,nqtot))

  !-----------------------------------------------------------------------
  !   Lecture de l'etat initial :
  !   ---------------------------

  !  lecture du fichier start.nc
  if (read_start) then
     ! we still need to run iniacademic to initialize some
     ! constants & fields, if we run the 'newtonian' or 'SW' cases:
     if (iflag_phys.ne.1) then
        CALL iniacademic_loc(vcov,ucov,teta,q,masse,ps,phis,time_0)
     endif

     !        if (planet_type.eq."earth") then
     ! Load an Earth-format start file
     CALL dynetat0_loc("start.nc",vcov,ucov, &
          teta,q,masse,ps,phis, time_0)
     !        endif ! of if (planet_type.eq."earth")

     !       write(73,*) 'ucov',ucov
     !       write(74,*) 'vcov',vcov
     !       write(75,*) 'teta',teta
     !       write(76,*) 'ps',ps
     !       write(77,*) 'q',q

  endif ! of if (read_start)

  ! le cas echeant, creation d un etat initial
  IF (prt_level > 9) WRITE(lunout,*) &
       'GCM: AVANT iniacademic AVANT AVANT AVANT AVANT'
  if (.not.read_start) then
     CALL iniacademic_loc(vcov,ucov,teta,q,masse,ps,phis,time_0)
  endif

  !-----------------------------------------------------------------------
  !   Lecture des parametres de controle pour la simulation :
  !   -------------------------------------------------------
  !  on recalcule eventuellement le pas de temps

  IF(MOD(day_step,iperiod).NE.0) THEN
     abort_message =  &
          'Il faut choisir un nb de pas par jour multiple de iperiod'
     call abort_gcm(modname,abort_message,1)
  ENDIF

  IF(MOD(day_step,iphysiq).NE.0) THEN
     abort_message =  &
          'Il faut choisir un nb de pas par jour multiple de iphysiq'
     call abort_gcm(modname,abort_message,1)
  ENDIF

  zdtvr    = daysec/REAL(day_step)
  IF(dtvr.NE.zdtvr) THEN
     WRITE(lunout,*) &
          'WARNING!!! changement de pas de temps',dtvr,'>',zdtvr
  ENDIF

  !
  ! on remet le calendrier \`a zero si demande
  !
  IF (start_time /= starttime) then
     WRITE(lunout,*)' GCM: Attention l''heure de depart lue dans le' &
          ,' fichier restart ne correspond pas a celle lue dans le run.def'
     IF (raz_date == 1) then
        WRITE(lunout,*)'Je prends l''heure lue dans run.def'
        start_time = starttime
     ELSE
        WRITE(lunout,*)'Je m''arrete'
        CALL abort
     ENDIF
  ENDIF
  IF (raz_date == 1) THEN
     annee_ref = anneeref
     day_ref = dayref
     day_ini = dayref
     itau_dyn = 0
     itau_phy = 0
     time_0 = 0.
     write(lunout,*) &
          'GCM: On reinitialise a la date lue dans gcm.def'
  ELSE IF (annee_ref .ne. anneeref .or. day_ref .ne. dayref) THEN
     write(lunout,*) &
          'GCM: Attention les dates initiales lues dans le fichier'
     write(lunout,*) &
          ' restart ne correspondent pas a celles lues dans '
     write(lunout,*)' gcm.def'
     write(lunout,*)' annee_ref=',annee_ref," anneeref=",anneeref
     write(lunout,*)' day_ref=',day_ref," dayref=",dayref
     write(lunout,*)' Pas de remise a zero'
  ENDIF
  !      if (annee_ref .ne. anneeref .or. day_ref .ne. dayref) then
  !        write(lunout,*)
  !     .  'GCM: Attention les dates initiales lues dans le fichier'
  !        write(lunout,*)
  !     .  ' restart ne correspondent pas a celles lues dans '
  !        write(lunout,*)' gcm.def'
  !        write(lunout,*)' annee_ref=',annee_ref," anneeref=",anneeref
  !        write(lunout,*)' day_ref=',day_ref," dayref=",dayref
  !        if (raz_date .ne. 1) then
  !          write(lunout,*)
  !     .    'GCM: On garde les dates du fichier restart'
  !        else
  !          annee_ref = anneeref
  !          day_ref = dayref
  !          day_ini = dayref
  !          itau_dyn = 0
  !          itau_phy = 0
  !          time_0 = 0.
  !          write(lunout,*)
  !     .   'GCM: On reinitialise a la date lue dans gcm.def'
  !        endif
  !      ELSE
  !        raz_date = 0
  !      endif

#ifdef CPP_IOIPSL
  mois = 1
  heure = 0.
  call ymds2ju(annee_ref, mois, day_ref, heure, jD_ref)
  jH_ref = jD_ref - int(jD_ref)
  jD_ref = int(jD_ref)

  call ioconf_startdate(INT(jD_ref), jH_ref)

  write(lunout,*)'DEBUG'
  write(lunout,*)'annee_ref, mois, day_ref, heure, jD_ref'
  write(lunout,*)annee_ref, mois, day_ref, heure, jD_ref
  call ju2ymds(jD_ref+jH_ref,an, mois, jour, heure)
  write(lunout,*)'jD_ref+jH_ref,an, mois, jour, heure'
  write(lunout,*)jD_ref+jH_ref,an, mois, jour, heure
#else
  ! Ehouarn: we still need to define JD_ref and JH_ref
  ! and since we don't know how many days there are in a year
  ! we set JD_ref to 0 (this should be improved ...)
  jD_ref=0
  jH_ref=0
#endif

  if (iflag_phys.eq.1) then
     ! these initialisations have already been done (via iniacademic)
     ! if running in SW or Newtonian mode
     !-----------------------------------------------------------------------
     !   Initialisation des constantes dynamiques :
     !   ------------------------------------------
     dtvr = zdtvr
     CALL iniconst

     !-----------------------------------------------------------------------
     !   Initialisation de la geometrie :
     !   --------------------------------
     CALL inigeom

     !-----------------------------------------------------------------------
     !   Initialisation du filtre :
     !   --------------------------
     CALL inifilr
  endif ! of if (iflag_phys.eq.1)
  !
  !-----------------------------------------------------------------------
  !   Initialisation de la dissipation :
  !   ----------------------------------

  CALL inidissip( lstardis, nitergdiv, nitergrot, niterh   , &
       tetagdiv, tetagrot , tetatemp, vert_prof_dissip)

  !-----------------------------------------------------------------------
  !   Initialisation de la physique :
  !   -------------------------------
  IF ((iflag_phys==1).or.(iflag_phys>=100)) THEN
     ! Physics:
#ifdef CPP_PHYS
     CALL iniphysiq(iim,jjm,llm, &
          distrib_phys(mpi_rank),comm_lmdz, &
          daysec,day_ini,dtphys/nsplit_phys, &
          rlatu,rlatv,rlonu,rlonv,aire,cu,cv,rad,g,r,cpp, &
          iflag_phys)
#endif
  ENDIF ! of IF ((iflag_phys==1).or.(iflag_phys>=100))


  !-----------------------------------------------------------------------
  !   Initialisation des I/O :
  !   ------------------------


  if (nday>=0) then
     day_end = day_ini + nday
  else
     day_end = day_ini - nday/day_step
  endif

  WRITE(lunout,300)day_ini,day_end
300 FORMAT('1'/,15x,'run du jour',i7,2x,'au jour',i7//)

#ifdef CPP_IOIPSL
  call ju2ymds(jD_ref + day_ini - day_ref, an, mois, jour, heure)
  write (lunout,301)jour, mois, an
  call ju2ymds(jD_ref + day_end - day_ref, an, mois, jour, heure)
  write (lunout,302)jour, mois, an
301 FORMAT('1'/,15x,'run du ', i2,'/',i2,'/',i4)
302 FORMAT('1'/,15x,'    au ', i2,'/',i2,'/',i4)
#endif

  !      if (planet_type.eq."earth") then
  ! Write an Earth-format restart file
  CALL dynredem0_loc("restart.nc", day_end, phis)
  !      endif

  ecripar = .TRUE.

#ifdef CPP_IOIPSL
  time_step = zdtvr
     if (ok_dyn_ins) then
        ! initialize output file for instantaneous outputs
        ! t_ops = iecri * daysec ! do operations every t_ops
        t_ops =((1.0*iecri)/day_step) * daysec  
        t_wrt = daysec ! iecri * daysec ! write output every t_wrt
        CALL inithist_loc(day_ref,annee_ref,time_step, &
             t_ops,t_wrt)
     endif

     IF (ok_dyn_ave) THEN 
        ! initialize output file for averaged outputs
        t_ops = iperiod * time_step ! do operations every t_ops
        t_wrt = periodav * daysec   ! write output every t_wrt
        CALL initdynav_loc(day_ref,annee_ref,time_step,t_ops,t_wrt)
     END IF
  dtav = iperiod*dtvr/daysec
#endif
  ! #endif of #ifdef CPP_IOIPSL

  !  Choix des frequences de stokage pour le offline
  !      istdyn=day_step/4     ! stockage toutes les 6h=1jour/4
  !      istdyn=day_step/12     ! stockage toutes les 2h=1jour/12
  istdyn=day_step/4     ! stockage toutes les 6h=1jour/12
  istphy=istdyn/iphysiq     


  !
  !-----------------------------------------------------------------------
  !   Integration temporelle du modele :
  !   ----------------------------------

  !       write(78,*) 'ucov',ucov
  !       write(78,*) 'vcov',vcov
  !       write(78,*) 'teta',teta
  !       write(78,*) 'ps',ps
  !       write(78,*) 'q',q

  !!$OMP PARALLEL DEFAULT(SHARED) COPYIN(/temps/,/logici/,/logicl/)
  !$OMP PARALLEL DEFAULT(SHARED) &
  !     Copy all threadprivate variables in temps_mod
  !$OMP COPYIN(dt,jD_ref,jH_ref,start_time,hour_ini,day_ini,day_end) &
  !$OMP COPYIN(annee_ref,day_ref,itau_dyn,itau_phy,itaufin,calend) &
  !     Copy all threadprivate variables from logic_mod
  !$OMP COPYIN(purmats,forward,leapf,apphys,statcl,conser,apdiss,apdelq) &
  !$OMP COPYIN(saison,ecripar,fxyhypb,ysinus,read_start,ok_guide) &
  !$OMP COPYIN(ok_strato,ok_gradsfile,ok_limit,ok_etat0) &
  !$OMP COPYIN(iflag_phys,iflag_trac)
  CALL leapfrog_loc(ucov,vcov,teta,ps,masse,phis,q,time_0)
  !$OMP END PARALLEL

  !      OPEN(unit=5487,file='ok_lmdz',status='replace')
  !      WRITE(5487,*) 'ok_lmdz'
  !      CLOSE(5487)
END PROGRAM gcm
