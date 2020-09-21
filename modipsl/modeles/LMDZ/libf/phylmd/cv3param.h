!------------------------------------------------------------
! Parameters for convectL, iflag_con=3:
! (includes - microphysical parameters,
!			- parameters that control the rate of approach
!               to quasi-equilibrium)
!			- noff & minorig (previously in input of convect1)
!------------------------------------------------------------

      logical ok_homo_tend
      logical ok_optim_yield
      logical ok_entrain
      logical ok_convstop
      logical ok_intermittent
      integer noff, minorig, nl, nlp, nlm
      integer cv_flag_feed
      integer flag_epKEorig,flag_wb
      real sigdz, spfac
      real pbcrit, ptcrit
      real elcrit, tlcrit
      real coef_peel
      real omtrain
      real dtovsh, dpbase, dttrig
      real dtcrit, tau, beta, alpha, alpha1
      real T_top_max
      real tau_stop, noconv_stop
      real wbmax
      real delta
      real betad

      COMMON /cv3param/ sigdz, spfac &
                      ,pbcrit, ptcrit &
                      ,elcrit, tlcrit &
                      ,coef_peel &
                      ,omtrain &
                      ,dtovsh, dpbase, dttrig &
                      ,dtcrit, tau, beta, alpha, alpha1 &
                      ,T_top_max &
                      ,tau_stop, noconv_stop &
                      ,wbmax &
                      ,delta, betad  &
                      ,flag_epKEorig &
                      ,flag_wb, cv_flag_feed &
                      ,noff, minorig, nl, nlp, nlm  &
                      ,ok_convstop, ok_intermittent &
                      ,ok_optim_yield &
                      ,ok_entrain &
                      ,ok_homo_tend
!$OMP THREADPRIVATE(/cv3param/)

