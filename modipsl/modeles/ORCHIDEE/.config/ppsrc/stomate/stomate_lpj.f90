










! ================================================================================================================================
! MODULE       : stomate_lpj
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Main entry point for daily processes in STOMATE and LPJ (phenology, 
!! allocation, npp_calc, kill, turn, light, establish, crown, cover, lcchange)
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/tags/ORCHIDEE_2_0/ORCHIDEE/src_stomate/stomate_lpj.f90 $
!! $Date: 2019-11-29 15:15:17 +0100 (Fri, 29 Nov 2019) $
!! $Revision: 6367 $
!! \n
!_ ================================================================================================================================

MODULE stomate_lpj

  ! modules used:

  USE ioipsl_para
  USE xios_orchidee
  USE grid
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE lpj_constraints
  USE lpj_pftinout
  USE lpj_kill
  USE lpj_crown
  USE lpj_fire
  USE lpj_gap
  USE lpj_light
  USE lpj_establish
  USE lpj_cover
  USE stomate_prescribe
  USE stomate_phenology
  USE stomate_alloc
  USE stomate_npp
  USE stomate_turnover
  USE stomate_litter
  USE stomate_soilcarbon
  USE stomate_vmax
  USE stomate_lcchange
  USE stomate_woodharvest


  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC StomateLpj,StomateLpj_clear

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : StomateLpj_clear
!!
!>\BRIEF        Re-initialisation of variable
!!
!! DESCRIPTION  : This subroutine reinitializes variables. To be used if we want to relaunch 
!! ORCHIDEE but the routine is not used in current version.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE StomateLpj_clear

    CALL prescribe_clear
    CALL phenology_clear
    CALL npp_calc_clear
    CALL turn_clear
    CALL soilcarbon_clear
    CALL constraints_clear
    CALL establish_clear
    CALL fire_clear
    CALL gap_clear
    CALL light_clear
    CALL pftinout_clear
    CALL alloc_clear
  END SUBROUTINE StomateLpj_clear


!! ================================================================================================================================
!! SUBROUTINE   : StomateLPJ
!!
!>\BRIEF        Main entry point for daily processes in STOMATE and LPJ, structures the call sequence 
!!              to the different processes such as dispersion, establishment, competition and mortality of PFT's.
!! 
!! DESCRIPTION  : This routine is the main entry point to all processes calculated on a 
!! daily time step. Is mainly devoted to call the different STOMATE and LPJ routines 
!! depending of the ok_dgvm (is dynamic veg used) and lpj_constant_mortality (is background mortality used).
!! It also prepares the cumulative 
!! fluxes or pools (e.g TOTAL_M TOTAL_BM_LITTER etc...)
!!
!! This routine makes frequent use of "weekly", "monthly" and "long term" variables. Quotion is used because
!! by default "weekly" denotes 7 days, by default "monthly" denotes 20 days and by default "Long term" denotes
!! 3 years. dtslow refers to 24 hours (1 day).
!!
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): All variables related to stomate and required for LPJ dynamic vegetation mode.
!!
!! REFERENCE(S) : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudré, J. Ogeé, J. Polcher, P. Friedlingstein, P. Ciais, S. Sitch, 
!! and I. C. Prentice. 2005. A dynamic global vegetation model for studies of the coupled atmosphere-biosphere 
!! system. Global Biogeochemical Cycles 19:GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, I. C. Prentice, A. Arneth, A. Bondeau, W. Cramer, J. O. Kaplan, S. Levis, W. Lucht, 
!! M. T. Sykes, K. Thonicke, and S. Venevsky. 2003. Evaluation of ecosystem dynamics, plant geography and 
!! terrestrial carbon cycling in the LPJ dynamic global vegetation model. Global Change Biology 9:161-185.
!!
!! FLOWCHART    : Update with existing flowchart from N Viovy (Jan 19, 2012)
!! \n
!_ ================================================================================================================================
 
  SUBROUTINE StomateLpj (npts, dt_days, &
       neighbours, resolution, &
       clay, herbivores, &
       tsurf_daily, tsoil_daily, t2m_daily, t2m_min_daily, &
       litterhum_daily, soilhum_daily, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       gdd0_lastyear, precip_lastyear, &
       moiavail_month, moiavail_week, t2m_longterm, t2m_month, t2m_week, &
       tsoil_month, soilhum_month, &
       gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
       turnover_longterm, gpp_daily, &
       time_hum_min, maxfpc_lastyear, resp_maint_part, &
       PFTpresent, age, fireindex, firelitter, &
       leaf_age, leaf_frac, biomass, ind, adapted, regenerate, &
       senescence, when_growthinit, &
       litterpart, litter, dead_leaves, carbon, lignin_struc, &
       veget_cov_max, veget_cov_max_new, woodharvest, fraclut, npp_longterm, lm_lastyearmax, veget_lastlight, &
       everywhere, need_adjacent, RIP_time, &
       lai, rprof,npp_daily, turnover_daily, turnover_time,&
       control_moist, control_temp, soilcarbon_input, &
       co2_to_bm, co2_fire, resp_hetero, resp_hetero_litter, resp_hetero_soil, resp_maint, resp_growth, &
       height, deadleaf_cover, vcmax, &
       bm_to_litter, &
       prod10,prod100,flux10, flux100, &
       convflux,cflux_prod10,cflux_prod100, &
       prod10_harvest,prod100_harvest,flux10_harvest, flux100_harvest, &
       convflux_harvest,cflux_prod10_harvest,cflux_prod100_harvest, woodharvestpft, &
       convfluxpft, fDeforestToProduct, fLulccResidue,fHarvestToProduct, &
       harvest_above, carb_mass_total, fpc_max, MatrixA, &
       Tseason, Tmin_spring_time, begin_leaves, onset_date)
    
  !! 0. Variable and parameter declaration

    !! 0.1 input

    INTEGER(i_std), INTENT(in)                                 :: npts                 !! Domain size (unitless)
    REAL(r_std), INTENT(in)                                    :: dt_days              !! Time step of Stomate (days)
    INTEGER(i_std), DIMENSION(npts,NbNeighb), INTENT(in)       :: neighbours           !! Indices of the 8 neighbours of each grid 
                                                                                       !! point [1=North and then clockwise] 
    REAL(r_std), DIMENSION(npts,2), INTENT(in)                 :: resolution           !! Resolution at each grid point (m)  
                                                                                       !! [1=E-W, 2=N-S] 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: clay                 !! Clay fraction (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: herbivores           !! Time constant of probability of a leaf to 
                                                                                       !! be eaten by a herbivore (days) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tsurf_daily          !! Daily surface temperatures (K)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: tsoil_daily          !! Daily soil temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_daily            !! Daily 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_min_daily        !! Daily minimum 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: litterhum_daily      !! Daily litter humidity (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: soilhum_daily        !! Daily soil humidity (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxmoiavail_lastyear !! Last year's maximum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: minmoiavail_lastyear !! Last year's minimum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: gdd0_lastyear        !! Last year's GDD0 (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: precip_lastyear      !! Lastyear's precipitation 
                                                                                       !! @tex $(mm year^{-1})$ @endtex
										       !! to determine if establishment possible
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_month       !! "Monthly" moisture availability (0 to 1, 
                                                                                       !! unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_week        !! "Weekly" moisture availability 
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_longterm         !! "Long term" 2 meter reference 
                                                                                       !! temperatures (K) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_month            !! "Monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_week             !! "Weekly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: Tseason              !! "seasonal" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: Tmin_spring_time     !! Number of days after begin_leaves (leaf onset) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: onset_date           !! Date in the year at when the leaves started to grow(begin_leaves)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: tsoil_month          !! "Monthly" soil temperatures (K)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: soilhum_month        !! "Monthly" soil humidity
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: gdd_m5_dormance      !! Growing degree days (K), threshold -5 deg 
                                                                                       !! C (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: gdd_from_growthinit  !! growing degree days, since growthinit for crops
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: gdd_midwinter        !! Growing degree days (K), since midwinter 
                                                                                       !! (for phenology) - this is written to the history files 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: ncd_dormance         !! Number of chilling days (days), since 
                                                                                       !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: ngd_minus5           !! Number of growing days (days), threshold 
                                                                                       !! -5 deg C (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)        :: turnover_longterm    !! "Long term" turnover rate  
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: gpp_daily            !! Daily gross primary productivity  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: time_hum_min         !! Time elapsed since strongest moisture 
                                                                                       !! availability (days) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxfpc_lastyear      !! Last year's maximum foliage projected
                                                                                       !! coverage for each natural PFT,
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)        :: resp_maint_part      !! Maintenance respiration of different 
                                                                                       !! plant parts  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: fpc_max              !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm),INTENT(in)                :: veget_cov_max_new    !! New "maximal" coverage fraction of a PFT 
    REAL(r_std), DIMENSION(npts),INTENT(in)                    :: woodharvest          !! Harvested wood biomass (gC m-2 yr-1)
    REAL(r_std),DIMENSION(npts, nlut),INTENT(in)               :: fraclut              !! Fraction of landuse tile

  !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: npp_daily            !! Net primary productivity 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(out) :: turnover_daily       !! Turnover rates 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: co2_to_bm            !! CO2 taken up from atmosphere when 
                                                                                       !! introducing a new PFT (introduced for 
                                                                                       !! carbon balance closure) 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: co2_fire             !! Carbon emitted into the atmosphere by 
                                                                                       !! fire (living and dead biomass)  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: resp_hetero          !! Heterotrophic respiration
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: resp_hetero_litter   !! Heterotrophic respiration from litter
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: resp_hetero_soil     !! Heterotrophic respiration from soil
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: resp_maint           !! Maintenance respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: resp_growth          !! Growth respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: deadleaf_cover       !! Fraction of soil covered by dead leaves 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vcmax                !! Maximum rate of carboxylation 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(out):: bm_to_litter      !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    LOGICAL, DIMENSION(npts,nvm), INTENT(out)                  :: begin_leaves         !! signal to start putting leaves on (true/false)

    !! 0.3 Modified variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: height               !! Height of vegetation (m) 
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(inout)          :: control_moist        !! Moisture control of heterotrophic 
                                                                                       !! respiration (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(inout)          :: control_temp         !! Temperature control of heterotrophic 
                                                                                       !! respiration, above and below 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: soilcarbon_input     !! Quantity of carbon going into carbon 
                                                                                       !! pools from litter decomposition  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: lai                  !! Leaf area index OF AN INDIVIDUAL PLANT,
										       !! where a PFT contains n indentical plants
										       !! i.e., using the mean individual approach 
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: rprof                !! Prescribed root depth (m) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: PFTpresent           !! Tab indicating which PFTs are present in 
                                                                                       !! each pixel 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: age                  !! Age (years)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: fireindex            !! Probability of fire (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: firelitter           !! Longer term litter above the ground that 
                                                                                       !! can be burned, @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_age             !! Leaf age (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_frac            !! Fraction of leaves in leaf age class, 
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass        !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: ind                  !! Density of individuals 
                                                                                       !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: adapted              !! Adaptation of PFT (killed if too cold) 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: regenerate           !! "Fitness": Winter sufficiently cold for 
                                                                                       !! PFT regeneration ? (0 to 1, unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: senescence           !! Flag for setting senescence stage (only 
                                                                                       !! for deciduous trees) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: when_growthinit      !! How many days ago was the beginning of 
                                                                                       !! the growing season (days) 
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)      :: litterpart           !! Fraction of litter above the ground 
                                                                                       !! belonging to different PFTs
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter     !! Metabolic and structural litter, above 
                                                                                       !! and below ground 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)      :: dead_leaves          !! Dead leaves on ground, per PFT, metabolic 
                                                                                       !! and structural,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon               !! Carbon pool: active, slow, or passive, 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)      :: lignin_struc         !! Ratio of Lignin/Carbon in structural 
                                                                                       !! litter, above and below ground,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: veget_cov_max        !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: npp_longterm         !! "Long term" mean yearly primary 
                                                                                       !! productivity 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: lm_lastyearmax       !! Last year's maximum leaf mass, for each 
                                                                                       !! PFT @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: veget_lastlight      !! Vegetation fractions (on ground) after 
                                                                                       !! last light competition  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: everywhere           !! Is the PFT everywhere in the grid box or 
                                                                                       !! very localized (after its introduction) 
                                                                                       !! (unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: need_adjacent        !! In order for this PFT to be introduced, 
                                                                                       !! does it have to be present in an 
                                                                                       !! adjacent grid box? 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: RIP_time             !! How much time ago was the PFT eliminated 
                                                                                       !! for the last time (y) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: turnover_time        !! Turnover_time of leaves for grasses 
                                                                                       !! (days)
    REAL(r_std),DIMENSION(npts,0:10), INTENT(inout)            :: prod10               !! Products remaining in the 10
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (10
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,0:100), INTENT(inout)           :: prod100              !! Products remaining in the 100 
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (100 
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,10), INTENT(inout)              :: flux10               !! Annual release from the 10
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,100), INTENT(inout)             :: flux100              !! Annual release from the 100 
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: convflux             !! Release during first year following land 
                                                                                       !! cover change @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: cflux_prod10         !! Total annual release from the 10 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: cflux_prod100        !! Total annual release from the 100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,0:10), INTENT(inout)            :: prod10_harvest       !! Products remaining in the 10
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (10
                                                                                       !! + 1 : input from year of wood harvest)
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,0:100), INTENT(inout)           :: prod100_harvest      !! Products remaining in the 100 
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (100 
                                                                                       !! + 1 : input from year of wood harvest)
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,10), INTENT(inout)              :: flux10_harvest       !! Annual release from the 10
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,100), INTENT(inout)             :: flux100_harvest      !! Annual release from the 100 
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: convflux_harvest     !! Release during first year following wood 
                                                                                       !! harvest @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: cflux_prod10_harvest !! Total annual release from the 10 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: cflux_prod100_harvest!! Total annual release from the 100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,nvm), INTENT(inout)             :: woodharvestpft       !! Harvested wood biomass (gC m-2 dt_stomate-1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: convfluxpft         !! Convflux per PFT
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: fDeforestToProduct  !! Deforested biomass into product pool due to anthropogenic 
                                                                                       !! land use change
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: fLulccResidue       !! Carbon mass flux into soil and litter due to anthropogenic land use or land cover change
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: fHarvestToProduct   !! Deforested biomass into product pool due to anthropogenic 
                                                                                       !! land use 

    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: harvest_above        !! Harvest above ground biomass for 
                                                                                       !! agriculture @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: carb_mass_total      !! Carbon Mass total (soil, litter, veg) 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL(r_std), DIMENSION(npts,nvm,nbpools,nbpools), INTENT(inout) :: MatrixA         !! Matrix containing the fluxes  
                                                                                       !! between the carbon pools
                                                                                       !! per sechiba time step 
                                                                                       !! @tex $(gC.m^2.day^{-1})$ @endtex

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: tot_bm_to_litter    !! Total conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: tot_live_biomass    !! Total living biomass  
                                                                                       !! @tex $(gC m{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: cOther              !! Carbon Mass in Vegetation Components 
                                                                                       !! other than Leaves, Stems and Roots
                                                                                       !! @tex $(gC m{-2})$ @endtex

    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements)           :: bm_alloc            !! Biomass increase, i.e. NPP per plant part 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: tot_turnover        !! Total turnover rate  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_litter_soil_carb!! Total soil and litter carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_litter_carb     !! Total litter carbon 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_soil_carb       !! Total soil carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts)                                :: sum_cLitterGrass    !! Carbon mass in litter on grass tiles
                                                                                       !! @tex $(kgC !m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts)                                :: sum_cLitterCrop     !! Carbon mass in litter on crop tiles
                                                                                       !! @tex $(kgC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                                :: sum_cSoilGrass      !! Carbon mass in soil on grass tiles 
                                                                                       !! @tex $(kgC !m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                                :: sum_cSoilCrop       !! Carbon mass in soil on crop tiles
                                                                                       !! @tex $(kgC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                                :: sum_cVegGrass       !! Carbon mass in vegetation on grass tiles 
                                                                                       !! @tex $(kgC !m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                                :: sum_cVegCrop        !! Carbon mass in vegetation on crop tiles
                                                                                       !! @tex $(kgC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                                :: sum_cLitterTree     !! Carbon mass in litter on tree tiles
                                                                                       !! @tex $(kgC !!m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts)                                :: sum_cSoilTree       !! Carbon mass in soil on tree tiles 
                                                                                       !! @tex $(kgC !m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                                :: sum_cVegTree        !! Carbon mass in vegetation on tree tiles 
                                                                                       !! @tex $(kgC !m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                                :: carb_mass_variation !! Carbon Mass variation  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: cn_ind              !! Crown area of individuals 
                                                                                       !! @tex $(m^{2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: woodmass_ind        !! Woodmass of individuals (gC) 
    REAL(r_std), DIMENSION(npts,nvm,nparts)                     :: f_alloc             !! Fraction that goes into plant part 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: avail_tree          !! Space availability for trees 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: avail_grass         !! Space availability for grasses 
                                                                                       !! (0 to 1, unitless) 
    INTEGER                                                     :: j,k
    REAL(r_std),DIMENSION(npts)                                 :: prod10_total        !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: prod100_total       !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: cflux_prod_total    !! Total flux from conflux and the 10/100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: prod10_harvest_total!! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: prod100_harvest_total!! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: cflux_prod_harvest_total!! Total flux from conflux and the 10/100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std),DIMENSION(npts,nvm)                             :: veget_cov_max_tmp   !! "Maximal" coverage fraction of a PFT  
                                                                                       !! (LAI-> infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(npts,nvm)                            :: mortality           !! Fraction of individual dying this time 
                                                                                       !! step (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: vartmp              !! Temporary variable used to add history
    REAL(r_std), DIMENSION(npts,nvm)                            :: histvar             !! History variables

    REAL(r_std), DIMENSION(npts,nlut)                           :: clitterlut          !! Litter carbon on landusetype4 (nlut)
    REAL(r_std), DIMENSION(npts,nlut)                           :: csoillut            !! Soil carbon on landusetype4 (nlut)
    REAL(r_std), DIMENSION(npts,nlut)                           :: cveglut             !! Carbon in vegetation on landusetype4 (nlut)
    REAL(r_std), DIMENSION(npts,nlut)                           :: lailut              !! LAI on landusetype4 (nlut)
    REAL(r_std), DIMENSION(npts,nlut)                           :: ralut               !! Autotrophic respiration on landusetype4 (nlut) 
    REAL(r_std), DIMENSION(npts,nlut)                           :: rhlut               !! Heterotrophic respiration on landusetype4 (nlut) 
    REAL(r_std), DIMENSION(npts,nlut)                           :: npplut              !! Net Primary Productivity on landusetype4 (nlut)    
    REAL(r_std), DIMENSION(npts,nlut)                           :: ctotfirelut         !! Fire CO2 emission on landusetype4 (nlut) 
    REAL(r_std), DIMENSION(npts,nlut)                           :: cproductlut
    REAL(r_std), DIMENSION(npts,nlut)                           :: flulccatmlut
    REAL(r_std), DIMENSION(npts,nlut)                           :: flulccproductlut
    REAL(r_std), DIMENSION(npts,nlut)                           :: flulccresiduelut
    REAL(r_std), DIMENSION(npts,ncarb)                          :: csoilpools          !! Diagnostics for carbon in soil pools
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering stomate_lpj'

  
  !! 1. Initializations
    
    !! 1.1 Initialize variables to zero
    co2_to_bm(:,:) = zero
    co2_fire(:,:) = zero
    npp_daily(:,:) = zero
    resp_maint(:,:) = zero
    resp_growth(:,:) = zero
    harvest_above(:) = zero
    bm_to_litter(:,:,:,:) = zero
    cn_ind(:,:) = zero
    woodmass_ind(:,:) = zero
    turnover_daily(:,:,:,:) = zero
    
    !! 1.2  Initialize variables to veget_cov_max
    veget_cov_max_tmp(:,:) = veget_cov_max(:,:)

    !! 1.3 Calculate some vegetation characteristics
    
    !! 1.3.1 Calculate some vegetation characteristics 
    !        Calculate cn_ind (individual crown mass) and individual height from
    !        state variables if running DGVM or dynamic mortality in static cover mode
    !??        Explain (maybe in the header once) why you mulitply with veget_cov_max in the DGVM
    !??        and why you don't multiply with veget_cov_max in stomate.
    IF ( ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN
       IF(ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon)) & 
                  *veget_cov_max(:,:))/ind(:,:)
          ENDWHERE
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF

       CALL crown (npts,  PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)
    ENDIF

    !! 1.3.2 Prescribe characteristics if the vegetation is not dynamic
    !        IF the DGVM is not activated, the density of individuals and their crown
    !        areas don't matter, but they should be defined for the case we switch on
    !        the DGVM afterwards. At the first call, if the DGVM is not activated, 
    !        impose a minimum biomass for prescribed PFTs and declare them present.
    CALL prescribe (npts, &
         veget_cov_max, dt_days, PFTpresent, everywhere, when_growthinit, &
         biomass, leaf_frac, ind, cn_ind, co2_to_bm)


  !! 2. Climatic constraints for PFT presence and regenerativeness

    !   Call this even when DGVM is not activated so that "adapted" and "regenerate"
    !   are kept up to date for the moment when the DGVM is activated.
    CALL constraints (npts, dt_days, &
         t2m_month, t2m_min_daily,when_growthinit, Tseason, &
         adapted, regenerate)

    
  !! 3. Determine introduction and elimination of PTS based on climate criteria
 
    IF ( ok_dgvm ) THEN
      
       !! 3.1 Calculate introduction and elimination
       CALL pftinout (npts, dt_days, adapted, regenerate, &
            neighbours, veget_cov_max, &
            biomass, ind, cn_ind, age, leaf_frac, npp_longterm, lm_lastyearmax, senescence, &
            PFTpresent, everywhere, when_growthinit, need_adjacent, RIP_time, &
            co2_to_bm, &
            avail_tree, avail_grass)

       !! 3.2 Reset attributes for eliminated PFTs.
       !     This also kills PFTs that had 0 leafmass during the last year. The message
       !     "... after pftinout" is misleading in this case.
       CALL kill (npts, 'pftinout  ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

       
       !! 3.3 Calculate woodmass of individual tree
       IF(ok_dgvm) THEN
          WHERE ((ind(:,:).GT.min_stomate))
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))*veget_cov_max(:,:))/ind(:,:)
          ENDWHERE
       ELSE
          WHERE ((ind(:,:).GT.min_stomate))
             woodmass_ind(:,:) =(biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF
       
       ! Calculate crown area and diameter for all PFTs (including the newly established)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)

    ENDIF
    
  !! 4. Phenology

    !! 4.1 Write values to history file
    !      Current values for ::when_growthinit 
    CALL xios_orchidee_send_field("WHEN_GROWTHINIT",when_growthinit)

    CALL histwrite_p (hist_id_stomate, 'WHEN_GROWTHINIT', itime, when_growthinit, npts*nvm, horipft_index)

    ! Set and write values for ::PFTpresent
    WHERE(PFTpresent)
       histvar=un
    ELSEWHERE
       histvar=zero
    ENDWHERE

    CALL xios_orchidee_send_field("PFTPRESENT",histvar)

    CALL histwrite_p (hist_id_stomate, 'PFTPRESENT', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_midwinter
    WHERE(gdd_midwinter.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_midwinter
    ENDWHERE

    CALL xios_orchidee_send_field("GDD_MIDWINTER",histvar)

    CALL histwrite_p (hist_id_stomate, 'GDD_MIDWINTER', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_m5_dormance
    WHERE(gdd_m5_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_m5_dormance
    ENDWHERE
    
    CALL xios_orchidee_send_field('GDD_M5_DORMANCE',histvar)
    CALL histwrite_p (hist_id_stomate, 'GDD_M5_DORMANCE', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for ncd_dormance
    WHERE(ncd_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=ncd_dormance
    ENDWHERE

    CALL xios_orchidee_send_field("NCD_DORMANCE",histvar)

    CALL histwrite_p (hist_id_stomate, 'NCD_DORMANCE', itime, histvar, npts*nvm, horipft_index)

    !! 4.2 Calculate phenology
    CALL phenology (npts, dt_days, PFTpresent, &
         veget_cov_max, &
         t2m_longterm, t2m_month, t2m_week, gpp_daily, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_month, moiavail_week, &
         gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
         senescence, time_hum_min, &
         biomass, leaf_frac, leaf_age, &
         when_growthinit, co2_to_bm, &
         begin_leaves)
    
  !! 5. Allocate C to different plant parts
    
    CALL alloc (npts, dt_days, &
         lai, veget_cov_max, senescence, when_growthinit, &
         moiavail_week, tsoil_month, soilhum_month, &
         biomass, age, leaf_age, leaf_frac, rprof, f_alloc)

  !! 6. NPP, maintenance and growth respiration

    !! 6.1 Calculate NPP and respiration terms
    CALL npp_calc (npts, dt_days, &
         PFTpresent, &
         t2m_daily, tsoil_daily, lai, rprof, &
         gpp_daily, f_alloc, bm_alloc, resp_maint_part,&
         biomass, leaf_age, leaf_frac, age, &
         resp_maint, resp_growth, npp_daily, co2_to_bm)

    !! 6.2 Kill slow growing PFTs in DGVM or STOMATE with constant mortality
    IF ( ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN
       CALL kill (npts, 'npp       ', lm_lastyearmax,  &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

       !! 6.2.1 Update wood biomass      
       !        For the DGVM
       IF(ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon)) & 
                  *veget_cov_max(:,:))/ind(:,:)
          ENDWHERE

       ! For all pixels with individuals
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF ! ok_dgvm

       !! 6.2.2 New crown area and maximum vegetation cover after growth
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind,&
            veget_cov_max, cn_ind, height)

    ENDIF ! ok_dgvm
    
  !! 7. fire

    !! 7.1. Burn PFTs
    CALL fire (npts, dt_days, litterpart, &
         litterhum_daily, t2m_daily, lignin_struc, veget_cov_max, &
         fireindex, firelitter, biomass, ind, &
         litter, dead_leaves, bm_to_litter, &
         co2_fire, MatrixA)

    !! 7.2 Kill PFTs in DGVM
    IF ( ok_dgvm ) THEN

       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'fire      ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

    ENDIF ! ok_dgvm
 
  !! 8. Tree mortality

    ! Does not depend on age, therefore does not change crown area.
    CALL gap (npts, dt_days, &
         npp_longterm, turnover_longterm, lm_lastyearmax, &
         PFTpresent, t2m_min_daily, Tmin_spring_time, &
         biomass, ind, bm_to_litter, mortality)


    IF ( ok_dgvm ) THEN

       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'gap       ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

    ENDIF

  !! 9. Leaf senescence, new lai and other turnover processes

    CALL turn (npts, dt_days, PFTpresent, &
         herbivores, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_week,  moiavail_month,t2m_longterm, t2m_month, t2m_week, veget_cov_max, &
         gdd_from_growthinit, leaf_age, leaf_frac, age, lai, biomass, &
         turnover_daily, senescence,turnover_time)

    !! 10. Light competition
    
    !! If not using constant mortality then kill with light competition
!    IF ( ok_dgvm .OR. .NOT.(lpj_gap_const_mort) ) THEN
    IF ( ok_dgvm ) THEN
 
       !! 10.1 Light competition
       CALL light (npts, dt_days, &
            veget_cov_max, fpc_max, PFTpresent, cn_ind, lai, maxfpc_lastyear, &
            lm_lastyearmax, ind, biomass, veget_lastlight, bm_to_litter, mortality)
       
       !! 10.2 Reset attributes for eliminated PFTs
       CALL kill (npts, 'light     ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_cov_max, bm_to_litter)

    ENDIF

    
  !! 11. Establishment of saplings
    
    IF ( ok_dgvm .OR. .NOT.lpj_gap_const_mort ) THEN

       !! 11.1 Establish new plants
       CALL establish (npts, dt_days, PFTpresent, regenerate, &
            neighbours, resolution, need_adjacent, herbivores, &
            precip_lastyear, gdd0_lastyear, lm_lastyearmax, &
            cn_ind, lai, avail_tree, avail_grass, npp_longterm, &
            leaf_age, leaf_frac, &
            ind, biomass, age, everywhere, co2_to_bm, veget_cov_max, woodmass_ind, &
            mortality, bm_to_litter)

       !! 11.2 Calculate new crown area (and maximum vegetation cover)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_cov_max, cn_ind, height)

    ENDIF

  !! 12. Calculate final LAI and vegetation cover
    
    CALL cover (npts, cn_ind, ind, biomass, &
         veget_cov_max, veget_cov_max_tmp, lai, &
         litter, carbon, turnover_daily, bm_to_litter, &
         co2_to_bm, co2_fire, resp_hetero, resp_hetero_litter, resp_hetero_soil, resp_maint, resp_growth, gpp_daily)

  !! 13. Update litter pools to account for harvest
 
    ! the whole litter stuff:
    !    litter update, lignin content, PFT parts, litter decay, 
    !    litter heterotrophic respiration, dead leaf soil cover.
    !    No vertical discretisation in the soil for litter decay.\n
    ! added by shilong for harvest
    IF(harvest_agri) THEN
       CALL harvest(npts, dt_days, veget_cov_max, &
            bm_to_litter, turnover_daily, &
            harvest_above)
    ENDIF

    !! 14. Land cover change if it is time to do so
    !! The flag do_now_lcchange is set in slowproc_main at the same time as the vegetation is read from file.
    !! The vegetation fractions are not updated yet and will be updated in the end of sechiba_main.
    IF (do_now_stomate_woodharvest) THEN 
       CALL woodharvest_main(npts, dt_days, veget_cov_max, &
            biomass, &
            flux10_harvest,flux100_harvest, prod10_harvest,prod100_harvest,&
            convflux_harvest,cflux_prod10_harvest,cflux_prod100_harvest,&
            woodharvest,woodharvestpft,fHarvestToProduct)
       do_now_stomate_woodharvest=.FALSE.
    ENDIF


    IF (do_now_stomate_lcchange) THEN
       CALL lcchange_main (npts, dt_days, veget_cov_max, veget_cov_max_new, &
            biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &
            co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, cn_ind,flux10,flux100, &
            prod10,prod100,convflux,cflux_prod10,cflux_prod100,leaf_frac,&
            npp_longterm, lm_lastyearmax, litter, carbon, &
            convfluxpft,  fDeforestToProduct, fLulccResidue)
       do_now_stomate_lcchange=.FALSE.

       ! Set the flag done_stomate_lcchange to be used in the end of sechiba_main to update the fractions.
       done_stomate_lcchange=.TRUE.
    ENDIF


    !! 15. Calculate vcmax 

    CALL vmax (npts, dt_days, &
         leaf_age, leaf_frac, &
         vcmax)

    !MM déplacement pour initialisation correcte des grandeurs cumulées :
    cflux_prod_total(:) = convflux(:) + cflux_prod10(:) + cflux_prod100(:)
    prod10_total(:)=SUM(prod10,dim=2)
    prod100_total(:)=SUM(prod100,dim=2)

    cflux_prod_harvest_total(:) = convflux_harvest(:) + cflux_prod10_harvest(:) + cflux_prod100_harvest(:)
    prod10_harvest_total(:)=SUM(prod10_harvest,dim=2)
    prod100_harvest_total(:)=SUM(prod100_harvest,dim=2)
    
  !! 16. Total heterotrophic respiration

    tot_soil_carb(:,:) = zero
    tot_litter_carb(:,:) = zero
    sum_cLitterGrass = zero
    sum_cLitterCrop = zero
    sum_cSoilGrass = zero
    sum_cSoilCrop = zero
    sum_cVegGrass = zero
    sum_cVegCrop = zero
    sum_cLitterTree = zero
    sum_cSoilTree = zero
    sum_cVegTree = zero

    DO j=2,nvm

       tot_litter_carb(:,j) = tot_litter_carb(:,j) + (litter(:,istructural,j,iabove,icarbon) + &
            &          litter(:,imetabolic,j,iabove,icarbon) + &
            &          litter(:,istructural,j,ibelow,icarbon) + litter(:,imetabolic,j,ibelow,icarbon))

      IF ((.NOT. is_tree(j))  .AND. natural(j)) THEN
         sum_cLitterGrass(:) = sum_cLitterGrass(:) + tot_litter_carb(:,j)*veget_cov_max(:,j)
      ELSE IF ((.NOT. is_tree(j))  .AND. (.NOT. natural(j)) ) THEN
         sum_cLitterCrop(:) = sum_cLitterCrop(:) + tot_litter_carb(:,j)*veget_cov_max(:,j)
      ELSE IF (is_tree(j)) THEN
         sum_cLitterTree(:) = sum_cLitterTree(:) + tot_litter_carb(:,j)*veget_cov_max(:,j)
      ENDIF

      tot_soil_carb(:,j) = tot_soil_carb(:,j) + (carbon(:,iactive,j) + &
           &          carbon(:,islow,j)+  carbon(:,ipassive,j))

      IF ((.NOT. is_tree(j)) .AND. natural(j)) THEN
         sum_cSoilGrass(:) = sum_cSoilGrass(:) + tot_soil_carb(:,j)*veget_cov_max(:,j)
      ELSE IF ((.NOT. is_tree(j))  .AND. (.NOT. natural(j)) ) THEN
         sum_cSoilCrop(:) = sum_cSoilCrop(:) + tot_soil_carb(:,j)*veget_cov_max(:,j)
      ELSE IF (is_tree(j)) THEN
         sum_cSoilTree(:) = sum_cSoilTree(:) + tot_soil_carb(:,j)*veget_cov_max(:,j)
      END IF
    ENDDO

    tot_litter_soil_carb(:,:) = tot_litter_carb(:,:) + tot_soil_carb(:,:)

!!$     DO k = 1, nelements ! Loop over # elements
!!$        tot_live_biomass(:,:,k) = biomass(:,:,ileaf,k) + biomass(:,:,isapabove,k) + biomass(:,:,isapbelow,k) +&
!!$             &                    biomass(:,:,iheartabove,k) + biomass(:,:,iheartbelow,k) + &
!!$             &                    biomass(:,:,iroot,k) + biomass(:,:,ifruit,k) + biomass(:,:,icarbres,k)
!!$    END DO ! Loop over # elements

    tot_live_biomass(:,:,:) = biomass(:,:,ileaf,:) + biomass(:,:,isapabove,:) + biomass(:,:,isapbelow,:) +&
             &                    biomass(:,:,iheartabove,:) + biomass(:,:,iheartbelow,:) + &
             &                    biomass(:,:,iroot,:) + biomass(:,:,ifruit,:) + biomass(:,:,icarbres,:)

    DO j= 1,nvm
       cOther(:,j,:) = biomass(:,j,ifruit,:) + biomass(:,j,icarbres,:)

       IF ((.NOT. is_tree(j))  .AND. natural(j)) THEN
          sum_cVegGrass(:) = sum_cVegGrass(:) + tot_live_biomass(:,j,icarbon)*veget_cov_max(:,j)
       ELSE IF ((.NOT. is_tree(j))  .AND. (.NOT. natural(j)) ) THEN
          sum_cVegCrop(:) = sum_cVegCrop(:) + tot_live_biomass(:,j,icarbon)*veget_cov_max(:,j)
       ELSE IF (is_tree(j)) THEN
          sum_cVegTree(:) = sum_cVegTree(:) + tot_live_biomass(:,j,icarbon)*veget_cov_max(:,j)
       ENDIF
    END DO


    tot_turnover(:,:,:) = turnover_daily(:,:,ileaf,:) + turnover_daily(:,:,isapabove,:) + &
         &         turnover_daily(:,:,isapbelow,:) + turnover_daily(:,:,iheartabove,:) + &
         &         turnover_daily(:,:,iheartbelow,:) + turnover_daily(:,:,iroot,:) + &
         &         turnover_daily(:,:,ifruit,:) + turnover_daily(:,:,icarbres,:)

    tot_bm_to_litter(:,:,:) = bm_to_litter(:,:,ileaf,:) + bm_to_litter(:,:,isapabove,:) +&
         &             bm_to_litter(:,:,isapbelow,:) + bm_to_litter(:,:,iheartbelow,:) +&
         &             bm_to_litter(:,:,iheartabove,:) + bm_to_litter(:,:,iroot,:) + &
         &             bm_to_litter(:,:,ifruit,:) + bm_to_litter(:,:,icarbres,:)

    carb_mass_variation(:)=-carb_mass_total(:)
    carb_mass_total(:)=SUM((tot_live_biomass(:,:,icarbon)+tot_litter_carb+tot_soil_carb)*veget_cov_max,dim=2) + &
         &                 (prod10_total + prod100_total) +  (prod10_harvest_total + prod100_harvest_total)
    carb_mass_variation(:)=carb_mass_total(:)+carb_mass_variation(:)
    
  !! 17. Write history

    CALL xios_orchidee_send_field("RESOLUTION_X",resolution(:,1))
    CALL xios_orchidee_send_field("RESOLUTION_Y",resolution(:,2))
    CALL xios_orchidee_send_field("CONTFRAC_STOMATE",contfrac(:))
    CALL xios_orchidee_send_field("T2M_MONTH",t2m_month)
    CALL xios_orchidee_send_field("T2M_WEEK",t2m_week)
    CALL xios_orchidee_send_field("TSEASON",Tseason)
    CALL xios_orchidee_send_field("TMIN_SPRING_TIME", Tmin_spring_time)
    CALL xios_orchidee_send_field("ONSET_DATE",onset_date)
    CALL xios_orchidee_send_field("FPC_MAX",fpc_max)
    CALL xios_orchidee_send_field("MAXFPC_LASTYEAR",maxfpc_lastyear)
    CALL xios_orchidee_send_field("HET_RESP",resp_hetero(:,:))
    CALL xios_orchidee_send_field("CO2_FIRE",co2_fire)
    CALL xios_orchidee_send_field("CO2_TAKEN",co2_to_bm)
    CALL xios_orchidee_send_field("LAI",lai)
    CALL xios_orchidee_send_field("VEGET_COV_MAX",veget_cov_max)
    CALL xios_orchidee_send_field("NPP_STOMATE",npp_daily)
    CALL xios_orchidee_send_field("GPP",gpp_daily)
    CALL xios_orchidee_send_field("IND",ind)
    CALL xios_orchidee_send_field("CN_IND",cn_ind)
    CALL xios_orchidee_send_field("WOODMASS_IND",woodmass_ind)
    CALL xios_orchidee_send_field("TOTAL_M",tot_live_biomass)
    CALL xios_orchidee_send_field("MOISTRESS",moiavail_week)
    CALL xios_orchidee_send_field("LEAF_M",biomass(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("SAP_M_AB",biomass(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("SAP_M_BE",biomass(:,:,isapbelow,icarbon))
    CALL xios_orchidee_send_field("HEART_M_AB",biomass(:,:,iheartabove,icarbon))
    CALL xios_orchidee_send_field("HEART_M_BE",biomass(:,:,iheartbelow,icarbon))
    CALL xios_orchidee_send_field("ROOT_M",biomass(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("FRUIT_M",biomass(:,:,ifruit,icarbon))
    CALL xios_orchidee_send_field("RESERVE_M",biomass(:,:,icarbres,icarbon))
    CALL xios_orchidee_send_field("TOTAL_TURN",tot_turnover)
    CALL xios_orchidee_send_field("LEAF_TURN",turnover_daily(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("MAINT_RESP",resp_maint)
    CALL xios_orchidee_send_field("GROWTH_RESP",resp_growth)
    CALL xios_orchidee_send_field("SAP_AB_TURN",turnover_daily(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("ROOT_TURN",turnover_daily(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("FRUIT_TURN",turnover_daily(:,:,ifruit,icarbon))
    CALL xios_orchidee_send_field("TOTAL_BM_LITTER",tot_bm_to_litter(:,:,icarbon))
    CALL xios_orchidee_send_field("LEAF_BM_LITTER",bm_to_litter(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("SAP_AB_BM_LITTER",bm_to_litter(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("SAP_BE_BM_LITTER",bm_to_litter(:,:,isapbelow,icarbon))
    CALL xios_orchidee_send_field("HEART_AB_BM_LITTER",bm_to_litter(:,:,iheartabove,icarbon))
    CALL xios_orchidee_send_field("HEART_BE_BM_LITTER",bm_to_litter(:,:,iheartbelow,icarbon))
    CALL xios_orchidee_send_field("ROOT_BM_LITTER",bm_to_litter(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("FRUIT_BM_LITTER",bm_to_litter(:,:,ifruit,icarbon))
    CALL xios_orchidee_send_field("RESERVE_BM_LITTER",bm_to_litter(:,:,icarbres,icarbon))
    CALL xios_orchidee_send_field("LITTER_STR_AB",litter(:,istructural,:,iabove,icarbon))
    CALL xios_orchidee_send_field("LITTER_MET_AB",litter(:,imetabolic,:,iabove,icarbon))
    CALL xios_orchidee_send_field("LITTER_STR_BE",litter(:,istructural,:,ibelow,icarbon))
    CALL xios_orchidee_send_field("LITTER_MET_BE",litter(:,imetabolic,:,ibelow,icarbon))
    CALL xios_orchidee_send_field("DEADLEAF_COVER",deadleaf_cover)
    CALL xios_orchidee_send_field("TOTAL_SOIL_CARB",tot_litter_soil_carb)
    CALL xios_orchidee_send_field("CARBON_ACTIVE",carbon(:,iactive,:))
    CALL xios_orchidee_send_field("CARBON_SLOW",carbon(:,islow,:))
    CALL xios_orchidee_send_field("CARBON_PASSIVE",carbon(:,ipassive,:))
    CALL xios_orchidee_send_field("LITTERHUM",litterhum_daily)
    CALL xios_orchidee_send_field("TURNOVER_TIME",turnover_time)
    CALL xios_orchidee_send_field("PROD10",prod10)
    CALL xios_orchidee_send_field("FLUX10",flux10)
    CALL xios_orchidee_send_field("PROD100",prod100)
    CALL xios_orchidee_send_field("FLUX100",flux100)
    CALL xios_orchidee_send_field("CONVFLUX",convflux)
    CALL xios_orchidee_send_field("CFLUX_PROD10",cflux_prod10)
    CALL xios_orchidee_send_field("CFLUX_PROD100",cflux_prod100)
    CALL xios_orchidee_send_field("PROD10_HARVEST",prod10_harvest)
    CALL xios_orchidee_send_field("FLUX10_HARVEST",flux10_harvest)
    CALL xios_orchidee_send_field("PROD100_HARVEST",prod100_harvest)
    CALL xios_orchidee_send_field("FLUX100_HARVEST",flux100_harvest)
    CALL xios_orchidee_send_field("CONVFLUX_HARVEST",convflux_harvest)
    CALL xios_orchidee_send_field("CFLUX_PROD10_HARVEST",cflux_prod10_harvest)
    CALL xios_orchidee_send_field("CFLUX_PROD100_HARVEST",cflux_prod100_harvest)
    CALL xios_orchidee_send_field("WOOD_HARVEST",woodharvest/one_year*dt_days)
    CALL xios_orchidee_send_field("WOOD_HARVEST_PFT",woodharvestpft)
    CALL xios_orchidee_send_field("HARVEST_ABOVE",harvest_above)
    CALL xios_orchidee_send_field("VCMAX",vcmax)
    CALL xios_orchidee_send_field("AGE",age)
    CALL xios_orchidee_send_field("HEIGHT",height)
    CALL xios_orchidee_send_field("FIREINDEX",fireindex(:,:))

    ! ipcc history
    ! Carbon stock transformed from gC/m2 into kgC/m2
    CALL xios_orchidee_send_field("cVeg",SUM(tot_live_biomass(:,:,icarbon)*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cVegGrass",sum_cVegGrass/1e3)
    CALL xios_orchidee_send_field("cVegCrop",sum_cVegCrop/1e3)
    CALL xios_orchidee_send_field("cVegTree",sum_cVegTree/1e3)
    CALL xios_orchidee_send_field("cOther",SUM(cOther(:,:,icarbon)*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cLitter",SUM(tot_litter_carb*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cLitterGrass",sum_cLitterGrass/1e3)
    CALL xios_orchidee_send_field("cLitterCrop",sum_cLitterCrop/1e3)
    CALL xios_orchidee_send_field("cLitterTree",sum_cLitterTree/1e3)
    CALL xios_orchidee_send_field("cSoil",SUM(tot_soil_carb*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cSoilGrass",sum_cSoilGrass/1e3)
    CALL xios_orchidee_send_field("cSoilCrop",sum_cSoilCrop/1e3)
    CALL xios_orchidee_send_field("cSoilTree",sum_cSoilTree/1e3)
    CALL xios_orchidee_send_field("cProduct",(prod10_total + prod100_total + prod10_harvest_total + prod100_harvest_total)/1e3)

    cproductlut(:,:)=xios_default_val
    WHERE(fraclut(:,id_psl) > min_sechiba)
       cproductlut(:,id_psl)= ( prod10_total + prod100_total + prod10_harvest_total + prod100_harvest_total)/1e3 &
            / fraclut(:,id_psl)
    ENDWHERE
    CALL xios_orchidee_send_field("cproductlut",cproductlut)

    CALL xios_orchidee_send_field("cMassVariation",carb_mass_variation/1e3/one_day)

    CALL xios_orchidee_send_field("lai_ipcc",SUM(lai*veget_cov_max,dim=2)) ! m2/m2
    
    ! Carbon fluxes transformed from gC/m2/d into kgC/m2/s
    CALL xios_orchidee_send_field("gpp_ipcc",SUM(gpp_daily*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("ra",SUM((resp_maint+resp_growth)*veget_cov_max,dim=2)/1e3/one_day)
    vartmp(:)=zero
    DO j = 2, nvm
       IF ( .NOT. is_tree(j) .AND. natural(j) ) THEN
          vartmp(:) = vartmp(:) + (resp_maint(:,j)+resp_growth(:,j))*veget_cov_max(:,j)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("raGrass",vartmp/1e3/one_day)

    vartmp(:)=zero
    DO j = 2, nvm
       IF (( .NOT. is_tree(j)) .AND. (.NOT. natural(j)) ) THEN
          vartmp(:) = vartmp(:) + (resp_maint(:,j)+resp_growth(:,j))*veget_cov_max(:,j)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("raCrop",vartmp/1e3/one_day)

    vartmp(:)=zero
    DO j = 2, nvm
       IF ( is_tree(j) ) THEN
          vartmp(:) = vartmp(:) + (resp_maint(:,j)+resp_growth(:,j))*veget_cov_max(:,j)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("raTree",vartmp/1e3/one_day)

    CALL xios_orchidee_send_field("npp_ipcc",SUM(npp_daily*veget_cov_max,dim=2)/1e3/one_day)
    vartmp(:)=zero
    DO j = 2, nvm
       IF ( .NOT. is_tree(j) .AND. natural(j) ) THEN
          vartmp(:) = vartmp(:) + npp_daily(:,j)*veget_cov_max(:,j)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("nppGrass",vartmp/1e3/one_day)

    DO j = 2, nvm
       IF ( (.NOT. is_tree(j)) .AND. (.NOT. natural(j)) ) THEN
          vartmp(:) = vartmp(:) + npp_daily(:,j)*veget_cov_max(:,j)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("nppCrop",vartmp/1e3/one_day)

    vartmp(:)=zero
    DO j = 2, nvm
       IF ( is_tree(j) ) THEN
          vartmp(:) = vartmp(:) + npp_daily(:,j)*veget_cov_max(:,j)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("nppTree",vartmp/1e3/one_day)

    CALL xios_orchidee_send_field("rh",SUM(resp_hetero*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("rhLitter",SUM(resp_hetero_litter*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("rhSoil",SUM(resp_hetero_soil*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("HET_RESP_SOIL",resp_hetero_soil(:,:))
    CALL xios_orchidee_send_field("HET_RESP_LITTER",resp_hetero_litter(:,:))

    CALL xios_orchidee_send_field("fFire",SUM(co2_fire*veget_cov_max,dim=2)/1e3/one_day)
    ctotfirelut(:,:)=xios_default_val
    WHERE(fraclut(:,id_psl) > min_sechiba)  
       ctotfirelut(:,id_psl)= SUM(co2_fire*veget_cov_max,dim=2)/1e3/one_day &
            / fraclut(:,id_psl)
    ENDWHERE
    CALL xios_orchidee_send_field("ctotfirelut",ctotfirelut)

    CALL xios_orchidee_send_field("fHarvest",harvest_above/1e3/one_day)
    Call xios_orchidee_send_field("fLuc",cflux_prod_total/1e3/one_day)
    CALL xios_orchidee_send_field("fWoodharvest",cflux_prod_harvest_total/1e3/one_day)
    CALL xios_orchidee_send_field("fDeforestToProduct",SUM((fDeforestToProduct-convfluxpft),dim=2)/1e3/one_day)
    
    flulccproductlut(:,:)=xios_default_val
    flulccproductlut(:,id_psl)=0.
    flulccproductlut(:,id_crp)=0.
    DO j=1,nvm
       IF(natural(j)) THEN
          flulccproductlut(:,id_psl)=flulccproductlut(:,id_psl) + fDeforestToProduct(:,j)
       ELSE
          flulccproductlut(:,id_crp)=flulccproductlut(:,id_crp) + fDeforestToProduct(:,j)
       ENDIF
    ENDDO

    WHERE(fraclut(:,id_psl) > min_sechiba)
       flulccproductlut(:,id_psl)= (flulccproductlut(:,id_psl) &
            +SUM(fHarvestToProduct,dim=2)) &
            /1e3/one_day / fraclut(:,id_psl)
    ELSEWHERE
       flulccproductlut(:,id_psl)= xios_default_val
    ENDWHERE
    WHERE(fraclut(:,id_crp) > min_sechiba)
       flulccproductlut(:,id_crp)= flulccproductlut(:,id_crp) &
            /1e3/one_day / fraclut(:,id_crp)
    ELSEWHERE
       flulccproductlut(:,id_crp)= xios_default_val
    ENDWHERE
    CALL xios_orchidee_send_field("flulccproductlut",flulccproductlut)

    flulccresiduelut(:,:)=xios_default_val
    flulccresiduelut(:,id_psl)=0.
    flulccresiduelut(:,id_crp)=0.
    DO j=1,nvm
       IF(natural(j)) THEN
          flulccresiduelut(:,id_psl) = flulccresiduelut(:,id_psl) + fLulccResidue(:,j)
       ELSE
          flulccresiduelut(:,id_crp) = flulccresiduelut(:,id_crp) + fLulccResidue(:,j)
       ENDIF
    ENDDO

    WHERE(fraclut(:,id_psl) > min_sechiba)
       flulccresiduelut(:,id_psl)= flulccresiduelut(:,id_psl)/1e3/one_day/fraclut(:,id_psl)
    ELSEWHERE
       flulccresiduelut(:,id_psl)= xios_default_val
    ENDWHERE
    WHERE(fraclut(:,id_crp) > min_sechiba)
       flulccresiduelut(:,id_crp)= flulccresiduelut(:,id_crp)/1e3/one_day /fraclut(:,id_crp)
    ELSEWHERE
       flulccresiduelut(:,id_crp)= xios_default_val
    ENDWHERE
    CALL xios_orchidee_send_field("flulccresiduelut",flulccresiduelut)
    CALL xios_orchidee_send_field("flulccresidue",flulccresidue)

    flulccatmlut(:,:)=xios_default_val
    flulccatmlut(:,id_psl)=0.
    flulccatmlut(:,id_crp)=0.
    DO j=1,nvm
       IF(natural(j)) THEN
          flulccatmlut(:,id_psl)=flulccatmlut(:,id_psl)+convfluxpft(:,j)
       ELSE
          flulccatmlut(:,id_crp)=flulccatmlut(:,id_crp)+convfluxpft(:,j)
       ENDIF
    ENDDO

    WHERE(fraclut(:,id_psl) > min_sechiba)
       flulccatmlut(:,id_psl)= (flulccatmlut(:,id_psl)+cflux_prod10+cflux_prod100+cflux_prod_harvest_total) &
            /1e3/one_day / fraclut(:,id_psl)
    ELSEWHERE
       flulccatmlut(:,id_psl)= xios_default_val
    ENDWHERE
    WHERE(fraclut(:,id_crp) > min_sechiba)
       flulccatmlut(:,id_crp)= (flulccatmlut(:,id_crp)+harvest_above)/1e3/one_day / fraclut(:,id_crp)
    ELSEWHERE
       flulccatmlut(:,id_crp)= xios_default_val
    ENDWHERE
    CALL xios_orchidee_send_field("flulccatmlut",flulccatmlut)

   ! co2_to_bm is not added as it is already included in gpp
    CALL xios_orchidee_send_field("nbp",(SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero)-co2_fire) * &
          veget_cov_max,dim=2)-cflux_prod_total-cflux_prod_harvest_total-harvest_above)/1e3/one_day)
    CALL xios_orchidee_send_field("fVegLitter",SUM((tot_bm_to_litter(:,:,icarbon) + tot_turnover(:,:,icarbon))*&
         veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("fLitterSoil",SUM(SUM(soilcarbon_input,dim=2)*veget_cov_max,dim=2)/1e3/one_day)

    ! Carbon stock transformed from gC/m2 into kgC/m2
    CALL xios_orchidee_send_field("cLeaf",SUM(biomass(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cStem",SUM((biomass(:,:,isapabove,icarbon)+biomass(:,:,iheartabove,icarbon))*&
          veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cWood",SUM((biomass(:,:,isapabove,icarbon)+biomass(:,:,iheartabove,icarbon) + &
         biomass(:,:,isapbelow,icarbon) + biomass(:,:,iheartbelow,icarbon))* &
         veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cRoot",SUM(( biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon) + &
         biomass(:,:,iheartbelow,icarbon) )*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cMisc",SUM(( biomass(:,:,icarbres,icarbon) + biomass(:,:,ifruit,icarbon))*&
         veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cLitterAbove",SUM((litter(:,istructural,:,iabove,icarbon)+&
         litter(:,imetabolic,:,iabove,icarbon))*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cLitterBelow",SUM((litter(:,istructural,:,ibelow,icarbon)+&
         litter(:,imetabolic,:,ibelow,icarbon))*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cSoilFast",SUM(carbon(:,iactive,:)*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cSoilMedium",SUM(carbon(:,islow,:)*veget_cov_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("cSoilSlow",SUM(carbon(:,ipassive,:)*veget_cov_max,dim=2)/1e3)

    DO k = 1, ncarb
       csoilpools(:,k) = SUM(carbon(:,k,:)*veget_cov_max,dim=2)/1e3
    END DO
    CALL xios_orchidee_send_field("cSoilPools",csoilpools)

    ! Vegetation fractions [0,100]
    CALL xios_orchidee_send_field("landCoverFrac",veget_cov_max*100)

    ! Carbon fluxes transformed from gC/m2/d into kgC/m2/s
    CALL xios_orchidee_send_field("rGrowth",SUM(resp_growth*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("rMaint",SUM(resp_maint*veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("nppLeaf",SUM(bm_alloc(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3/one_day)
    ! nppStem : from wood above surface
    CALL xios_orchidee_send_field("nppStem",SUM((bm_alloc(:,:,isapabove,icarbon) + bm_alloc(:,:,iheartabove,icarbon)) &
         *veget_cov_max,dim=2)/1e3/one_day)
    ! nppWood : from wood above and below surface 
    CALL xios_orchidee_send_field("nppWood",SUM((bm_alloc(:,:,isapabove,icarbon) + bm_alloc(:,:,iheartabove,icarbon) + &
         bm_alloc(:,:,isapbelow,icarbon) + bm_alloc(:,:,iheartbelow,icarbon)) &
         *veget_cov_max,dim=2)/1e3/one_day)
    ! nppRoot : from wood below surface and fine roots
    CALL xios_orchidee_send_field("nppRoot",SUM((bm_alloc(:,:,isapbelow,icarbon) + bm_alloc(:,:,iheartbelow,icarbon) + &
         bm_alloc(:,:,iroot,icarbon)) * veget_cov_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("nppOther",SUM(( bm_alloc(:,:,ifruit,icarbon) + bm_alloc(:,:,icarbres,icarbon) ) * &
         veget_cov_max,dim=2)/1e3/one_day)

    ! Calculate variables according to LUMIP specifications. 
    clitterlut(:,:) = 0 
    csoillut(:,:) = 0 
    cveglut(:,:) = 0
    lailut(:,:) = 0
    ralut(:,:) = 0
    rhlut(:,:) = 0
    npplut(:,:) = 0

    DO j=1,nvm
       IF (natural(j)) THEN
          clitterlut(:,id_psl) = clitterlut(:,id_psl) + tot_litter_carb(:,j)*veget_cov_max(:,j)/1e3
          csoillut(:,id_psl) = csoillut(:,id_psl) + tot_soil_carb(:,j)*veget_cov_max(:,j)/1e3
          cveglut(:,id_psl) = cveglut(:,id_psl) + tot_live_biomass(:,j,icarbon)*veget_cov_max(:,j)/1e3
          lailut(:,id_psl) = lailut(:,id_psl) + lai(:,j)*veget_cov_max(:,j)
          ralut(:,id_psl) = ralut(:,id_psl) + &
               (resp_maint(:,j)+resp_growth(:,j))*veget_cov_max(:,j)/1e3/one_day
          rhlut(:,id_psl) = rhlut(:,id_psl) + &
               resp_hetero(:,j)*veget_cov_max(:,j)/1e3/one_day
          npplut(:,id_psl) = npplut(:,id_psl) + &
               npp_daily(:,j)*veget_cov_max(:,j)/1e3/one_day
       ELSE
          clitterlut(:,id_crp) = clitterlut(:,id_crp) + tot_litter_carb(:,j)*veget_cov_max(:,j)/1e3
          csoillut(:,id_crp) = csoillut(:,id_crp) + tot_soil_carb(:,j)*veget_cov_max(:,j)/1e3
          cveglut(:,id_crp) = cveglut(:,id_crp) + tot_live_biomass(:,j,icarbon)*veget_cov_max(:,j)/1e3
          lailut(:,id_crp) = lailut(:,id_crp) + lai(:,j)*veget_cov_max(:,j)
          ralut(:,id_crp) = ralut(:,id_crp) + &
               (resp_maint(:,j)+resp_growth(:,j))*veget_cov_max(:,j)/1e3/one_day
          rhlut(:,id_crp) = rhlut(:,id_crp) + &
               resp_hetero(:,j)*veget_cov_max(:,j)/1e3/one_day
          npplut(:,id_crp) = npplut(:,id_crp) + &
               npp_daily(:,j)*veget_cov_max(:,j)/1e3/one_day
       END IF
    END DO

    WHERE (fraclut(:,id_psl)>min_sechiba)
       clitterlut(:,id_psl) = clitterlut(:,id_psl)/fraclut(:,id_psl)
       csoillut(:,id_psl) = csoillut(:,id_psl)/fraclut(:,id_psl)
       cveglut(:,id_psl) = cveglut(:,id_psl)/fraclut(:,id_psl)
       lailut(:,id_psl) = lailut(:,id_psl)/fraclut(:,id_psl)
       ralut(:,id_psl) = ralut(:,id_psl)/fraclut(:,id_psl)
       rhlut(:,id_psl) = rhlut(:,id_psl)/fraclut(:,id_psl)
       npplut(:,id_psl) = npplut(:,id_psl)/fraclut(:,id_psl)
    ELSEWHERE
       clitterlut(:,id_psl) = xios_default_val
       csoillut(:,id_psl) = xios_default_val
       cveglut(:,id_psl) = xios_default_val
       lailut(:,id_psl) = xios_default_val
       ralut(:,id_psl) = xios_default_val
       rhlut(:,id_psl) = xios_default_val
       npplut(:,id_psl) = xios_default_val
    END WHERE

    WHERE (fraclut(:,id_crp)>min_sechiba)
       clitterlut(:,id_crp) = clitterlut(:,id_crp)/fraclut(:,id_crp)
       csoillut(:,id_crp) = csoillut(:,id_crp)/fraclut(:,id_crp)
       cveglut(:,id_crp) = cveglut(:,id_crp)/fraclut(:,id_crp)
       lailut(:,id_crp) = lailut(:,id_crp)/fraclut(:,id_crp)
       ralut(:,id_crp) = ralut(:,id_crp)/fraclut(:,id_crp)
       rhlut(:,id_crp) = rhlut(:,id_crp)/fraclut(:,id_crp)
       npplut(:,id_crp) = npplut(:,id_crp)/fraclut(:,id_crp)
    ELSEWHERE
       clitterlut(:,id_crp) = xios_default_val
       csoillut(:,id_crp) = xios_default_val
       cveglut(:,id_crp) = xios_default_val
       lailut(:,id_crp) = xios_default_val
       ralut(:,id_crp) = xios_default_val
       rhlut(:,id_crp) = xios_default_val
       npplut(:,id_crp) = xios_default_val
    END WHERE

    clitterlut(:,id_pst) = xios_default_val
    clitterlut(:,id_urb) = xios_default_val
    csoillut(:,id_pst)   = xios_default_val
    csoillut(:,id_urb)   = xios_default_val
    cveglut(:,id_pst)    = xios_default_val
    cveglut(:,id_urb)    = xios_default_val
    lailut(:,id_pst)     = xios_default_val
    lailut(:,id_urb)     = xios_default_val
    ralut(:,id_pst)      = xios_default_val
    ralut(:,id_urb)      = xios_default_val
    rhlut(:,id_pst)      = xios_default_val
    rhlut(:,id_urb)      = xios_default_val
    npplut(:,id_pst)     = xios_default_val
    npplut(:,id_urb)     = xios_default_val

    CALL xios_orchidee_send_field("clitterlut",clitterlut)
    CALL xios_orchidee_send_field("csoillut",csoillut)
    CALL xios_orchidee_send_field("cveglut",cveglut)
    CALL xios_orchidee_send_field("lailut",lailut)
    CALL xios_orchidee_send_field("ralut",ralut)
    CALL xios_orchidee_send_field("rhlut",rhlut)
    CALL xios_orchidee_send_field("npplut",npplut)

    CALL histwrite_p (hist_id_stomate, 'RESOLUTION_X', itime, &
         resolution(:,1), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'RESOLUTION_Y', itime, &
         resolution(:,2), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CONTFRAC', itime, &
         contfrac(:), npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AB', itime, &
         litter(:,istructural,:,iabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AB', itime, &
         litter(:,imetabolic,:,iabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_BE', itime, &
         litter(:,istructural,:,ibelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_BE', itime, &
         litter(:,imetabolic,:,ibelow,icarbon), npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'DEADLEAF_COVER', itime, &
         deadleaf_cover, npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'TOTAL_SOIL_CARB', itime, &
         tot_litter_soil_carb, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_ACTIVE', itime, &
         carbon(:,iactive,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_SLOW', itime, &
         carbon(:,islow,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_PASSIVE', itime, &
         carbon(:,ipassive,:), npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'T2M_MONTH', itime, &
         t2m_month, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'T2M_WEEK', itime, &
         t2m_week, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'TSEASON', itime, &
         Tseason, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'TMIN_SPRING_TIME', itime, &
         Tmin_spring_time, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ONSET_DATE', itime, &
         onset_date(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FPC_MAX', itime, &
         fpc_max, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MAXFPC_LASTYEAR', itime, &
         maxfpc_lastyear, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HET_RESP', itime, &
         resp_hetero(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FIREINDEX', itime, &
         fireindex(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTERHUM', itime, &
         litterhum_daily, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CO2_FIRE', itime, &
         co2_fire, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CO2_TAKEN', itime, &
         co2_to_bm, npts*nvm, horipft_index)
    ! land cover change
    CALL histwrite_p (hist_id_stomate, 'CONVFLUX', itime, &
         convflux, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD10', itime, &
         cflux_prod10, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD100', itime, &
         cflux_prod100, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CONVFLUX_HARVEST', itime, &
         convflux_harvest, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD10_HARVEST', itime, &
         cflux_prod10_harvest, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD100_HARVES', itime, &
         cflux_prod100_harvest, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'WOOD_HARVEST', itime, &
         woodharvest/one_year*dt_days, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'WOOD_HARVEST_PFT', itime, &
         woodharvestpft, npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'HARVEST_ABOVE', itime, &
         harvest_above, npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'LAI', itime, &
         lai, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'VEGET_COV_MAX', itime, &
         veget_cov_max, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NPP', itime, &
         npp_daily, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'GPP', itime, &
         gpp_daily, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'IND', itime, &
         ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CN_IND', itime, &
         cn_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'WOODMASS_IND', itime, &
         woodmass_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TOTAL_M', itime, &
         tot_live_biomass(:,:,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_M', itime, &
         biomass(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_M_AB', itime, &
         biomass(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_M_BE', itime, &
         biomass(:,:,isapbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_M_AB', itime, &
         biomass(:,:,iheartabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_M_BE', itime, &
         biomass(:,:,iheartbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_M', itime, &
         biomass(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FRUIT_M', itime, &
         biomass(:,:,ifruit,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'RESERVE_M', itime, &
         biomass(:,:,icarbres,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TOTAL_TURN', itime, &
         tot_turnover(:,:,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_TURN', itime, &
         turnover_daily(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_AB_TURN', itime, &
         turnover_daily(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_TURN', itime, &
         turnover_daily(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FRUIT_TURN', itime, &
         turnover_daily(:,:,ifruit,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TOTAL_BM_LITTER', itime, &
         tot_bm_to_litter(:,:,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_BM_LITTER', itime, &
         bm_to_litter(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_AB_BM_LITTER', itime, &
         bm_to_litter(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_BE_BM_LITTER', itime, &
         bm_to_litter(:,:,isapbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_AB_BM_LITTER', itime, &
         bm_to_litter(:,:,iheartabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_BE_BM_LITTER', itime, &
         bm_to_litter(:,:,iheartbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_BM_LITTER', itime, &
         bm_to_litter(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FRUIT_BM_LITTER', itime, &
         bm_to_litter(:,:,ifruit,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'RESERVE_BM_LITTER', itime, &
         bm_to_litter(:,:,icarbres,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MAINT_RESP', itime, &
         resp_maint, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'GROWTH_RESP', itime, &
         resp_growth, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'AGE', itime, &
         age, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEIGHT', itime, &
         height, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MOISTRESS', itime, &
         moiavail_week, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'VCMAX', itime, &
         vcmax, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TURNOVER_TIME', itime, &
         turnover_time, npts*nvm, horipft_index)
    ! land cover change
    CALL histwrite_p (hist_id_stomate, 'PROD10', itime, &
         prod10, npts*11, horip11_index)
    CALL histwrite_p (hist_id_stomate, 'PROD100', itime, &
         prod100, npts*101, horip101_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX10', itime, &
         flux10, npts*10, horip10_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX100', itime, &
         flux100, npts*100, horip100_index)
    CALL histwrite_p (hist_id_stomate, 'PROD10_HARVEST', itime, &
         prod10_harvest, npts*11, horip11_index)
    CALL histwrite_p (hist_id_stomate, 'PROD100_HARVEST', itime, &
         prod100_harvest, npts*101, horip101_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX10_HARVEST', itime, &
         flux10_harvest, npts*10, horip10_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX100_HARVEST', itime, &
         flux100_harvest, npts*100, horip100_index)

    IF ( hist_id_stomate_IPCC > 0 ) THEN
       vartmp(:)=SUM(tot_live_biomass(:,:,icarbon)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cVeg", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_litter_carb*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_soil_carb*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=(prod10_total + prod100_total + prod10_harvest_total + prod100_harvest_total)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cProduct", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=carb_mass_variation/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "cMassVariation", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(lai*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate_IPCC, "lai", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(gpp_daily*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "gpp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((resp_maint+resp_growth)*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "ra", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(npp_daily*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "npp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_hetero*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "rh", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(co2_fire*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fFire", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=harvest_above/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fHarvest", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=cflux_prod_total/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fLuc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=cflux_prod_harvest_total/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fWoodharvest", itime, &
            vartmp, npts, hori_index)
       ! co2_to_bm is not added as it is already included in gpp
       vartmp(:)=(SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero)-co2_fire) &
            &        *veget_cov_max,dim=2)-cflux_prod_total-cflux_prod_harvest_total-harvest_above)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "nbp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((tot_bm_to_litter(:,:,icarbon) + tot_turnover(:,:,icarbon))*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fVegLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(SUM(soilcarbon_input,dim=2)*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "fLitterSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(biomass(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLeaf", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((biomass(:,:,isapabove,icarbon)+biomass(:,:,iheartabove,icarbon))*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cStem", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon) + biomass(:,:,iheartbelow,icarbon) ) &
            &        *veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cRoot", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,icarbres,icarbon) + biomass(:,:,ifruit,icarbon))*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cMisc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,iabove,icarbon)+litter(:,imetabolic,:,iabove,icarbon))*&
            veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitterAbove", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,ibelow,icarbon)+litter(:,imetabolic,:,ibelow,icarbon))*&
            veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitterBelow", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,iactive,:)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilFast", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,islow,:)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilMedium", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,ipassive,:)*veget_cov_max,dim=2)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilSlow", itime, &
            vartmp, npts, hori_index)
       DO j=1,nvm
          histvar(:,j)=veget_cov_max(:,j)*100
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "landCoverFrac", itime, &
            histvar, npts*nvm, horipft_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_deciduous(j)) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "treeFracPrimDec", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_evergreen(j)) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "treeFracPrimEver", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( .NOT.(is_c4(j)) ) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "c3PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( is_c4(j) ) THEN
             vartmp(:) = vartmp(:) + veget_cov_max(:,j)*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "c4PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=SUM(resp_growth*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "rGrowth", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_maint*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "rMaint", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(bm_alloc(:,:,ileaf,icarbon)*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "nppLeaf", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((bm_alloc(:,:,isapabove,icarbon) + bm_alloc(:,:,iheartabove,icarbon))*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "nppStem", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( bm_alloc(:,:,isapbelow,icarbon) + bm_alloc(:,:,iroot,icarbon) )*veget_cov_max,dim=2)/1e3/one_day
       CALL histwrite_p (hist_id_stomate_IPCC, "nppRoot", itime, &
            vartmp, npts, hori_index)

       CALL histwrite_p (hist_id_stomate_IPCC, 'RESOLUTION_X', itime, &
            resolution(:,1), npts, hori_index)
       CALL histwrite_p (hist_id_stomate_IPCC, 'RESOLUTION_Y', itime, &
            resolution(:,2), npts, hori_index)
       CALL histwrite_p (hist_id_stomate_IPCC, 'CONTFRAC', itime, &
            contfrac(:), npts, hori_index)

    ENDIF

    IF (printlev>=4) WRITE(numout,*) 'Leaving stomate_lpj'

  END SUBROUTINE StomateLpj


!! ================================================================================================================================
!! SUBROUTINE   : harvest
!!
!>\BRIEF        Harvest of croplands
!!
!! DESCRIPTION  : To take into account biomass harvest from crop (mainly to take 
!! into account for the reduced litter input and then decreased soil carbon. it is a 
!! constant (40\%) fraction of above ground biomass.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::harvest_above the harvested biomass
!!
!! REFERENCE(S) :
!! - Piao, S., P. Ciais, P. Friedlingstein, N. de Noblet-Ducoudre, P. Cadule, N. Viovy, and T. Wang. 2009. 
!!   Spatiotemporal patterns of terrestrial carbon cycle during the 20th century. Global Biogeochemical 
!!   Cycles 23:doi:10.1029/2008GB003339.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE harvest(npts, dt_days, veget_cov_max, &
       bm_to_litter, turnover_daily, &
       harvest_above)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                    :: npts            !! Domain size (unitless) 
    REAL(r_std), INTENT(in)                                :: dt_days         !! Time step (days)                               
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: veget_cov_max       !! new "maximal" coverage fraction of a PFT (LAI -> 
                                                                              !! infinity) on ground @tex $(m^2 m^{-2})$ @endtex 
    
   !! 0.2 Output variables
   
   !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! [DISPENSABLE] conversion of biomass to litter 
                                                                                     !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily   !! Turnover rates 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: harvest_above    !! harvest above ground biomass for agriculture 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    !! 0.4 Local variables

    INTEGER(i_std)                                         :: i, j, k, l, m    !! indices                       
    REAL(r_std)                                            :: above_old        !! biomass of previous time step 
                                                                               !! @tex $(gC m^{-2})$ @endtex 
!_ ================================================================================================================================

  !! 1. Yearly initialisation

    above_old             = zero
    harvest_above         = zero

    DO i = 1, npts
       DO j = 1,nvm
          IF (.NOT. natural(j)) THEN
             above_old = turnover_daily(i,j,ileaf,icarbon) + turnover_daily(i,j,isapabove,icarbon) + &
                  &       turnover_daily(i,j,iheartabove,icarbon) + turnover_daily(i,j,ifruit,icarbon) + &
                  &       turnover_daily(i,j,icarbres,icarbon) + turnover_daily(i,j,isapbelow,icarbon) + &
                  &       turnover_daily(i,j,iheartbelow,icarbon) + turnover_daily(i,j,iroot,icarbon)

             turnover_daily(i,j,ileaf,icarbon) = turnover_daily(i,j,ileaf,icarbon)*frac_turnover_daily
             turnover_daily(i,j,isapabove,icarbon) = turnover_daily(i,j,isapabove,icarbon)*frac_turnover_daily
             turnover_daily(i,j,isapbelow,icarbon) = turnover_daily(i,j,isapbelow,icarbon)*frac_turnover_daily
             turnover_daily(i,j,iheartabove,icarbon) = turnover_daily(i,j,iheartabove,icarbon)*frac_turnover_daily
             turnover_daily(i,j,iheartbelow,icarbon) = turnover_daily(i,j,iheartbelow,icarbon)*frac_turnover_daily
             turnover_daily(i,j,iroot,icarbon) = turnover_daily(i,j,iroot,icarbon)*frac_turnover_daily
             turnover_daily(i,j,ifruit,icarbon) = turnover_daily(i,j,ifruit,icarbon)*frac_turnover_daily
             turnover_daily(i,j,icarbres,icarbon) = turnover_daily(i,j,icarbres,icarbon)*frac_turnover_daily
             harvest_above(i)  = harvest_above(i) + veget_cov_max(i,j) * above_old *(un - frac_turnover_daily)
          ENDIF
       ENDDO
    ENDDO

!!$    harvest_above = harvest_above
  END SUBROUTINE harvest
END MODULE stomate_lpj
