subroutine Each_Year
!**** Processes that are required each year, prior to any process-based simulation steps.
!****
!**** R. Boone   Last modified: September 13, 2014
  use Parameter_Vars
  use Structures
  implicit none

  real shrub_c_sum, tree_c_sum
  integer icell, iunit, ipart, ifacet, in_year, unit_id, unit_cnt

  do icell=1,range_cells                 ! Process all of the cells classed as rangeland only
    iunit = Rng(icell)%range_type

    Rng(icell)%annual_evapotranspiration = 0.0
    Rng(icell)%neg_error_count = 0
    Rng(icell)%large_error_count = 0

    ! Fill the tree carbon allocation array.  This is dynamic in Century 4.5, static in 4.0, and static here.  (Calculated each year, but no matter)
    ! This will be approximate, as all five pieces are not specified (but could be).
    shrub_c_sum = 0.0
    tree_c_sum = 0.0
    do ipart = 1, WOODY_PARTS
      shrub_c_sum = shrub_c_sum + Parms(iunit)%shrub_carbon(ipart, ALIVE) + Parms(iunit)%shrub_carbon(ipart, DEAD)
      tree_c_sum = tree_c_sum + Parms(iunit)%tree_carbon(ipart, ALIVE) + Parms(iunit)%tree_carbon(ipart, DEAD)
    end do
    if (shrub_c_sum .gt. 0.0) then
      do ipart = 1, WOODY_PARTS
        Rng(icell)%carbon_allocation(S_FACET, ipart) = ( Parms(iunit)%shrub_carbon(ipart, ALIVE) + &
                                                         Parms(iunit)%shrub_carbon(ipart, DEAD) )/ shrub_c_sum
      end do
    else
      do ipart = 1, WOODY_PARTS
        Rng(icell)%carbon_allocation(S_FACET, ipart) = 0.0
      end do
    end if
    if (tree_c_sum .gt. 0.0) then
      do ipart = 1, WOODY_PARTS
        Rng(icell)%carbon_allocation(T_FACET, ipart) = ( Parms(iunit)%tree_carbon(ipart, ALIVE) + &
                                                         Parms(iunit)%shrub_carbon(ipart, DEAD) ) / tree_c_sum
      end do
    else
      do ipart = 1, WOODY_PARTS
        Rng(icell)%carbon_allocation(T_FACET, ipart) = 0.0
      end do
    end if

    ! Clearing-out the shrub_carbon and tree_carbon variables, such that they represent the contribution of new carbon
    ! to woody parts for the year in question.
    do ipart = 1, WOODY_PARTS
      Rng(icell)%shrub_carbon(ipart) = 0.0
      Rng(icell)%tree_carbon(ipart) = 0.0
      Rng(icell)%shrub_nitrogen(ipart) = 0.0
      Rng(icell)%tree_nitrogen(ipart) = 0.0
    end do

    do ifacet = 1, FACETS
      ! Some of the dead plant materials are long-term storage (e.g., dead coarse roots, standing dead)
      ! Others don't exist in Century, but are incorporated here for completeness.  They are storage holders for use in annual
      ! model tracking.  They need to be cleared-out prior to each year's simulation.
      ! (Note: Arrays can be zeroed out with a single call, but not as clear)
      ! Totals summed across facets are to be reset each year
      Rng(icell)%dead_total_fine_branch_carbon = 0.0
      Rng(icell)%dead_total_fine_branch_nitrogen = 0.0
      Rng(icell)%dead_total_coarse_branch_carbon = 0.0
      Rng(icell)%dead_total_coarse_branch_nitrogen = 0.0
      Rng(icell)%dead_total_coarse_root_carbon = 0.0
      Rng(icell)%dead_total_coarse_root_nitrogen = 0.0
      ! Other dead placeholders may be zeroed-out as well.  The only reason they would store anything after the 
      ! call to DECOMPOSITION would be from rounding errors.  And DEAD_LEAF_CARBON and NITROGEN plus FINE_BRANCH_CARBON and NITROGEN
      ! aren't used in modeling at all, they are only accumulators.   The death of those elements go to STANDING_DEAD_CARBON and NITROGEN
      ! where they are partitioned.  STANDING_DEAD* SHOULD NOT be zeroed out, as that material can accumulate over more than one season.
      

      ! Storage areas that serve as accumulators only should be reset to 0 each year, to avoid a continual accumulation of
      ! values through an entire simulation.  There may be a case or two where such accumulations are helpful, but most
      ! will be reset to zero.   Some of these don't include facets and will be reset more than required, but no matter here.
      Rng(icell)%evaporation = 0.0
      Rng(icell)%maintain_respiration(ifacet) = 0.0
      Rng(icell)%respiration_annual = 0.0
      Rng(icell)%nitrogen_fixed(ifacet) = 0.0

      ! Recalculate the proportion of residue that is lignin, which follows from annual precipitation (CMPLIG.F in Century.  No equilvalent in Savanna)
      Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX) = Rng(icell)%lignin_leaf(ifacet) + &
                              Parms(iunit)%lignin_content_fraction_and_precip(1,SURFACE_INDEX) + &
                            ( Parms(iunit)%lignin_content_fraction_and_precip(2,SURFACE_INDEX) * &
                              Globe(Rng(icell)%x, Rng(icell)%y)%precip_average ) / 2.0
      Rng(icell)%plant_lignin_fraction(ifacet, SOIL_INDEX) = Rng(icell)%lignin_fine_root(ifacet) + &
                              Parms(iunit)%lignin_content_fraction_and_precip(1,SOIL_INDEX) + &
                            ( Parms(iunit)%lignin_content_fraction_and_precip(2,SOIL_INDEX) * &
                              Globe(Rng(icell)%x, Rng(icell)%y)%precip_average ) / 2.0
      Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX) = max(0.02, Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX) )
      Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX) = min(0.50, Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX) )
      Rng(icell)%plant_lignin_fraction(ifacet, SOIL_INDEX) = max(0.02, Rng(icell)%plant_lignin_fraction(ifacet, SOIL_INDEX) )
      Rng(icell)%plant_lignin_fraction(ifacet, SOIL_INDEX) = min(0.50, Rng(icell)%plant_lignin_fraction(ifacet, SOIL_INDEX) )

      ! Fire modeling
      Rng(icell)%burned_carbon = 0.0
      Rng(icell)%burned_nitrogen = 0.0
      Rng(icell)%fire_severity = 0.0
   end do                                          ! Facet loop

  end do                                           ! Cell loop

  ! Opening the CO2 effects file each month, just for simplicity
  open(SHORT_USE_FILE, FILE=parm_path(1:len_trim(parm_path))//Sim_Parm%co2effect_file_name, ACTION='READ', IOSTAT=ioerr)
  if (ioerr == 0) then
    in_year = -9999
    read(SHORT_USE_FILE,*) unit_cnt            ! GRange never knows the number of landscape units, so required at the top of file or some other pathway
    read(SHORT_USE_FILE,*)                     ! Skip the header information
    do while (in_year .ne. year)
      read(SHORT_USE_FILE,*) in_year, (Parms(unit_id)%effect_of_co2_on_production(H_FACET), &
       Parms(unit_id)%effect_of_co2_on_production(S_FACET), Parms(unit_id)%effect_of_co2_on_production(T_FACET), &
       unit_id=1,unit_cnt)
    end do
    do icell=1,range_cells                 ! Process all of the cells classed as rangeland only
      iunit = Rng(icell)%range_type
      Rng(icell)%co2_effect_on_production(H_FACET) = Parms(iunit)%effect_of_co2_on_production(H_FACET) 
      Rng(icell)%co2_effect_on_production(S_FACET) = Parms(iunit)%effect_of_co2_on_production(S_FACET) 
      Rng(icell)%co2_effect_on_production(T_FACET) = Parms(iunit)%effect_of_co2_on_production(T_FACET) 
    end do
    close(SHORT_USE_FILE)
  else
    write(*,*) 'There is a problem updating the CO2 effect on production values'
    stop
  end if

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'EACH_YR')

end subroutine


subroutine Update_Vegetation (icell)
!**** Update metrics that summarize vegetation.  Adapted from Century or from scratch
!****
!**** R. Boone   Last modified: July 18, 2014
  use Parameter_Vars
  use Structures
  implicit none

  real wdbmas, basf, data_val(10), linear
  real avg_a_live_biomass, avg_b_live_biomass, total_cover, frac_cover
  real biomass_live_per_layer(V_LYRS)
  integer icell, iunit, ifacet, i, ilyr

  iunit = Rng(icell)%range_type

  ! Update tree basal area
  wdbmas = (Rng(icell)%fine_branch_carbon(T_FACET) + Rng(icell)%coarse_branch_carbon(T_FACET)) * 2.0
  if (wdbmas .le. 0.0) then
    wdbmas = 50.                                                                                       ! Adjusted to allow trees to grow from zero and avoid underflows.
  end if
  basf = (wdbmas / (0.88 * ((wdbmas * 0.01)**0.635)))                                                  ! Division by zero avoided above.
  if (basf .lt. 250.0) then
    basf = basf * Parms(iunit)%tree_basal_area_to_wood_biomass
  end if
  Rng(icell)%tree_basal_area = wdbmas / basf                                                           ! Setting WDBMAS and the use of a parameter should avoid division by 0.

  Rng(icell)%water_function = 1. / (1. + 4. * exp(-6. * Rng(icell)%relative_water_content(1)))         ! 1+ etc. should avoid division by zero error.
  Rng(icell)%total_litter_carbon(SURFACE_INDEX) = Rng(icell)%litter_structural_carbon(SURFACE_INDEX) + &
                                                 Rng(icell)%litter_metabolic_carbon(SURFACE_INDEX)
  if (Rng(icell)%total_litter_carbon(SURFACE_INDEX) .lt. 0.0) then
    Rng(icell)%total_litter_carbon(SURFACE_INDEX) = 0.0
  end if
  Rng(icell)%total_litter_carbon(SOIL_INDEX) = Rng(icell)%litter_structural_carbon(SOIL_INDEX) + &
                                               Rng(icell)%litter_metabolic_carbon(SOIL_INDEX)
  if (Rng(icell)%total_litter_carbon(SOIL_INDEX) .lt. 0.0) then
    Rng(icell)%total_litter_carbon(SOIL_INDEX) = 0.0
  end if
  Rng(icell)%total_litter_nitrogen(SURFACE_INDEX) = Rng(icell)%litter_structural_nitrogen(SURFACE_INDEX) + &
                                                 Rng(icell)%litter_metabolic_nitrogen(SURFACE_INDEX)
  if (Rng(icell)%total_litter_nitrogen(SURFACE_INDEX) .lt. 0.0) then
    Rng(icell)%total_litter_nitrogen(SURFACE_INDEX) = 0.0
  end if
  Rng(icell)%total_litter_nitrogen(SOIL_INDEX) = Rng(icell)%litter_structural_nitrogen(SOIL_INDEX) + &
                                               Rng(icell)%litter_metabolic_nitrogen(SOIL_INDEX)
  if (Rng(icell)%total_litter_nitrogen(SOIL_INDEX) .lt. 0.0) then
    Rng(icell)%total_litter_nitrogen(SOIL_INDEX) = 0.0
  end if

  ! Update the phenology of the plants
  do ifacet = 1, FACETS
    do i = 1,10
      data_val(i) = Parms(iunit)%degree_days_phen(ifacet,i)
    end do
    Rng(icell)%phenology(ifacet) = linear(Rng(icell)%heat_accumulation, data_val, 5)
    if (Rng(icell)%heat_accumulation .ge. Parms(iunit)%degree_days_reset(ifacet)) then
      Rng(icell)%phenology(ifacet) = 0.0
    end if
  end do
  
  ! The following is likely not very good for herbs, using the same values as trees, but still a helpful index.  Could be expanded to include LAI to biomass relationship specific to herbs
  ! The following only considers the tallest facet.  Could (should?) include shrubs and herbs in tree facet, etc?   Separated out from a facet loop here incase that change is made.
  do ifacet = 1, FACETS 
    Rng(icell)%leaf_area_index(ifacet) = Rng(icell)%leaf_carbon(ifacet) * (2.5 * Parms(iunit)%biomass_to_leaf_area_index_factor)
  end do

  Rng(icell)%soil_total_carbon = Rng(icell)%fast_soil_carbon(SOIL_INDEX) + Rng(icell)%intermediate_soil_carbon + &
                                 Rng(icell)%passive_soil_carbon + Rng(icell)%litter_structural_carbon(SOIL_INDEX) + &
                                 Rng(icell)%litter_metabolic_carbon(SOIL_INDEX)
  Rng(icell)%carbon_nitrogen_ratio = ( Rng(icell)%fast_soil_carbon(SOIL_INDEX) + Rng(icell)%intermediate_soil_carbon + &
                                 Rng(icell)%passive_soil_carbon + Rng(icell)%litter_structural_carbon(SOIL_INDEX) + &
                                 Rng(icell)%litter_metabolic_carbon(SOIL_INDEX) ) / &
                                     ( Rng(icell)%fast_soil_nitrogen(SOIL_INDEX) + Rng(icell)%intermediate_soil_nitrogen + &
                                 Rng(icell)%passive_soil_nitrogen + Rng(icell)%litter_structural_nitrogen(SOIL_INDEX) + &
                                 Rng(icell)%litter_metabolic_nitrogen(SOIL_INDEX) ) 

  ! Moving the following from WATER_LOSS, since it was only being done with no snow present.
  avg_a_live_biomass = 0.0
  avg_b_live_biomass = 0.0
 
  ! ABOVEGROUND
  ! Using method used in productivity.  Does not use plant populations in layers, but uses the facets instead.  Not at precise but less prone to vast swings. 
  total_cover =   Rng(icell)%facet_cover(H_FACET) + Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(T_FACET)
  if (total_cover .gt. 0.000001) then
    ! HERBS
    frac_cover   = Rng(icell)%facet_cover(H_FACET) / total_cover
    biomass_live_per_layer(H_LYR) = ( Rng(icell)%leaf_carbon(H_FACET) + Rng(icell)%seed_carbon(H_FACET) ) * 2.5 * frac_cover
    frac_cover   = Rng(icell)%facet_cover(S_FACET) / total_cover
    biomass_live_per_layer(H_S_LYR) = ( Rng(icell)%leaf_carbon(H_FACET) + Rng(icell)%seed_carbon(H_FACET) ) * 2.5 * frac_cover
    frac_cover   = Rng(icell)%facet_cover(T_FACET) / total_cover
    biomass_live_per_layer(H_T_LYR) = ( Rng(icell)%leaf_carbon(H_FACET) + Rng(icell)%seed_carbon(H_FACET) ) * 2.5 * frac_cover
    ! SHRUBS
    if ( ( Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(T_FACET) ) .gt. 0.000001 ) then
      frac_cover   = Rng(icell)%facet_cover(S_FACET) / ( Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(T_FACET) )
      biomass_live_per_layer(S_LYR) =( Rng(icell)%leaf_carbon(S_FACET) + Rng(icell)%seed_carbon(S_FACET) + &
                          Rng(icell)%fine_branch_carbon(S_FACET) + Rng(icell)%coarse_branch_carbon(S_FACET) ) * 2.5 * frac_cover
      frac_cover   = Rng(icell)%facet_cover(T_FACET) / ( Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(T_FACET) )
      biomass_live_per_layer(S_T_LYR) = ( Rng(icell)%leaf_carbon(S_FACET) + Rng(icell)%seed_carbon(S_FACET) + &
                          Rng(icell)%fine_branch_carbon(S_FACET) + Rng(icell)%coarse_branch_carbon(S_FACET) ) * 2.5 * frac_cover
    else
      biomass_live_per_layer(S_LYR) = 0.0
      biomass_live_per_layer(T_LYR) = 0.0
    end if
    ! TREES
    frac_cover = Rng(icell)%facet_cover(T_FACET)
    biomass_live_per_layer(T_LYR) = ( Rng(icell)%leaf_carbon(T_FACET) + Rng(icell)%seed_carbon(T_FACET) + &
                          Rng(icell)%fine_branch_carbon(T_FACET) + Rng(icell)%coarse_branch_carbon(T_FACET) ) * 2.5 * frac_cover
   
    do ilyr = 1, V_LYRS
      avg_a_live_biomass = avg_a_live_biomass + biomass_live_per_layer(ilyr)
    end do
  else
    avg_a_live_biomass = 0.0                 ! There is no cover on the cell
  end if

  ! BELOWGROUND
  if (total_cover .gt. 0.000001) then
    ! HERBS
    frac_cover   = Rng(icell)%facet_cover(H_FACET) / total_cover
    biomass_live_per_layer(H_LYR) = Rng(icell)%fine_root_carbon(H_FACET) * 2.5 * frac_cover
    frac_cover   = Rng(icell)%facet_cover(S_FACET) / total_cover
    biomass_live_per_layer(H_S_LYR) = Rng(icell)%fine_root_carbon(H_FACET) * 2.5 * frac_cover
    frac_cover   = Rng(icell)%facet_cover(T_FACET) / total_cover
    biomass_live_per_layer(H_T_LYR) = Rng(icell)%fine_root_carbon(H_FACET) * 2.5 * frac_cover
    ! SHRUBS
    if ( ( Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(T_FACET) ) .gt. 0.000001 ) then
      frac_cover   = Rng(icell)%facet_cover(S_FACET) / ( Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(T_FACET) )
      biomass_live_per_layer(S_LYR) = &
                         ( Rng(icell)%fine_root_carbon(S_FACET) + Rng(icell)%coarse_root_carbon(S_FACET) ) * 2.5 * frac_cover
      frac_cover   = Rng(icell)%facet_cover(T_FACET) / ( Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(T_FACET) )
      biomass_live_per_layer(S_T_LYR) = &
                         ( Rng(icell)%fine_root_carbon(S_FACET) + Rng(icell)%coarse_root_carbon(S_FACET) ) * 2.5 * frac_cover
    else
      biomass_live_per_layer(S_LYR) = 0.0
      biomass_live_per_layer(T_LYR) = 0.0
    end if
    ! TREES
    frac_cover = Rng(icell)%facet_cover(T_FACET)
    biomass_live_per_layer(T_LYR) = &
                         ( Rng(icell)%fine_root_carbon(T_FACET) + Rng(icell)%coarse_root_carbon(T_FACET) ) * 2.5 * frac_cover
   
    do ilyr = 1, V_LYRS
      avg_b_live_biomass = avg_b_live_biomass + biomass_live_per_layer(ilyr)
    end do
  else
    avg_b_live_biomass = 0.0                 ! There is no cover on the cell
  end if

  ! The following is used in Decomposition, and perhaps elsewhere.
  Rng(icell)%total_aground_live_biomass = avg_a_live_biomass
  Rng(icell)%total_bground_live_biomass = avg_b_live_biomass

  ! Calculate monthly net primary productivity
  Rng(icell)%monthly_net_primary_production = &
     ( Rng(icell)%total_pot_prod_limited_by_n(H_LYR) * Rng(icell)%facet_cover(H_FACET) ) + &
     ( Rng(icell)%total_pot_prod_limited_by_n(H_S_LYR) * Rng(icell)%facet_cover(S_FACET) ) + &
     ( Rng(icell)%total_pot_prod_limited_by_n(H_T_LYR) * Rng(icell)%facet_cover(T_FACET) ) + &
     ( Rng(icell)%total_pot_prod_limited_by_n(S_LYR) * Rng(icell)%facet_cover(S_FACET) ) + &
     ( Rng(icell)%total_pot_prod_limited_by_n(S_T_LYR) * Rng(icell)%facet_cover(T_FACET) ) + &
     ( Rng(icell)%total_pot_prod_limited_by_n(T_LYR) * Rng(icell)%facet_cover(T_FACET) )

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'UPD_VEG')

end subroutine


subroutine Wrap_Up
!**** Tasks to end a simulation.
!****
!**** R. Boone   Last modified: April 24, 2011
  use Parameter_Vars
  use Structures
  implicit none

  integer yr_cnt, Vals(2)
  real end_time

  !***** Handle State Variable Request
  if (Sim_Parm%state_var_flag .eq. 1 .or. Sim_Parm%state_var_flag .eq. 3) then
    write(*,*) ' '
    write(*,*) 'Writing the state of the model to: ', app_path(1:len_trim(app_path))//Sim_Parm%state_var_file_out
    call Write_State
  end if

  open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//EXCEED_NAME, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)
  write(SHORT_USE_FILE,*) 'Number_of_times_each_variable_cell-month_combination_exceeded_a_threshold_and_was_reset Very_Large Zero'
  if (ioerr == 0) then
    ! The following list could be placed in a file, but hardwiring the list seems more secure, that is hidden from users
    Vals(1) = Exceed%last_month_day_length(1)  ;  Vals(2) = Exceed%last_month_day_length(2)
    call One_Out ('last_month_day_length', 21, Vals)
    Vals(1) = Exceed%day_length(1)  ;  Vals(2) = Exceed%day_length(2)
    call One_Out ('day_length', 10, Vals)
    Vals(1) = Exceed%heat_accumulation(1)  ;  Vals(2) = Exceed%heat_accumulation(2)
    call One_Out ('heat_accumulation', 17, Vals)
    Vals(1) = Exceed%facet_cover(1)  ;  Vals(2) = Exceed%facet_cover(2)
    call One_Out ('facet_cover', 11, Vals)
    Vals(1) = Exceed%total_population(1)  ;  Vals(2) = Exceed%total_population(2)
    call One_Out ('total_population', 16, Vals)
    Vals(1) = Exceed%bare_cover(1)  ;  Vals(2) = Exceed%bare_cover(2)
    call One_Out ('bare_cover', 10, Vals)
    Vals(1) = Exceed%prop_annual_decid(1)  ;  Vals(2) = Exceed%prop_annual_decid(2)
    call One_Out ('prop_annual_decid', 17, Vals)
    Vals(1) = Exceed%pot_evap(1)  ;  Vals(2) = Exceed%pot_evap(2)
    call One_Out ('pot_evap', 8, Vals)
    Vals(1) = Exceed%evaporation(1)  ;  Vals(2) = Exceed%evaporation(2)
    call One_Out ('evaporation', 11, Vals)
    Vals(1) = Exceed%snow(1)  ;  Vals(2) = Exceed%snow(2)
    call One_Out ('snow', 4, Vals)
    Vals(1) = Exceed%snow_liquid(1)  ;  Vals(2) = Exceed%snow_liquid(2)
    call One_Out ('snow_liquid', 11, Vals)
    Vals(1) = Exceed%melt(1)  ;  Vals(2) = Exceed%melt(2)
    call One_Out ('melt', 4, Vals)
    Vals(1) = Exceed%pet_remaining(1)  ;  Vals(2) = Exceed%pet_remaining(2)
    call One_Out ('pet_remaining', 13, Vals)
    Vals(1) = Exceed%ppt_soil(1)  ;  Vals(2) = Exceed%ppt_soil(2)
    call One_Out ('ppt_soil', 8, Vals)
    Vals(1) = Exceed%runoff(1)  ;  Vals(2) = Exceed%runoff(2)
    call One_Out ('runoff', 6, Vals)
    Vals(1) = Exceed%ratio_water_pet(1)  ;  Vals(2) = Exceed%ratio_water_pet(2)
    call One_Out ('ratio_water_pet', 15, Vals)
!    Vals(1) = Exceed%co2_value(1)  ;  Vals(2) = Exceed%co2_value(2)
!    call One_Out ('co2_value', 9, Vals)
    Vals(1) = Exceed%pet_top_soil(1)  ;  Vals(2) = Exceed%pet_top_soil(2)
    call One_Out ('pet_top_soil', 12, Vals)
    Vals(1) = Exceed%n_leached(1)  ;  Vals(2) = Exceed%n_leached(2)
    call One_Out ('n_leached', 9, Vals)
    Vals(1) = Exceed%asmos(1)  ;  Vals(2) = Exceed%asmos(2)
    call One_Out ('asmos', 5, Vals)
    Vals(1) = Exceed%amov(1)  ;  Vals(2) = Exceed%amov(2)
    call One_Out ('amov', 4, Vals)
    Vals(1) = Exceed%storm_flow(1)  ;  Vals(2) = Exceed%storm_flow(2)
    call One_Out ('storm_flow', 10, Vals)
    Vals(1) = Exceed%holding_tank(1)  ;  Vals(2) = Exceed%holding_tank(2)
    call One_Out ('holding_tank', 12, Vals)
    Vals(1) = Exceed%transpiration(1)  ;  Vals(2) = Exceed%transpiration(2)
    call One_Out ('transpiration', 13, Vals)
    Vals(1) = Exceed%relative_water_content(1)  ;  Vals(2) = Exceed%relative_water_content(2)
    call One_Out ('relative_water_content', 22, Vals)
    Vals(1) = Exceed%water_available(1)  ;  Vals(2) = Exceed%water_available(2)
    call One_Out ('water_available', 15, Vals)
    Vals(1) = Exceed%annual_evapotranspiration(1)  ;  Vals(2) = Exceed%annual_evapotranspiration(2)
    call One_Out ('annual_evapotranspiration', 25, Vals)
    Vals(1) = Exceed%total_aground_live_biomass(1)  ;  Vals(2) = Exceed%total_aground_live_biomass(2)
    call One_Out ('total_aground_live_biomass', 26, Vals)
    Vals(1) = Exceed%total_litter_carbon(1)  ;  Vals(2) = Exceed%total_litter_carbon(2)
    call One_Out ('total_litter_carbon', 19, Vals)
    Vals(1) = Exceed%total_litter_nitrogen(1)  ;  Vals(2) = Exceed%total_litter_nitrogen(2)
    call One_Out ('total_litter_nitrogen', 21, Vals)
    Vals(1) = Exceed%root_shoot_ratio(1)  ;  Vals(2) = Exceed%root_shoot_ratio(2)
    call One_Out ('root_shoot_ratio', 16, Vals)
    Vals(1) = Exceed%tree_basal_area(1)  ;  Vals(2) = Exceed%tree_basal_area(2)
    call One_Out ('tree_basal_area', 15, Vals)
    Vals(1) = Exceed%soil_surface_temperature(1)  ;  Vals(2) = Exceed%soil_surface_temperature(2)
    call One_Out ('soil_surface_temperature', 24, Vals)
    Vals(1) = Exceed%mineral_nitrogen(1)  ;  Vals(2) = Exceed%mineral_nitrogen(2)
    call One_Out ('mineral_nitrogen', 16, Vals)
    Vals(1) = Exceed%field_capacity(1)  ;  Vals(2) = Exceed%field_capacity(2)
    call One_Out ('field_capacity', 14, Vals)
    Vals(1) = Exceed%wilting_point(1)  ;  Vals(2) = Exceed%wilting_point(2)
    call One_Out ('wilting_point', 13, Vals)
    Vals(1) = Exceed%soil_total_carbon(1)  ;  Vals(2) = Exceed%soil_total_carbon(2)
    call One_Out ('soil_total_carbon', 17, Vals)
    Vals(1) = Exceed%tree_carbon(1)  ;  Vals(2) = Exceed%tree_carbon(2)
    call One_Out ('tree_carbon', 11, Vals)
    Vals(1) = Exceed%tree_nitrogen(1)  ;  Vals(2) = Exceed%tree_nitrogen(2)
    call One_Out ('tree_nitrogen', 13, Vals)
    Vals(1) = Exceed%shrub_carbon(1)  ;  Vals(2) = Exceed%shrub_carbon(2)
    call One_Out ('shrub_carbon', 12, Vals)
    Vals(1) = Exceed%shrub_nitrogen(1)  ;  Vals(2) = Exceed%shrub_nitrogen(2)
    call One_Out ('shrub_nitrogen', 14, Vals)
    Vals(1) = Exceed%carbon_nitrogen_ratio(1)  ;  Vals(2) = Exceed%carbon_nitrogen_ratio(2)
    call One_Out ('carbon_nitrogen_ratio', 21, Vals)
    Vals(1) = Exceed%fast_soil_carbon(1)  ;  Vals(2) = Exceed%fast_soil_carbon(2)
    call One_Out ('fast_soil_carbon', 16, Vals)
    Vals(1) = Exceed%intermediate_soil_carbon(1)  ;  Vals(2) = Exceed%intermediate_soil_carbon(2)
    call One_Out ('intermediate_soil_carbon', 24, Vals)
    Vals(1) = Exceed%passive_soil_carbon(1)  ;  Vals(2) = Exceed%passive_soil_carbon(2)
    call One_Out ('passive_soil_carbon', 19, Vals)
    Vals(1) = Exceed%fast_soil_nitrogen(1)  ;  Vals(2) = Exceed%fast_soil_nitrogen(2)
    call One_Out ('fast_soil_nitrogen', 18, Vals)
    Vals(1) = Exceed%intermediate_soil_nitrogen(1)  ;  Vals(2) = Exceed%intermediate_soil_nitrogen(2)
    call One_Out ('intermediate_soil_nitrogen', 26, Vals)
    Vals(1) = Exceed%passive_soil_nitrogen(1)  ;  Vals(2) = Exceed%passive_soil_nitrogen(2)
    call One_Out ('passive_soil_nitrogen', 21, Vals)
    Vals(1) = Exceed%potential_production(1)  ;  Vals(2) = Exceed%potential_production(2)
    call One_Out ('potential_production', 20, Vals)
    Vals(1) = Exceed%belowground_pot_production(1)  ;  Vals(2) = Exceed%belowground_pot_production(2)
    call One_Out ('belowground_pot_production', 26, Vals)
    Vals(1) = Exceed%aboveground_pot_production(1)  ;  Vals(2) = Exceed%aboveground_pot_production(2)
    call One_Out ('aboveground_pot_production', 26, Vals)
    Vals(1) = Exceed%total_pot_production(1)  ;  Vals(2) = Exceed%total_pot_production(2)
    call One_Out ('total_pot_production', 20, Vals)
    Vals(1) = Exceed%co2_effect_on_production(1)  ;  Vals(2) = Exceed%co2_effect_on_production(2)
    call One_Out ('co2_effect_on_production', 24, Vals)
    Vals(1) = Exceed%total_pot_prod_limited_by_n(1)  ;  Vals(2) = Exceed%total_pot_prod_limited_by_n(2)
    call One_Out ('total_pot_prod_limited_by_n', 27, Vals)
    Vals(1) = Exceed%monthly_net_primary_production(1)  ;  Vals(2) = Exceed%monthly_net_primary_production(2)
    call One_Out ('monthly_net_primary_production', 30, Vals)
    Vals(1) = Exceed%fraction_live_removed_grazing(1)  ;  Vals(2) = Exceed%fraction_live_removed_grazing(2)
    call One_Out ('fraction_live_removed_grazing', 29, Vals)
    Vals(1) = Exceed%fraction_dead_removed_grazing(1)  ;  Vals(2) = Exceed%fraction_dead_removed_grazing(2)
    call One_Out ('fraction_dead_removed_grazing', 29, Vals)
    Vals(1) = Exceed%temp_effect_on_decomp(1)  ;  Vals(2) = Exceed%temp_effect_on_decomp(2)
    call One_Out ('temp_effect_on_decomp', 21, Vals)
    Vals(1) = Exceed%water_effect_on_decomp(1)  ;  Vals(2) = Exceed%water_effect_on_decomp(2)
    call One_Out ('water_effect_on_decomp', 22, Vals)
    Vals(1) = Exceed%anerobic_effect_on_decomp(1)  ;  Vals(2) = Exceed%anerobic_effect_on_decomp(2)
    call One_Out ('anerobic_effect_on_decomp', 25, Vals)
    Vals(1) = Exceed%all_effects_on_decomp(1)  ;  Vals(2) = Exceed%all_effects_on_decomp(2)
    call One_Out ('all_effects_on_decomp', 21, Vals)
    Vals(1) = Exceed%dead_fine_root_carbon(1)  ;  Vals(2) = Exceed%dead_fine_root_carbon(2)
    call One_Out ('dead_fine_root_carbon', 21, Vals)
    Vals(1) = Exceed%dead_fine_root_nitrogen(1)  ;  Vals(2) = Exceed%dead_fine_root_nitrogen(2)
    call One_Out ('dead_fine_root_nitrogen', 23, Vals)
    Vals(1) = Exceed%dead_standing_carbon(1)  ;  Vals(2) = Exceed%dead_standing_carbon(2)
    call One_Out ('dead_standing_carbon', 20, Vals)
    Vals(1) = Exceed%dead_standing_nitrogen(1)  ;  Vals(2) = Exceed%dead_standing_nitrogen(2)
    call One_Out ('dead_standing_nitrogen', 22, Vals)
    Vals(1) = Exceed%dead_seed_carbon(1)  ;  Vals(2) = Exceed%dead_seed_carbon(2)
    call One_Out ('dead_seed_carbon', 16, Vals)
    Vals(1) = Exceed%dead_seed_nitrogen(1)  ;  Vals(2) = Exceed%dead_seed_nitrogen(2)
    call One_Out ('dead_seed_nitrogen', 18, Vals)
    Vals(1) = Exceed%dead_leaf_carbon(1)  ;  Vals(2) = Exceed%dead_leaf_carbon(2)
    call One_Out ('dead_leaf_carbon', 16, Vals)
    Vals(1) = Exceed%dead_leaf_nitrogen(1)  ;  Vals(2) = Exceed%dead_leaf_nitrogen(2)
    call One_Out ('dead_leaf_nitrogen', 18, Vals)
    Vals(1) = Exceed%dead_fine_branch_carbon(1)  ;  Vals(2) = Exceed%dead_fine_branch_carbon(2)
    call One_Out ('dead_fine_branch_carbon', 23, Vals)
    Vals(1) = Exceed%dead_total_fine_branch_carbon(1)  ;  Vals(2) = Exceed%dead_total_fine_branch_carbon(2)
    call One_Out ('dead_total_fine_branch_carbon', 29, Vals)
    Vals(1) = Exceed%dead_fine_branch_nitrogen(1)  ;  Vals(2) = Exceed%dead_fine_branch_nitrogen(2)
    call One_Out ('dead_fine_branch_nitrogen', 25, Vals)
    Vals(1) = Exceed%dead_total_fine_branch_nitrogen(1)  ;  Vals(2) = Exceed%dead_total_fine_branch_nitrogen(2)
    call One_Out ('dead_total_fine_branch_nitrogen', 31, Vals)
    Vals(1) = Exceed%dead_coarse_root_carbon(1)  ;  Vals(2) = Exceed%dead_coarse_root_carbon(2)
    call One_Out ('dead_coarse_root_carbon', 23, Vals)
    Vals(1) = Exceed%dead_total_coarse_root_carbon(1)  ;  Vals(2) = Exceed%dead_total_coarse_root_carbon(2)
    call One_Out ('dead_total_coarse_root_carbon', 29, Vals)
    Vals(1) = Exceed%dead_coarse_root_nitrogen(1)  ;  Vals(2) = Exceed%dead_coarse_root_nitrogen(2)
    call One_Out ('dead_coarse_root_nitrogen', 25, Vals)
    Vals(1) = Exceed%dead_total_coarse_root_nitrogen(1)  ;  Vals(2) = Exceed%dead_total_coarse_root_nitrogen(2)
    call One_Out ('dead_total_coarse_root_nitrogen', 31, Vals)
    Vals(1) = Exceed%dead_coarse_branch_carbon(1)  ;  Vals(2) = Exceed%dead_coarse_branch_carbon(2)
    call One_Out ('dead_coarse_branch_carbon', 25, Vals)
    Vals(1) = Exceed%dead_total_coarse_branch_carbon(1)  ;  Vals(2) = Exceed%dead_total_coarse_branch_carbon(2)
    call One_Out ('dead_total_coarse_branch_carbon', 31, Vals)
    Vals(1) = Exceed%dead_coarse_branch_nitrogen(1)  ;  Vals(2) = Exceed%dead_coarse_branch_nitrogen(2)
    call One_Out ('dead_coarse_branch_nitrogen', 27, Vals)
    Vals(1) = Exceed%dead_total_coarse_branch_nitrogen(1)  ;  Vals(2) = Exceed%dead_total_coarse_branch_nitrogen(2)
    call One_Out ('dead_total_coarse_branch_nitrogen', 33, Vals)
    Vals(1) = Exceed%lignin_fine_root(1)  ;  Vals(2) = Exceed%lignin_fine_root(2)
    call One_Out ('lignin_fine_root', 16, Vals)
    Vals(1) = Exceed%lignin_coarse_root(1)  ;  Vals(2) = Exceed%lignin_coarse_root(2)
    call One_Out ('lignin_coarse_root', 18, Vals)
    Vals(1) = Exceed%lignin_fine_branch(1)  ;  Vals(2) = Exceed%lignin_fine_branch(2)
    call One_Out ('lignin_fine_branch', 18, Vals)
    Vals(1) = Exceed%lignin_coarse_branch(1)  ;  Vals(2) = Exceed%lignin_coarse_branch(2)
    call One_Out ('lignin_coarse_branch', 20, Vals)
    Vals(1) = Exceed%lignin_leaf(1)  ;  Vals(2) = Exceed%lignin_leaf(2)
    call One_Out ('lignin_leaf', 11, Vals)
    Vals(1) = Exceed%plant_lignin_fraction(1)  ;  Vals(2) = Exceed%plant_lignin_fraction(2)
    call One_Out ('plant_lignin_fraction', 21, Vals)
    Vals(1) = Exceed%litter_structural_carbon(1)  ;  Vals(2) = Exceed%litter_structural_carbon(2)
    call One_Out ('litter_structural_carbon', 24, Vals)
    Vals(1) = Exceed%litter_metabolic_carbon(1)  ;  Vals(2) = Exceed%litter_metabolic_carbon(2)
    call One_Out ('litter_metabolic_carbon', 23, Vals)
    Vals(1) = Exceed%litter_structural_nitrogen(1)  ;  Vals(2) = Exceed%litter_structural_nitrogen(2)
    call One_Out ('litter_structural_nitrogen', 26, Vals)
    Vals(1) = Exceed%litter_metabolic_nitrogen(1)  ;  Vals(2) = Exceed%litter_metabolic_nitrogen(2)
    call One_Out ('litter_metabolic_nitrogen', 25, Vals)
    Vals(1) = Exceed%tnetmin(1)  ;  Vals(2) = Exceed%tnetmin(2)
    call One_Out ('tnetmin', 7, Vals)
    Vals(1) = Exceed%tminup(1)  ;  Vals(2) = Exceed%tminup(2)
    call One_Out ('tminup', 6, Vals)
    Vals(1) = Exceed%grossmin(1)  ;  Vals(2) = Exceed%grossmin(2)
    call One_Out ('grossmin', 8, Vals)
    Vals(1) = Exceed%volitn(1)  ;  Vals(2) = Exceed%volitn(2)
    call One_Out ('volitn', 6, Vals)
    Vals(1) = Exceed%fixnit(1)  ;  Vals(2) = Exceed%fixnit(2)
    call One_Out ('fixnit', 6, Vals)
    Vals(1) = Exceed%runoffn(1)  ;  Vals(2) = Exceed%runoffn(2)
    call One_Out ('runoffn', 7, Vals)
    Vals(1) = Exceed%e_up(1)  ;  Vals(2) = Exceed%e_up(2)
    call One_Out ('e_up', 4, Vals)
    Vals(1) = Exceed%volatized_n(1)  ;  Vals(2) = Exceed%volatized_n(2)
    call One_Out ('volatized_n', 11, Vals)
    Vals(1) = Exceed%maintain_respiration(1)  ;  Vals(2) = Exceed%maintain_respiration(2)
    call One_Out ('maintain_respiration', 20, Vals)
    Vals(1) = Exceed%phenology(1)  ;  Vals(2) = Exceed%phenology(2)
    call One_Out ('phenology', 19, Vals)
    Vals(1) = Exceed%fine_root_carbon(1)  ;  Vals(2) = Exceed%fine_root_carbon(2)
    call One_Out ('fine_root_carbon', 16, Vals)
    Vals(1) = Exceed%fine_root_nitrogen(1)  ;  Vals(2) = Exceed%fine_root_nitrogen(2)
    call One_Out ('fine_root_nitrogen', 18, Vals)
    Vals(1) = Exceed%seed_carbon(1)  ;  Vals(2) = Exceed%seed_carbon(2)
    call One_Out ('seed_carbon', 11, Vals)
    Vals(1) = Exceed%seed_nitrogen(1)  ;  Vals(2) = Exceed%seed_nitrogen(2)
    call One_Out ('seed_nitrogen', 13, Vals)
    Vals(1) = Exceed%leaf_carbon(1)  ;  Vals(2) = Exceed%leaf_carbon(2)
    call One_Out ('leaf_carbon', 11, Vals)
    Vals(1) = Exceed%leaf_nitrogen(1)  ;  Vals(2) = Exceed%leaf_nitrogen(2)
    call One_Out ('leaf_nitrogen', 13, Vals)
    Vals(1) = Exceed%fine_branch_carbon(1)  ;  Vals(2) = Exceed%fine_branch_carbon(2)
    call One_Out ('fine_branch_carbon', 18, Vals)
    Vals(1) = Exceed%fine_branch_nitrogen(1)  ;  Vals(2) = Exceed%fine_branch_nitrogen(2)
    call One_Out ('fine_branch_nitrogen', 20, Vals)
    Vals(1) = Exceed%coarse_root_carbon(1)  ;  Vals(2) = Exceed%coarse_root_carbon(2)
    call One_Out ('coarse_root_carbon', 18, Vals)
    Vals(1) = Exceed%coarse_root_nitrogen(1)  ;  Vals(2) = Exceed%coarse_root_nitrogen(2)
    call One_Out ('coarse_root_nitrogen', 20, Vals)
    Vals(1) = Exceed%coarse_branch_carbon(1)  ;  Vals(2) = Exceed%coarse_branch_carbon(2)
    call One_Out ('coarse_branch_carbon', 20, Vals)
    Vals(1) = Exceed%coarse_branch_nitrogen(1)  ;  Vals(2) = Exceed%coarse_branch_nitrogen(2)
    call One_Out ('coarse_branch_nitrogen', 22, Vals)
    Vals(1) = Exceed%stored_nitrogen(1)  ;  Vals(2) = Exceed%stored_nitrogen(2)
    call One_Out ('stored_nitrogen', 15, Vals)
    Vals(1) = Exceed%plant_nitrogen_fixed(1)  ;  Vals(2) = Exceed%plant_nitrogen_fixed(2)
    call One_Out ('plant_nitrogen_fixed', 20, Vals)
    Vals(1) = Exceed%nitrogen_fixed(1)  ;  Vals(2) = Exceed%nitrogen_fixed(2)
    call One_Out ('nitrogen_fixed', 14, Vals)
    Vals(1) = Exceed%respiration_flows(1)  ;  Vals(2) = Exceed%respiration_flows(2)
    call One_Out ('respiration_flows', 17, Vals)
    Vals(1) = Exceed%respiration_annual(1)  ;  Vals(2) = Exceed%respiration_annual(2)
    call One_Out ('respiration_annual', 18, Vals)
    Vals(1) = Exceed%carbon_source_sink(1)  ;  Vals(2) = Exceed%carbon_source_sink(2)
    call One_Out ('carbon_source_sink', 18, Vals)
    Vals(1) = Exceed%nitrogen_source_sink(1)  ;  Vals(2) = Exceed%nitrogen_source_sink(2)
    call One_Out ('nitrogen_source_sink', 20, Vals)
    Vals(1) = Exceed%carbon_allocation(1)  ;  Vals(2) = Exceed%carbon_allocation(2)
    call One_Out ('carbon_allocation', 17, Vals)
    Vals(1) = Exceed%optimum_leaf_area_index(1)  ;  Vals(2) = Exceed%optimum_leaf_area_index(2)
    call One_Out ('optimum_leaf_area_index', 23, Vals)
    Vals(1) = Exceed%leaf_area_index(1)  ;  Vals(2) = Exceed%leaf_area_index(2)
    call One_Out ('leaf_area_index', 15, Vals)
    Vals(1) = Exceed%water_function(1)  ;  Vals(2) = Exceed%water_function(2)
    call One_Out ('water_function', 14, Vals)
  else
    write(ECHO_FILE,*) 'Unable to open the exceed.gof'
  end if
  close(SHORT_USE_FILE)

  open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//RUN_TIME_NAME, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)
  if (ioerr == 0) then
    call cpu_time (end_time)
    write(SHORT_USE_FILE,'(A44, F10.1)') 'The elapsed time, in seconds:             ', end_time - start_time
    write(SHORT_USE_FILE,'(A44, F11.2)') 'The elapsed time, in minutes:             ', (end_time - start_time) / 60.0
    write(SHORT_USE_FILE,'(A44, F12.3)') 'The elapsed time, in hours:               ', ((end_time - start_time) / 60.0) / 60.0
    write(SHORT_USE_FILE,*) '  '
    write(SHORT_USE_FILE,'(A44, I10)') 'The number of rangeland cells processed:  ', range_cells
    yr_cnt = ( Sim_Parm%end_Yr - Sim_Parm%start_yr ) + 1
    write(SHORT_USE_FILE,'(A44, I5)') 'The number of years simulated:            ', yr_cnt
    write(SHORT_USE_FILE,'(A44, F16.7)') 'The seconds to simulate a cell-month:     ', (end_time - start_time) / &
                                                                                       ( yr_cnt * 12 * range_cells )
    write(SHORT_USE_FILE,'(A44, F160.7)') 'The seconds to simulate 1000 cell-months: ', ( (end_time - start_time) / &
                                                                         ( yr_cnt * 12 * range_cells ) ) * 1000.0
  else
    write(ECHO_FILE,*) 'Unable to open the run_time.gof'
  endif

  close(SHORT_USE_FILE)
  close(ECHO_FILE)

end subroutine


subroutine Progress
!**** Report progress of the simulation.  Not using CPU_TIME here, using real time so as to be useful to the user.
!****
!**** R. Boone   Last modified: April 24, 2011
  use Parameter_Vars
  use Structures
  implicit none

  real  progress_to_date, months_to_date, total_months
  months_to_date = ( ( year - Sim_Parm%start_yr ) * 12 ) + month
  total_months = ( ( Sim_Parm%end_yr - Sim_Parm%start_yr ) + 1 ) * 12
  progress_to_date = ( months_to_date / total_months ) * 100.0
  write(*, 1099) 'Year: ', year, ' Month: ', month, '  Progress: ', int(progress_to_date), '%'
1099  format(A6,I4,A8,I2,A12,I3,A2)


end subroutine



subroutine One_Out (out_string, str_n, Vals)
  ! Output a line in a file that reports the number of times parameter values were exceeded.
  ! R. B. Boone          Last modified:  March 23, 2011
  use Parameter_Vars
  use Structures
  implicit none

  integer Vals(2), str_n
  character(40) out_String, out2

  out2 = adjustl(out_string(1:str_n))
  write(SHORT_USE_FILE,'(A40,A4,I10,A4,I10)') out2, '    ', Vals(1), '    ', Vals(2)

end subroutine


real function linear (x, data_val, imx)
  !**** A linear interpolation routine, from Savanna ALINT.
  ! X is the value that needs a Y
  ! Data_vals are pairs of values that define the relationship, x1, y1, x2, y2, etc.
  ! imx is the number of pairs given, often 2
  ! Rewritten March 12, 2013 to ensure that values are going in the placed needed.
  ! Typo fixed that had do i = l versus do i = 1 prior to the exit call.   
  ! Last edited:  June 17, 2014.
    real data_val(10)                                          ! Hardwired maximum vector size, the largest is associated with phenology
    real data_v(2,imx)
    integer k, m, n

    do m=1,imx
      n = m * 2
      data_v(1, m) = data_val(n - 1)
      data_v(2, m) = data_val(n)
    end do

    if (x .le. data_v(1,1)) then
      linear = data_v(2,1)
      return
    end if
    if (x .ge. data_v(1,imx)) then
      linear = data_v(2,imx)
      return
    end if

    do i = 1, imx - 1
      if (x .le. data_v(1,i+1)) then
        k = i
        exit
      end if
    end do

    linear = data_v(2,k) + (data_v(2,k+1)-data_v(2,k))/(data_v(1,k+1)-data_v(1,k))*(x-data_v(1,k))

end function


subroutine Each_Month
!**** Do processing steps that must be done each month.  The main steps are a long series of tests to ensure that
!**** values aren't exceeding a very large value, or moving negative.  Errors will cause tallying of counts of errors,
!**** both spatially and per entry.  That said, they won't be stored spatially for each individual entry, as that
!**** would almost double memory.
!****
!**** This routine includes a simple assignment of grazing fraction.  That logic is placed here to allow for it to
!**** be made more dynamic in the future.
!****
!**** R. Boone   Last modified: February 15, 2013
!****                           February 15, 2013 changes incorporate monthly clear-out of dead material stores
  use Parameter_Vars
  use Structures
  implicit none

  real    live_carbon, dead_carbon
  integer i, j, icell, iunit

 do icell = 1, range_cells

  !  last_month_day_length
  if (Rng(icell)%last_month_day_length .gt. V_LARGE) then
    Rng(icell)%last_month_day_length = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%last_month_day_length(1) = Exceed%last_month_day_length(1) + 1
  end if
  if (Rng(icell)%last_month_day_length .lt.   0.0  ) then
    Rng(icell)%last_month_day_length = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%last_month_day_length(2) = Exceed%last_month_day_length(2) + 1
  end if

  !  day_length
  if (Rng(icell)%day_length .gt. V_LARGE) then
    Rng(icell)%day_length = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%day_length(1) = Exceed%day_length(1) + 1
  end if
  if (Rng(icell)%day_length .lt.   0.0  ) then
    Rng(icell)%day_length = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%day_length(2) = Exceed%day_length(2) + 1
  end if

  !  heat_accumulation
  if (Rng(icell)%heat_accumulation .gt. V_LARGE) then
    Rng(icell)%heat_accumulation = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%heat_accumulation(1) = Exceed%heat_accumulation(1) + 1
  end if
  if (Rng(icell)%heat_accumulation .lt.   0.0  ) then
    Rng(icell)%heat_accumulation = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%heat_accumulation(2) = Exceed%heat_accumulation(2) + 1
  end if

  !  facet_cover
  do i = 1, FACETS
    if (Rng(icell)%facet_cover(i) .gt. V_LARGE) then
      Rng(icell)%facet_cover(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%facet_cover(1) = Exceed%facet_cover(1) + 1
    end if
    if (Rng(icell)%facet_cover(i) .lt.   0.0  ) then
      Rng(icell)%facet_cover(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%facet_cover(2) = Exceed%facet_cover(2) + 1
    end if
  end do

  !  total_population
  do i = 1, V_LYRS
    if (Rng(icell)%total_population(i) .gt. V_LARGE) then
      Rng(icell)%total_population(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%total_population(1) = Exceed%total_population(1) + 1
    end if
    if (Rng(icell)%total_population(i) .lt.   0.0  ) then
      Rng(icell)%total_population(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%total_population(2) = Exceed%total_population(2) + 1
    end if
  end do

  !  bare_cover
  if (Rng(icell)%bare_cover .gt. V_LARGE) then
    Rng(icell)%bare_cover = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%bare_cover(1) = Exceed%bare_cover(1) + 1
  end if
  if (Rng(icell)%bare_cover .lt.   0.0  ) then
    Rng(icell)%bare_cover = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%bare_cover(2) = Exceed%bare_cover(2) + 1
  end if

  !  prop_annual_decid
  do i = 1, FACETS
    if (Rng(icell)%prop_annual_decid(i) .gt. V_LARGE) then
      Rng(icell)%prop_annual_decid(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%prop_annual_decid(1) = Exceed%prop_annual_decid(1) + 1
    end if
    if (Rng(icell)%prop_annual_decid(i) .lt.   0.0  ) then
      Rng(icell)%prop_annual_decid(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%prop_annual_decid(2) = Exceed%prop_annual_decid(2) + 1
    end if
  end do

  !  pot_evap
  if (Rng(icell)%pot_evap .gt. V_LARGE) then
    Rng(icell)%pot_evap = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%pot_evap(1) = Exceed%pot_evap(1) + 1
  end if
  if (Rng(icell)%pot_evap .lt.   0.0  ) then
    Rng(icell)%pot_evap = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%pot_evap(2) = Exceed%pot_evap(2) + 1
  end if

  !  evaporation
  if (Rng(icell)%evaporation .gt. V_LARGE) then
    Rng(icell)%evaporation = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%evaporation(1) = Exceed%evaporation(1) + 1
  end if
  if (Rng(icell)%evaporation .lt.   0.0  ) then
    Rng(icell)%evaporation = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%evaporation(2) = Exceed%evaporation(2) + 1
  end if

  !  snow
  if (Rng(icell)%snow .gt. V_LARGE) then
    Rng(icell)%snow = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%snow(1) = Exceed%snow(1) + 1
  end if
  if (Rng(icell)%snow .lt.   0.0  ) then
    Rng(icell)%snow = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%snow(2) = Exceed%snow(2) + 1
  end if

  !  snow_liquid
  if (Rng(icell)%snow_liquid .gt. V_LARGE) then
    Rng(icell)%snow_liquid = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%snow_liquid(1) = Exceed%snow_liquid(1) + 1
  end if
  if (Rng(icell)%snow_liquid .lt.   0.0  ) then
    Rng(icell)%snow_liquid = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%snow_liquid(2) = Exceed%snow_liquid(2) + 1
  end if

  !  melt
  if (Rng(icell)%melt .gt. V_LARGE) then
    Rng(icell)%melt = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%melt(1) = Exceed%melt(1) + 1
  end if
  if (Rng(icell)%melt .lt.   0.0  ) then
    Rng(icell)%melt = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%melt(2) = Exceed%melt(2) + 1
  end if

  !  pet_remaining
  if (Rng(icell)%pet_remaining .gt. V_LARGE) then
    Rng(icell)%pet_remaining = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%pet_remaining(1) = Exceed%pet_remaining(1) + 1
  end if
  if (Rng(icell)%pet_remaining .lt.   0.0  ) then
    Rng(icell)%pet_remaining = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%pet_remaining(2) = Exceed%pet_remaining(2) + 1
  end if

  !  ppt_soil
  if (Rng(icell)%ppt_soil .gt. V_LARGE) then
    Rng(icell)%ppt_soil = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%ppt_soil(1) = Exceed%ppt_soil(1) + 1
  end if
  if (Rng(icell)%ppt_soil .lt.   0.0  ) then
    Rng(icell)%ppt_soil = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%ppt_soil(2) = Exceed%ppt_soil(2) + 1
  end if

  !  runoff
  if (Rng(icell)%runoff .gt. V_LARGE) then
    Rng(icell)%runoff = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%runoff(1) = Exceed%runoff(1) + 1
  end if
  if (Rng(icell)%runoff .lt.   0.0  ) then
    Rng(icell)%runoff = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%runoff(2) = Exceed%runoff(2) + 1
  end if

  !  ratio_water_pet
  if (Rng(icell)%ratio_water_pet .gt. V_LARGE) then
    Rng(icell)%ratio_water_pet = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%ratio_water_pet(1) = Exceed%ratio_water_pet(1) + 1
  end if
  if (Rng(icell)%ratio_water_pet .lt.   0.0  ) then
    Rng(icell)%ratio_water_pet = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%ratio_water_pet(2) = Exceed%ratio_water_pet(2) + 1
  end if

!  !  co2_value
!  if (Rng(icell)%co2_value .gt. V_LARGE) then
!    Rng(icell)%co2_value = V_LARGE
!    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
!    Exceed%co2_value(1) = Exceed%co2_value(1) + 1
!  end if
!  if (Rng(icell)%co2_value .lt.   0.0  ) then
!    Rng(icell)%co2_value = 0.0
!    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
!    Exceed%co2_value(2) = Exceed%co2_value(2) + 1
!  end if

  !  pet_top_soil
  if (Rng(icell)%pet_top_soil .gt. V_LARGE) then
    Rng(icell)%pet_top_soil = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%pet_top_soil(1) = Exceed%pet_top_soil(1) + 1
  end if
  if (Rng(icell)%pet_top_soil .lt.   0.0  ) then
    Rng(icell)%pet_top_soil = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%pet_top_soil(2) = Exceed%pet_top_soil(2) + 1
  end if

  !  n_leached
  do i = 1, SOIL_LAYERS
    if (Rng(icell)%n_leached(i) .gt. V_LARGE) then
      Rng(icell)%n_leached(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%n_leached(1) = Exceed%n_leached(1) + 1
    end if
    if (Rng(icell)%n_leached(i) .lt.   0.0  ) then
      Rng(icell)%n_leached(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%n_leached(2) = Exceed%n_leached(2) + 1
    end if
  end do

  !  asmos
  do i = 1, SOIL_LAYERS
    if (Rng(icell)%asmos(i) .gt. V_LARGE) then
      Rng(icell)%asmos(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%asmos(1) = Exceed%asmos(1) + 1
    end if
    if (Rng(icell)%asmos(i) .lt.   0.0  ) then
      Rng(icell)%asmos(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%asmos(2) = Exceed%asmos(2) + 1
    end if
  end do

  !  amov
  do i = 1, SOIL_LAYERS
    if (Rng(icell)%amov(i) .gt. V_LARGE) then
      Rng(icell)%amov(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%amov(1) = Exceed%amov(1) + 1
    end if
    if (Rng(icell)%amov(i) .lt.   0.0  ) then
      Rng(icell)%amov(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%amov(2) = Exceed%amov(2) + 1
    end if
  end do

  !  storm_flow
  if (Rng(icell)%storm_flow .gt. V_LARGE) then
    Rng(icell)%storm_flow = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%storm_flow(1) = Exceed%storm_flow(1) + 1
  end if
  if (Rng(icell)%storm_flow .lt.   0.0  ) then
    Rng(icell)%storm_flow = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%storm_flow(2) = Exceed%storm_flow(2) + 1
  end if

  !  holding_tank
  if (Rng(icell)%holding_tank .gt. V_LARGE) then
    Rng(icell)%holding_tank = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%holding_tank(1) = Exceed%holding_tank(1) + 1
  end if
  if (Rng(icell)%holding_tank .lt.   0.0  ) then
    Rng(icell)%holding_tank = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%holding_tank(2) = Exceed%holding_tank(2) + 1
  end if

  !  transpiration
  if (Rng(icell)%transpiration .gt. V_LARGE) then
    Rng(icell)%transpiration = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%transpiration(1) = Exceed%transpiration(1) + 1
  end if
  if (Rng(icell)%transpiration .lt.   0.0  ) then
    Rng(icell)%transpiration = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%transpiration(2) = Exceed%transpiration(2) + 1
  end if

  !  relative_water_content
  do i = 1, SOIL_LAYERS
    if (Rng(icell)%relative_water_content(i) .gt. V_LARGE) then
      Rng(icell)%relative_water_content(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%relative_water_content(1) = Exceed%relative_water_content(1) + 1
    end if
    if (Rng(icell)%relative_water_content(i) .lt.   0.0  ) then
      Rng(icell)%relative_water_content(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%relative_water_content(2) = Exceed%relative_water_content(2) + 1
    end if
  end do

  !  water_available
  do i = 1, 3
    if (Rng(icell)%water_available(i) .gt. V_LARGE) then
      Rng(icell)%water_available(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%water_available(1) = Exceed%water_available(1) + 1
    end if
    if (Rng(icell)%water_available(i) .lt.   0.0  ) then
      Rng(icell)%water_available(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%water_available(2) = Exceed%water_available(2) + 1
    end if
  end do

  !  annual_evapotranspiration
  if (Rng(icell)%annual_evapotranspiration .gt. V_LARGE) then
    Rng(icell)%annual_evapotranspiration = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%annual_evapotranspiration(1) = Exceed%annual_evapotranspiration(1) + 1
  end if
  if (Rng(icell)%annual_evapotranspiration .lt.   0.0  ) then
    Rng(icell)%annual_evapotranspiration = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%annual_evapotranspiration(2) = Exceed%annual_evapotranspiration(2) + 1
  end if

  !  total_aground_live_biomass
  if (Rng(icell)%total_aground_live_biomass .gt. V_LARGE) then
    Rng(icell)%total_aground_live_biomass = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%total_aground_live_biomass(1) = Exceed%total_aground_live_biomass(1) + 1
  end if
  if (Rng(icell)%total_aground_live_biomass .lt.   0.0  ) then
    Rng(icell)%total_aground_live_biomass = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%total_aground_live_biomass(2) = Exceed%total_aground_live_biomass(2) + 1
  end if

  !  total_bground_live_biomass
  if (Rng(icell)%total_bground_live_biomass .gt. V_LARGE) then
    Rng(icell)%total_bground_live_biomass = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%total_bground_live_biomass(1) = Exceed%total_bground_live_biomass(1) + 1
  end if
  if (Rng(icell)%total_bground_live_biomass .lt.   0.0  ) then
    Rng(icell)%total_bground_live_biomass = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%total_bground_live_biomass(2) = Exceed%total_bground_live_biomass(2) + 1
  end if

  !  total_litter_carbon
  do i = 1, 2
    if (Rng(icell)%total_litter_carbon(i) .gt. V_LARGE) then
      Rng(icell)%total_litter_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%total_litter_carbon(1) = Exceed%total_litter_carbon(1) + 1
    end if
    if (Rng(icell)%total_litter_carbon(i) .lt.   0.0  ) then
      Rng(icell)%total_litter_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%total_litter_carbon(2) = Exceed%total_litter_carbon(2) + 1
    end if
  end do

  !  total_litter_nitrogen
  do i = 1, 2
    if (Rng(icell)%total_litter_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%total_litter_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%total_litter_nitrogen(1) = Exceed%total_litter_nitrogen(1) + 1
    end if
    if (Rng(icell)%total_litter_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%total_litter_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%total_litter_nitrogen(2) = Exceed%total_litter_nitrogen(2) + 1
    end if
  end do

  !  root_shoot_ratio
  do i = 1, FACETS
    if (Rng(icell)%root_shoot_ratio(i) .gt. V_LARGE) then
      Rng(icell)%root_shoot_ratio(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%root_shoot_ratio(1) = Exceed%root_shoot_ratio(1) + 1
    end if
    if (Rng(icell)%root_shoot_ratio(i) .lt.   0.0  ) then
      Rng(icell)%root_shoot_ratio(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%root_shoot_ratio(2) = Exceed%root_shoot_ratio(2) + 1
    end if
  end do

  !  tree_basal_area
  if (Rng(icell)%tree_basal_area .gt. V_LARGE) then
    Rng(icell)%tree_basal_area = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%tree_basal_area(1) = Exceed%tree_basal_area(1) + 1
  end if
  if (Rng(icell)%tree_basal_area .lt.   0.0  ) then
    Rng(icell)%tree_basal_area = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%tree_basal_area(2) = Exceed%tree_basal_area(2) + 1
  end if

  !  soil_surface_temperature
  if (Rng(icell)%soil_surface_temperature .gt. V_LARGE) then
    Rng(icell)%soil_surface_temperature = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%soil_surface_temperature(1) = Exceed%soil_surface_temperature(1) + 1
  end if
! Soil surface temperature may occassionally be below 0.0, so commented out.
!  if (Rng(icell)%soil_surface_temperature .lt.   0.0  ) then
!    Rng(icell)%soil_surface_temperature = 0.0
!    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
!    Exceed%soil_surface_temperature(2) = Exceed%soil_surface_temperature(2) + 1
!  end if

  !  mineral_nitrogen
  do i = 1, 4
    if (Rng(icell)%mineral_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%mineral_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%mineral_nitrogen(1) = Exceed%mineral_nitrogen(1) + 1
    end if
    if (Rng(icell)%mineral_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%mineral_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%mineral_nitrogen(2) = Exceed%mineral_nitrogen(2) + 1
    end if
  end do

  !  field_capacity
  do i = 1, 4
    if (Rng(icell)%field_capacity(i) .gt. V_LARGE) then
      Rng(icell)%field_capacity(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%field_capacity(1) = Exceed%field_capacity(1) + 1
    end if
    if (Rng(icell)%field_capacity(i) .lt.   0.0  ) then
      Rng(icell)%field_capacity(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%field_capacity(2) = Exceed%field_capacity(2) + 1
    end if
  end do

  !  wilting_point
  do i = 1, 4
    if (Rng(icell)%wilting_point(i) .gt. V_LARGE) then
      Rng(icell)%wilting_point(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%wilting_point(1) = Exceed%wilting_point(1) + 1
    end if
    if (Rng(icell)%wilting_point(i) .lt.   0.0  ) then
      Rng(icell)%wilting_point(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%wilting_point(2) = Exceed%wilting_point(2) + 1
    end if
  end do

  !  soil_total_carbon
  if (Rng(icell)%soil_total_carbon .gt. V_LARGE) then
    Rng(icell)%soil_total_carbon = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%soil_total_carbon(1) = Exceed%soil_total_carbon(1) + 1
  end if
  if (Rng(icell)%soil_total_carbon .lt.   0.0  ) then
    Rng(icell)%soil_total_carbon = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%soil_total_carbon(2) = Exceed%soil_total_carbon(2) + 1
  end if

  !  tree_carbon
  do i = 1, WOODY_PARTS
    if (Rng(icell)%tree_carbon(i) .gt. V_LARGE) then
      Rng(icell)%tree_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%tree_carbon(1) = Exceed%tree_carbon(1) + 1
    end if
    if (Rng(icell)%tree_carbon(i) .lt.   0.0  ) then
      Rng(icell)%tree_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%tree_carbon(2) = Exceed%tree_carbon(2) + 1
    end if
  end do

  !  tree_nitrogen
  do i = 1, WOODY_PARTS
    if (Rng(icell)%tree_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%tree_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%tree_nitrogen(1) = Exceed%tree_nitrogen(1) + 1
    end if
    if (Rng(icell)%tree_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%tree_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%tree_nitrogen(2) = Exceed%tree_nitrogen(2) + 1
    end if
  end do

  !  shrub_carbon
  do i = 1, WOODY_PARTS
    if (Rng(icell)%shrub_carbon(i) .gt. V_LARGE) then
      Rng(icell)%shrub_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%shrub_carbon(1) = Exceed%shrub_carbon(1) + 1
    end if
    if (Rng(icell)%shrub_carbon(i) .lt.   0.0  ) then
      Rng(icell)%shrub_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%shrub_carbon(2) = Exceed%shrub_carbon(2) + 1
    end if
  end do

  !  shrub_nitrogen
  do i = 1, WOODY_PARTS
    if (Rng(icell)%shrub_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%shrub_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%shrub_nitrogen(1) = Exceed%shrub_nitrogen(1) + 1
    end if
    if (Rng(icell)%shrub_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%shrub_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%shrub_nitrogen(2) = Exceed%shrub_nitrogen(2) + 1
    end if
  end do

  !  carbon_nitrogen_ratio
  do i = 1, 2
    if (Rng(icell)%carbon_nitrogen_ratio(i) .gt. V_LARGE) then
      Rng(icell)%carbon_nitrogen_ratio(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%carbon_nitrogen_ratio(1) = Exceed%carbon_nitrogen_ratio(1) + 1
    end if
    if (Rng(icell)%carbon_nitrogen_ratio(i) .lt.   0.0  ) then
      Rng(icell)%carbon_nitrogen_ratio(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%carbon_nitrogen_ratio(2) = Exceed%carbon_nitrogen_ratio(2) + 1
    end if
  end do

  !  fast_soil_carbon
  do i = 1, 2
    if (Rng(icell)%fast_soil_carbon(i) .gt. V_LARGE) then
      Rng(icell)%fast_soil_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%fast_soil_carbon(1) = Exceed%fast_soil_carbon(1) + 1
    end if
    if (Rng(icell)%fast_soil_carbon(i) .lt.   0.0  ) then
      Rng(icell)%fast_soil_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%fast_soil_carbon(2) = Exceed%fast_soil_carbon(2) + 1
    end if
  end do

  !  intermediate_soil_carbon
  if (Rng(icell)%intermediate_soil_carbon .gt. V_LARGE) then
    Rng(icell)%intermediate_soil_carbon = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%intermediate_soil_carbon(1) = Exceed%intermediate_soil_carbon(1) + 1
  end if
  if (Rng(icell)%intermediate_soil_carbon .lt.   0.0  ) then
    Rng(icell)%intermediate_soil_carbon = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%intermediate_soil_carbon(2) = Exceed%intermediate_soil_carbon(2) + 1
  end if

  !  passive_soil_carbon
  if (Rng(icell)%passive_soil_carbon .gt. V_LARGE) then
    Rng(icell)%passive_soil_carbon = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%passive_soil_carbon(1) = Exceed%passive_soil_carbon(1) + 1
  end if
  if (Rng(icell)%passive_soil_carbon .lt.   0.0  ) then
    Rng(icell)%passive_soil_carbon = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%passive_soil_carbon(2) = Exceed%passive_soil_carbon(2) + 1
  end if

  !  fast_soil_nitrogen
  do i = 1, 2
    if (Rng(icell)%fast_soil_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%fast_soil_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%fast_soil_nitrogen(1) = Exceed%fast_soil_nitrogen(1) + 1
    end if
    if (Rng(icell)%fast_soil_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%fast_soil_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%fast_soil_nitrogen(2) = Exceed%fast_soil_nitrogen(2) + 1
    end if
  end do

  !  intermediate_soil_nitrogen
  if (Rng(icell)%intermediate_soil_nitrogen .gt. V_LARGE) then
    Rng(icell)%intermediate_soil_nitrogen = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%intermediate_soil_nitrogen(1) = Exceed%intermediate_soil_nitrogen(1) + 1
  end if
  if (Rng(icell)%intermediate_soil_nitrogen .lt.   0.0  ) then
    Rng(icell)%intermediate_soil_nitrogen = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%intermediate_soil_nitrogen(2) = Exceed%intermediate_soil_nitrogen(2) + 1
  end if

  !  passive_soil_nitrogen
  if (Rng(icell)%passive_soil_nitrogen .gt. V_LARGE) then
    Rng(icell)%passive_soil_nitrogen = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%passive_soil_nitrogen(1) = Exceed%passive_soil_nitrogen(1) + 1
  end if
  if (Rng(icell)%passive_soil_nitrogen .lt.   0.0  ) then
    Rng(icell)%passive_soil_nitrogen = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%passive_soil_nitrogen(2) = Exceed%passive_soil_nitrogen(2) + 1
  end if

  !  potential_production
  if (Rng(icell)%potential_production .gt. V_LARGE) then
    Rng(icell)%potential_production = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%potential_production(1) = Exceed%potential_production(1) + 1
  end if
  if (Rng(icell)%potential_production .lt.   0.0  ) then
    Rng(icell)%potential_production = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%potential_production(2) = Exceed%potential_production(2) + 1
  end if

  !  belowground_pot_production
  do i = 1, V_LYRS
    if (Rng(icell)%belowground_pot_production(i) .gt. V_LARGE) then
      Rng(icell)%belowground_pot_production(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%belowground_pot_production(1) = Exceed%belowground_pot_production(1) + 1
    end if
    if (Rng(icell)%belowground_pot_production(i) .lt.   0.0  ) then
      Rng(icell)%belowground_pot_production(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%belowground_pot_production(2) = Exceed%belowground_pot_production(2) + 1
    end if
  end do

  !  aboveground_pot_production
  do i = 1, V_LYRS
    if (Rng(icell)%aboveground_pot_production(i) .gt. V_LARGE) then
      Rng(icell)%aboveground_pot_production(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%aboveground_pot_production(1) = Exceed%aboveground_pot_production(1) + 1
    end if
    if (Rng(icell)%aboveground_pot_production(i) .lt.   0.0  ) then
      Rng(icell)%aboveground_pot_production(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%aboveground_pot_production(2) = Exceed%aboveground_pot_production(2) + 1
    end if
  end do

  !  total_pot_production
  do i = 1, V_LYRS
    if (Rng(icell)%total_pot_production(i) .gt. V_LARGE) then
      Rng(icell)%total_pot_production(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%total_pot_production(1) = Exceed%total_pot_production(1) + 1
    end if
    if (Rng(icell)%total_pot_production(i) .lt.   0.0  ) then
      Rng(icell)%total_pot_production(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%total_pot_production(2) = Exceed%total_pot_production(2) + 1
    end if
  end do

  do i = 1, FACETS
    !  co2_effect_on_production
    if (Rng(icell)%co2_effect_on_production(i) .gt. V_LARGE) then
      Rng(icell)%co2_effect_on_production(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%co2_effect_on_production(1) = Exceed%co2_effect_on_production(1) + 1
    end if
    if (Rng(icell)%co2_effect_on_production(i) .lt.   0.0  ) then
      Rng(icell)%co2_effect_on_production(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%co2_effect_on_production(2) = Exceed%co2_effect_on_production(2) + 1
    end if
  end do

  !  total_pot_prod_limited_by_n
  do i = 1, V_LYRS
    if (Rng(icell)%total_pot_prod_limited_by_n(i) .gt. V_LARGE) then
      Rng(icell)%total_pot_prod_limited_by_n(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%total_pot_prod_limited_by_n(1) = Exceed%total_pot_prod_limited_by_n(1) + 1
    end if
    if (Rng(icell)%total_pot_prod_limited_by_n(i) .lt.   0.0  ) then
      Rng(icell)%total_pot_prod_limited_by_n(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%total_pot_prod_limited_by_n(2) = Exceed%total_pot_prod_limited_by_n(2) + 1
    end if
  end do

  !  Monthly net primary production
  if (Rng(icell)%monthly_net_primary_production .gt. V_LARGE) then
    Rng(icell)%monthly_net_primary_production = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%monthly_net_primary_production(1) = Exceed%monthly_net_primary_production(1) + 1
  end if
  if (Rng(icell)%monthly_net_primary_production .lt.   0.0  ) then
    Rng(icell)%monthly_net_primary_production = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%monthly_net_primary_production(2) = Exceed%monthly_net_primary_production(2) + 1
  end if
  
  !  fraction_live_removed_grazing
  if (Rng(icell)%fraction_live_removed_grazing .gt. V_LARGE) then
    Rng(icell)%fraction_live_removed_grazing = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%fraction_live_removed_grazing(1) = Exceed%fraction_live_removed_grazing(1) + 1
  end if
  if (Rng(icell)%fraction_live_removed_grazing .lt.   0.0  ) then
    Rng(icell)%fraction_live_removed_grazing = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%fraction_live_removed_grazing(2) = Exceed%fraction_live_removed_grazing(2) + 1
  end if

  !  fraction_dead_removed_grazing
  if (Rng(icell)%fraction_dead_removed_grazing .gt. V_LARGE) then
    Rng(icell)%fraction_dead_removed_grazing = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%fraction_dead_removed_grazing(1) = Exceed%fraction_dead_removed_grazing(1) + 1
  end if
  if (Rng(icell)%fraction_dead_removed_grazing .lt.   0.0  ) then
    Rng(icell)%fraction_dead_removed_grazing = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%fraction_dead_removed_grazing(2) = Exceed%fraction_dead_removed_grazing(2) + 1
  end if

  !  temp_effect_on_decomp
  if (Rng(icell)%temp_effect_on_decomp .gt. V_LARGE) then
    Rng(icell)%temp_effect_on_decomp = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%temp_effect_on_decomp(1) = Exceed%temp_effect_on_decomp(1) + 1
  end if
  if (Rng(icell)%temp_effect_on_decomp .lt.   0.0  ) then
    Rng(icell)%temp_effect_on_decomp = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%temp_effect_on_decomp(2) = Exceed%temp_effect_on_decomp(2) + 1
  end if

  !  water_effect_on_decomp
  if (Rng(icell)%water_effect_on_decomp .gt. V_LARGE) then
    Rng(icell)%water_effect_on_decomp = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%water_effect_on_decomp(1) = Exceed%water_effect_on_decomp(1) + 1
  end if
  if (Rng(icell)%water_effect_on_decomp .lt.   0.0  ) then
    Rng(icell)%water_effect_on_decomp = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%water_effect_on_decomp(2) = Exceed%water_effect_on_decomp(2) + 1
  end if

  !  anerobic_effect_on_decomp
  if (Rng(icell)%anerobic_effect_on_decomp .gt. V_LARGE) then
    Rng(icell)%anerobic_effect_on_decomp = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%anerobic_effect_on_decomp(1) = Exceed%anerobic_effect_on_decomp(1) + 1
  end if
  if (Rng(icell)%anerobic_effect_on_decomp .lt.   0.0  ) then
    Rng(icell)%anerobic_effect_on_decomp = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%anerobic_effect_on_decomp(2) = Exceed%anerobic_effect_on_decomp(2) + 1
  end if

  !  all_effects_on_decomp
  if (Rng(icell)%all_effects_on_decomp .gt. V_LARGE) then
    Rng(icell)%all_effects_on_decomp = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%all_effects_on_decomp(1) = Exceed%all_effects_on_decomp(1) + 1
  end if
  if (Rng(icell)%all_effects_on_decomp .lt.   0.0  ) then
    Rng(icell)%all_effects_on_decomp = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%all_effects_on_decomp(2) = Exceed%all_effects_on_decomp(2) + 1
  end if

  !  dead_fine_root_carbon
  do i = 1, FACETS
    if (Rng(icell)%dead_fine_root_carbon(i) .gt. V_LARGE) then
      Rng(icell)%dead_fine_root_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_fine_root_carbon(1) = Exceed%dead_fine_root_carbon(1) + 1
    end if
    if (Rng(icell)%dead_fine_root_carbon(i) .lt.   0.0  ) then
      Rng(icell)%dead_fine_root_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_fine_root_carbon(2) = Exceed%dead_fine_root_carbon(2) + 1
    end if
  end do

  !  dead_fine_root_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%dead_fine_root_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%dead_fine_root_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_fine_root_nitrogen(1) = Exceed%dead_fine_root_nitrogen(1) + 1
    end if
    if (Rng(icell)%dead_fine_root_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%dead_fine_root_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_fine_root_nitrogen(2) = Exceed%dead_fine_root_nitrogen(2) + 1
    end if
  end do

  !  dead_standing_carbon
  do i = 1, FACETS
    if (Rng(icell)%dead_standing_carbon(i) .gt. V_LARGE) then
      Rng(icell)%dead_standing_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_standing_carbon(1) = Exceed%dead_standing_carbon(1) + 1
    end if
    if (Rng(icell)%dead_standing_carbon(i) .lt.   0.0  ) then
      Rng(icell)%dead_standing_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_standing_carbon(2) = Exceed%dead_standing_carbon(2) + 1
    end if
  end do

  !  dead_standing_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%dead_standing_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%dead_standing_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_standing_nitrogen(1) = Exceed%dead_standing_nitrogen(1) + 1
    end if
    if (Rng(icell)%dead_standing_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%dead_standing_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_standing_nitrogen(2) = Exceed%dead_standing_nitrogen(2) + 1
    end if
  end do

  !  dead_seed_carbon
  do i = 1, FACETS
    if (Rng(icell)%dead_seed_carbon(i) .gt. V_LARGE) then
      Rng(icell)%dead_seed_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_seed_carbon(1) = Exceed%dead_seed_carbon(1) + 1
    end if
    if (Rng(icell)%dead_seed_carbon(i) .lt.   0.0  ) then
      Rng(icell)%dead_seed_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_seed_carbon(2) = Exceed%dead_seed_carbon(2) + 1
    end if
  end do

  !  dead_seed_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%dead_seed_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%dead_seed_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_seed_nitrogen(1) = Exceed%dead_seed_nitrogen(1) + 1
    end if
    if (Rng(icell)%dead_seed_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%dead_seed_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_seed_nitrogen(2) = Exceed%dead_seed_nitrogen(2) + 1
    end if
  end do

  !  dead_leaf_carbon
  do i = 1, FACETS
    if (Rng(icell)%dead_leaf_carbon(i) .gt. V_LARGE) then
      Rng(icell)%dead_leaf_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_leaf_carbon(1) = Exceed%dead_leaf_carbon(1) + 1
    end if
    if (Rng(icell)%dead_leaf_carbon(i) .lt.   0.0  ) then
      Rng(icell)%dead_leaf_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_leaf_carbon(2) = Exceed%dead_leaf_carbon(2) + 1
    end if
  end do

  !  dead_leaf_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%dead_leaf_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%dead_leaf_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_leaf_nitrogen(1) = Exceed%dead_leaf_nitrogen(1) + 1
    end if
    if (Rng(icell)%dead_leaf_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%dead_leaf_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_leaf_nitrogen(2) = Exceed%dead_leaf_nitrogen(2) + 1
    end if
  end do

  !  dead_fine_branch_carbon
  do i = 1, FACETS
    if (Rng(icell)%dead_fine_branch_carbon(i) .gt. V_LARGE) then
      Rng(icell)%dead_fine_branch_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_fine_branch_carbon(1) = Exceed%dead_fine_branch_carbon(1) + 1
    end if
    if (Rng(icell)%dead_fine_branch_carbon(i) .lt.   0.0  ) then
      Rng(icell)%dead_fine_branch_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_fine_branch_carbon(2) = Exceed%dead_fine_branch_carbon(2) + 1
    end if
  end do

  !  dead_total_fine_branch_carbon
  if (Rng(icell)%dead_total_fine_branch_carbon .gt. V_LARGE) then
    Rng(icell)%dead_total_fine_branch_carbon = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%dead_total_fine_branch_carbon(1) = Exceed%dead_total_fine_branch_carbon(1) + 1
  end if
  if (Rng(icell)%dead_total_fine_branch_carbon .lt.   0.0  ) then
    Rng(icell)%dead_total_fine_branch_carbon = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%dead_total_fine_branch_carbon(2) = Exceed%dead_total_fine_branch_carbon(2) + 1
  end if

  !  dead_fine_branch_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%dead_fine_branch_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%dead_fine_branch_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_fine_branch_nitrogen(1) = Exceed%dead_fine_branch_nitrogen(1) + 1
    end if
    if (Rng(icell)%dead_fine_branch_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%dead_fine_branch_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_fine_branch_nitrogen(2) = Exceed%dead_fine_branch_nitrogen(2) + 1
    end if
  end do

  !  dead_total_fine_branch_nitrogen
  if (Rng(icell)%dead_total_fine_branch_nitrogen .gt. V_LARGE) then
    Rng(icell)%dead_total_fine_branch_nitrogen = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%dead_total_fine_branch_nitrogen(1) = Exceed%dead_total_fine_branch_nitrogen(1) + 1
  end if
  if (Rng(icell)%dead_total_fine_branch_nitrogen .lt.   0.0  ) then
    Rng(icell)%dead_total_fine_branch_nitrogen = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%dead_total_fine_branch_nitrogen(2) = Exceed%dead_total_fine_branch_nitrogen(2) + 1
  end if

  !  dead_coarse_root_carbon
  do i = 1, FACETS
    if (Rng(icell)%dead_coarse_root_carbon(i) .gt. V_LARGE) then
      Rng(icell)%dead_coarse_root_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_coarse_root_carbon(1) = Exceed%dead_coarse_root_carbon(1) + 1
    end if
    if (Rng(icell)%dead_coarse_root_carbon(i) .lt.   0.0  ) then
      Rng(icell)%dead_coarse_root_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_coarse_root_carbon(2) = Exceed%dead_coarse_root_carbon(2) + 1
    end if
  end do

  !  dead_total_coarse_root_carbon
  if (Rng(icell)%dead_total_coarse_root_carbon .gt. V_LARGE) then
    Rng(icell)%dead_total_coarse_root_carbon = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%dead_total_coarse_root_carbon(1) = Exceed%dead_total_coarse_root_carbon(1) + 1
  end if
  if (Rng(icell)%dead_total_coarse_root_carbon .lt.   0.0  ) then
    Rng(icell)%dead_total_coarse_root_carbon = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%dead_total_coarse_root_carbon(2) = Exceed%dead_total_coarse_root_carbon(2) + 1
  end if

  !  dead_coarse_root_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%dead_coarse_root_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%dead_coarse_root_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_coarse_root_nitrogen(1) = Exceed%dead_coarse_root_nitrogen(1) + 1
    end if
    if (Rng(icell)%dead_coarse_root_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%dead_coarse_root_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_coarse_root_nitrogen(2) = Exceed%dead_coarse_root_nitrogen(2) + 1
    end if
  end do

  !  dead_total_coarse_root_nitrogen
  if (Rng(icell)%dead_total_coarse_root_nitrogen .gt. V_LARGE) then
    Rng(icell)%dead_total_coarse_root_nitrogen = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%dead_total_coarse_root_nitrogen(1) = Exceed%dead_total_coarse_root_nitrogen(1) + 1
  end if
  if (Rng(icell)%dead_total_coarse_root_nitrogen .lt.   0.0  ) then
    Rng(icell)%dead_total_coarse_root_nitrogen = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%dead_total_coarse_root_nitrogen(2) = Exceed%dead_total_coarse_root_nitrogen(2) + 1
  end if

  !  dead_coarse_branch_carbon
  do i = 1, FACETS
    if (Rng(icell)%dead_coarse_branch_carbon(i) .gt. V_LARGE) then
      Rng(icell)%dead_coarse_branch_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_coarse_branch_carbon(1) = Exceed%dead_coarse_branch_carbon(1) + 1
    end if
    if (Rng(icell)%dead_coarse_branch_carbon(i) .lt.   0.0  ) then
      Rng(icell)%dead_coarse_branch_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_coarse_branch_carbon(2) = Exceed%dead_coarse_branch_carbon(2) + 1
    end if
  end do

  !  dead_total_coarse_branch_carbon
  if (Rng(icell)%dead_total_coarse_branch_carbon .gt. V_LARGE) then
    Rng(icell)%dead_total_coarse_branch_carbon = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%dead_total_coarse_branch_carbon(1) = Exceed%dead_total_coarse_branch_carbon(1) + 1
  end if
  if (Rng(icell)%dead_total_coarse_branch_carbon .lt.   0.0  ) then
    Rng(icell)%dead_total_coarse_branch_carbon = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%dead_total_coarse_branch_carbon(2) = Exceed%dead_total_coarse_branch_carbon(2) + 1
  end if

  !  dead_coarse_branch_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%dead_coarse_branch_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%dead_coarse_branch_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%dead_coarse_branch_nitrogen(1) = Exceed%dead_coarse_branch_nitrogen(1) + 1
    end if
    if (Rng(icell)%dead_coarse_branch_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%dead_coarse_branch_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%dead_coarse_branch_nitrogen(2) = Exceed%dead_coarse_branch_nitrogen(2) + 1
    end if
  end do

  !  dead_total_coarse_branch_nitrogen
  if (Rng(icell)%dead_total_coarse_branch_nitrogen .gt. V_LARGE) then
    Rng(icell)%dead_total_coarse_branch_nitrogen = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%dead_total_coarse_branch_nitrogen(1) = Exceed%dead_total_coarse_branch_nitrogen(1) + 1
  end if
  if (Rng(icell)%dead_total_coarse_branch_nitrogen .lt.   0.0  ) then
    Rng(icell)%dead_total_coarse_branch_nitrogen = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%dead_total_coarse_branch_nitrogen(2) = Exceed%dead_total_coarse_branch_nitrogen(2) + 1
  end if

  !  lignin_fine_root
  do i = 1, FACETS
    if (Rng(icell)%lignin_fine_root(i) .gt. V_LARGE) then
      Rng(icell)%lignin_fine_root(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%lignin_fine_root(1) = Exceed%lignin_fine_root(1) + 1
    end if
    if (Rng(icell)%lignin_fine_root(i) .lt.   0.0  ) then
      Rng(icell)%lignin_fine_root(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%lignin_fine_root(2) = Exceed%lignin_fine_root(2) + 1
    end if
  end do

  !  lignin_coarse_root
  do i = 1, FACETS
    if (Rng(icell)%lignin_coarse_root(i) .gt. V_LARGE) then
      Rng(icell)%lignin_coarse_root(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%lignin_coarse_root(1) = Exceed%lignin_coarse_root(1) + 1
    end if
    if (Rng(icell)%lignin_coarse_root(i) .lt.   0.0  ) then
      Rng(icell)%lignin_coarse_root(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%lignin_coarse_root(2) = Exceed%lignin_coarse_root(2) + 1
    end if
  end do

  !  lignin_fine_branch
  do i = 1, FACETS
    if (Rng(icell)%lignin_fine_branch(i) .gt. V_LARGE) then
      Rng(icell)%lignin_fine_branch(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%lignin_fine_branch(1) = Exceed%lignin_fine_branch(1) + 1
    end if
    if (Rng(icell)%lignin_fine_branch(i) .lt.   0.0  ) then
      Rng(icell)%lignin_fine_branch(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%lignin_fine_branch(2) = Exceed%lignin_fine_branch(2) + 1
    end if
  end do

  !  lignin_coarse_branch
  do i = 1, FACETS
    if (Rng(icell)%lignin_coarse_branch(i) .gt. V_LARGE) then
      Rng(icell)%lignin_coarse_branch(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%lignin_coarse_branch(1) = Exceed%lignin_coarse_branch(1) + 1
    end if
    if (Rng(icell)%lignin_coarse_branch(i) .lt.   0.0  ) then
      Rng(icell)%lignin_coarse_branch(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%lignin_coarse_branch(2) = Exceed%lignin_coarse_branch(2) + 1
    end if
  end do

  !  lignin_leaf
  do i = 1, FACETS
    if (Rng(icell)%lignin_leaf(i) .gt. V_LARGE) then
      Rng(icell)%lignin_leaf(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%lignin_leaf(1) = Exceed%lignin_leaf(1) + 1
    end if
    if (Rng(icell)%lignin_leaf(i) .lt.   0.0  ) then
      Rng(icell)%lignin_leaf(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%lignin_leaf(2) = Exceed%lignin_leaf(2) + 1
    end if
  end do

  !  plant_lignin_fraction
  do i = 1, FACETS
    do j = 1, 2
      if (Rng(icell)%plant_lignin_fraction(i,j) .gt. V_LARGE) then
        Rng(icell)%plant_lignin_fraction(i,j) = V_LARGE
        Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
        Exceed%plant_lignin_fraction(1) = Exceed%plant_lignin_fraction(1) + 1
      end if
      if (Rng(icell)%plant_lignin_fraction(i,j) .lt.   0.0  ) then
        Rng(icell)%plant_lignin_fraction(i,j) = 0.0
        Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
        Exceed%plant_lignin_fraction(2) = Exceed%plant_lignin_fraction(2) + 1
      end if
    end do
  end do

  !  litter_structural_carbon
  do i = 1, 2
    if (Rng(icell)%litter_structural_carbon(i) .gt. V_LARGE) then
      Rng(icell)%litter_structural_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%litter_structural_carbon(1) = Exceed%litter_structural_carbon(1) + 1
    end if
    if (Rng(icell)%litter_structural_carbon(i) .lt.   0.0  ) then
      Rng(icell)%litter_structural_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%litter_structural_carbon(2) = Exceed%litter_structural_carbon(2) + 1
    end if
  end do

  !  litter_metabolic_carbon
  do i = 1, 2
    if (Rng(icell)%litter_metabolic_carbon(i) .gt. V_LARGE) then
      Rng(icell)%litter_metabolic_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%litter_metabolic_carbon(1) = Exceed%litter_metabolic_carbon(1) + 1
    end if
    if (Rng(icell)%litter_metabolic_carbon(i) .lt.   0.0  ) then
      Rng(icell)%litter_metabolic_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%litter_metabolic_carbon(2) = Exceed%litter_metabolic_carbon(2) + 1
    end if
  end do

  !  litter_structural_nitrogen
  do i = 1, 2
    if (Rng(icell)%litter_structural_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%litter_structural_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%litter_structural_nitrogen(1) = Exceed%litter_structural_nitrogen(1) + 1
    end if
    if (Rng(icell)%litter_structural_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%litter_structural_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%litter_structural_nitrogen(2) = Exceed%litter_structural_nitrogen(2) + 1
    end if
  end do

  !  litter_metabolic_nitrogen
  do i = 1, 2
    if (Rng(icell)%litter_metabolic_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%litter_metabolic_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%litter_metabolic_nitrogen(1) = Exceed%litter_metabolic_nitrogen(1) + 1
    end if
    if (Rng(icell)%litter_metabolic_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%litter_metabolic_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%litter_metabolic_nitrogen(2) = Exceed%litter_metabolic_nitrogen(2) + 1
    end if
  end do

  !  tnetmin
  do i = 1, 2
    if (Rng(icell)%tnetmin(i) .gt. V_LARGE) then
      Rng(icell)%tnetmin(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%tnetmin(1) = Exceed%tnetmin(1) + 1
    end if
    if (Rng(icell)%tnetmin(i) .lt.   0.0  ) then
      Rng(icell)%tnetmin(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%tnetmin(2) = Exceed%tnetmin(2) + 1
    end if
  end do

  !  tminup
  do i = 1, 2
    if (Rng(icell)%tminup(i) .gt. V_LARGE) then
      Rng(icell)%tminup(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%tminup(1) = Exceed%tminup(1) + 1
    end if
    if (Rng(icell)%tminup(i) .lt.   0.0  ) then
      Rng(icell)%tminup(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%tminup(2) = Exceed%tminup(2) + 1
    end if
  end do

  !  grossmin
  do i = 1, 2
    if (Rng(icell)%grossmin(i) .gt. V_LARGE) then
      Rng(icell)%grossmin(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%grossmin(1) = Exceed%grossmin(1) + 1
    end if
    if (Rng(icell)%grossmin(i) .lt.   0.0  ) then
      Rng(icell)%grossmin(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%grossmin(2) = Exceed%grossmin(2) + 1
    end if
  end do

  !  volitn
  do i = 1, 2
    if (Rng(icell)%volitn(i) .gt. V_LARGE) then
      Rng(icell)%volitn(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%volitn(1) = Exceed%volitn(1) + 1
    end if
    if (Rng(icell)%volitn(i) .lt.   0.0  ) then
      Rng(icell)%volitn(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%volitn(2) = Exceed%volitn(2) + 1
    end if
  end do

  !  fixnit
  if (Rng(icell)%fixnit .gt. V_LARGE) then
    Rng(icell)%fixnit = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%fixnit(1) = Exceed%fixnit(1) + 1
  end if
  if (Rng(icell)%fixnit .lt.   0.0  ) then
    Rng(icell)%fixnit = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%fixnit(2) = Exceed%fixnit(2) + 1
  end if

  !  runoffn
  if (Rng(icell)%runoffn .gt. V_LARGE) then
    Rng(icell)%runoffn = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%runoffn(1) = Exceed%runoffn(1) + 1
  end if
  if (Rng(icell)%runoffn .lt.   0.0  ) then
    Rng(icell)%runoffn = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%runoffn(2) = Exceed%runoffn(2) + 1
  end if

  !  e_up
  do i = 1, FACETS
    do j = 1, WOODY_PARTS
      if (Rng(icell)%e_up(i,j) .gt. V_LARGE) then
        Rng(icell)%e_up(i,j) = V_LARGE
        Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
        Exceed%e_up(1) = Exceed%e_up(1) + 1
      end if
      if (Rng(icell)%e_up(i,j) .lt.   0.0  ) then
        Rng(icell)%e_up(i,j) = 0.0
        Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
        Exceed%e_up(2) = Exceed%e_up(2) + 1
      end if
    end do
  end do

  !  volatized_n
  if (Rng(icell)%volatized_n .gt. V_LARGE) then
    Rng(icell)%volatized_n = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%volatized_n(1) = Exceed%volatized_n(1) + 1
  end if
  if (Rng(icell)%volatized_n .lt.   0.0  ) then
    Rng(icell)%volatized_n = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%volatized_n(2) = Exceed%volatized_n(2) + 1
  end if

  !  maintain_respiration
  do i = 1, FACETS
    if (Rng(icell)%maintain_respiration(i) .gt. V_LARGE) then
      Rng(icell)%maintain_respiration(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%maintain_respiration(1) = Exceed%maintain_respiration(1) + 1
    end if
    if (Rng(icell)%maintain_respiration(i) .lt.   0.0  ) then
      Rng(icell)%maintain_respiration(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%maintain_respiration(2) = Exceed%maintain_respiration(2) + 1
    end if
  end do

  !  phenology
  do i = 1, FACETS
    if (Rng(icell)%phenology(i) .gt. V_LARGE) then
      Rng(icell)%phenology(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%phenology(1) = Exceed%phenology(1) + 1
    end if
    if (Rng(icell)%phenology(i) .lt.   0.0  ) then
      Rng(icell)%phenology(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%phenology(2) = Exceed%phenology(2) + 1
    end if
  end do

  !  fine_root_carbon
  do i = 1, FACETS
    if (Rng(icell)%fine_root_carbon(i) .gt. V_LARGE) then
      Rng(icell)%fine_root_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%fine_root_carbon(1) = Exceed%fine_root_carbon(1) + 1
    end if
    if (Rng(icell)%fine_root_carbon(i) .lt.   0.0  ) then
      Rng(icell)%fine_root_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%fine_root_carbon(2) = Exceed%fine_root_carbon(2) + 1
    end if
  end do

  !  fine_root_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%fine_root_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%fine_root_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%fine_root_nitrogen(1) = Exceed%fine_root_nitrogen(1) + 1
    end if
    if (Rng(icell)%fine_root_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%fine_root_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%fine_root_nitrogen(2) = Exceed%fine_root_nitrogen(2) + 1
    end if
  end do

  !  seed_carbon
  do i = 1, FACETS
    if (Rng(icell)%seed_carbon(i) .gt. V_LARGE) then
      Rng(icell)%seed_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%seed_carbon(1) = Exceed%seed_carbon(1) + 1
    end if
    if (Rng(icell)%seed_carbon(i) .lt.   0.0  ) then
      Rng(icell)%seed_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%seed_carbon(2) = Exceed%seed_carbon(2) + 1
    end if
  end do

  !  seed_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%seed_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%seed_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%seed_nitrogen(1) = Exceed%seed_nitrogen(1) + 1
    end if
    if (Rng(icell)%seed_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%seed_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%seed_nitrogen(2) = Exceed%seed_nitrogen(2) + 1
    end if
  end do

  !  leaf_carbon
  do i = 1, FACETS
    if (Rng(icell)%leaf_carbon(i) .gt. V_LARGE) then
      Rng(icell)%leaf_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%leaf_carbon(1) = Exceed%leaf_carbon(1) + 1
    end if
    if (Rng(icell)%leaf_carbon(i) .lt.   0.0  ) then
      Rng(icell)%leaf_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%leaf_carbon(2) = Exceed%leaf_carbon(2) + 1
    end if
  end do

  !  leaf_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%leaf_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%leaf_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%leaf_nitrogen(1) = Exceed%leaf_nitrogen(1) + 1
    end if
    if (Rng(icell)%leaf_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%leaf_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%leaf_nitrogen(2) = Exceed%leaf_nitrogen(2) + 1
    end if
  end do

  !  fine_branch_carbon
  do i = 1, FACETS
    if (Rng(icell)%fine_branch_carbon(i) .gt. V_LARGE) then
      Rng(icell)%fine_branch_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%fine_branch_carbon(1) = Exceed%fine_branch_carbon(1) + 1
    end if
    if (Rng(icell)%fine_branch_carbon(i) .lt.   0.0  ) then
      Rng(icell)%fine_branch_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%fine_branch_carbon(2) = Exceed%fine_branch_carbon(2) + 1
    end if
  end do

  !  fine_branch_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%fine_branch_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%fine_branch_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%fine_branch_nitrogen(1) = Exceed%fine_branch_nitrogen(1) + 1
    end if
    if (Rng(icell)%fine_branch_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%fine_branch_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%fine_branch_nitrogen(2) = Exceed%fine_branch_nitrogen(2) + 1
    end if
  end do

  !  coarse_root_carbon
  do i = 1, FACETS
    if (Rng(icell)%coarse_root_carbon(i) .gt. V_LARGE) then
      Rng(icell)%coarse_root_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%coarse_root_carbon(1) = Exceed%coarse_root_carbon(1) + 1
    end if
    if (Rng(icell)%coarse_root_carbon(i) .lt.   0.0  ) then
      Rng(icell)%coarse_root_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%coarse_root_carbon(2) = Exceed%coarse_root_carbon(2) + 1
    end if
  end do

  !  coarse_root_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%coarse_root_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%coarse_root_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%coarse_root_nitrogen(1) = Exceed%coarse_root_nitrogen(1) + 1
    end if
    if (Rng(icell)%coarse_root_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%coarse_root_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%coarse_root_nitrogen(2) = Exceed%coarse_root_nitrogen(2) + 1
    end if
  end do

  !  coarse_branch_carbon
  do i = 1, FACETS
    if (Rng(icell)%coarse_branch_carbon(i) .gt. V_LARGE) then
      Rng(icell)%coarse_branch_carbon(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%coarse_branch_carbon(1) = Exceed%coarse_branch_carbon(1) + 1
    end if
    if (Rng(icell)%coarse_branch_carbon(i) .lt.   0.0  ) then
      Rng(icell)%coarse_branch_carbon(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%coarse_branch_carbon(2) = Exceed%coarse_branch_carbon(2) + 1
    end if
  end do

  !  coarse_branch_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%coarse_branch_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%coarse_branch_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%coarse_branch_nitrogen(1) = Exceed%coarse_branch_nitrogen(1) + 1
    end if
    if (Rng(icell)%coarse_branch_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%coarse_branch_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%coarse_branch_nitrogen(2) = Exceed%coarse_branch_nitrogen(2) + 1
    end if
  end do

  !  stored_nitrogen
  do i = 1, FACETS
    if (Rng(icell)%stored_nitrogen(i) .gt. V_LARGE) then
      Rng(icell)%stored_nitrogen(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%stored_nitrogen(1) = Exceed%stored_nitrogen(1) + 1
    end if
    if (Rng(icell)%stored_nitrogen(i) .lt.   0.0  ) then
      Rng(icell)%stored_nitrogen(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%stored_nitrogen(2) = Exceed%stored_nitrogen(2) + 1
    end if
  end do

  !  plant_nitrogen_fixed
  do i = 1, FACETS
    if (Rng(icell)%plant_nitrogen_fixed(i) .gt. V_LARGE) then
      Rng(icell)%plant_nitrogen_fixed(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%plant_nitrogen_fixed(1) = Exceed%plant_nitrogen_fixed(1) + 1
    end if
    if (Rng(icell)%plant_nitrogen_fixed(i) .lt.   0.0  ) then
      Rng(icell)%plant_nitrogen_fixed(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%plant_nitrogen_fixed(2) = Exceed%plant_nitrogen_fixed(2) + 1
    end if
  end do

  !  nitrogen_fixed
  do i = 1, FACETS
    if (Rng(icell)%nitrogen_fixed(i) .gt. V_LARGE) then
      Rng(icell)%nitrogen_fixed(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%nitrogen_fixed(1) = Exceed%nitrogen_fixed(1) + 1
    end if
    if (Rng(icell)%nitrogen_fixed(i) .lt.   0.0  ) then
      Rng(icell)%nitrogen_fixed(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%nitrogen_fixed(2) = Exceed%nitrogen_fixed(2) + 1
    end if
  end do

  !  respiration_flows
  do i = 1, FACETS
    if (Rng(icell)%respiration_flows(i) .gt. V_LARGE) then
      Rng(icell)%respiration_flows(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%respiration_flows(1) = Exceed%respiration_flows(1) + 1
    end if
    if (Rng(icell)%respiration_flows(i) .lt.   0.0  ) then
      Rng(icell)%respiration_flows(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%respiration_flows(2) = Exceed%respiration_flows(2) + 1
    end if
  end do

  !  respiration_annual
  do i = 1, FACETS
    if (Rng(icell)%respiration_annual(i) .gt. V_LARGE) then
      Rng(icell)%respiration_annual(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%respiration_annual(1) = Exceed%respiration_annual(1) + 1
    end if
    if (Rng(icell)%respiration_annual(i) .lt.   0.0  ) then
      Rng(icell)%respiration_annual(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%respiration_annual(2) = Exceed%respiration_annual(2) + 1
    end if
  end do

  !  carbon_source_sink
  if (Rng(icell)%carbon_source_sink .gt. V_LARGE) then
    Rng(icell)%carbon_source_sink = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%carbon_source_sink(1) = Exceed%carbon_source_sink(1) + 1
  end if
  if (Rng(icell)%carbon_source_sink .lt.   0.0  ) then
    Rng(icell)%carbon_source_sink = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%carbon_source_sink(2) = Exceed%carbon_source_sink(2) + 1
  end if

  !  nitrogen_source_sink
  if (Rng(icell)%nitrogen_source_sink .gt. V_LARGE) then
    Rng(icell)%nitrogen_source_sink = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%nitrogen_source_sink(1) = Exceed%nitrogen_source_sink(1) + 1
  end if
  if (Rng(icell)%nitrogen_source_sink .lt.   0.0  ) then
    Rng(icell)%nitrogen_source_sink = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%nitrogen_source_sink(2) = Exceed%nitrogen_source_sink(2) + 1
  end if

  !  carbon_allocation
  do i = 1, FACETS
    do j = 1, WOODY_PARTS
      if (Rng(icell)%carbon_allocation(i,j) .gt. V_LARGE) then
        Rng(icell)%carbon_allocation(i,j) = V_LARGE
        Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
        Exceed%carbon_allocation(1) = Exceed%carbon_allocation(1) + 1
      end if
      if (Rng(icell)%carbon_allocation(i,j) .lt.   0.0  ) then
        Rng(icell)%carbon_allocation(i,j) = 0.0
        Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
        Exceed%carbon_allocation(2) = Exceed%carbon_allocation(2) + 1
      end if
    end do
  end do

  !  optimum_leaf_area_index
  do i = 1, FACETS
    if (Rng(icell)%optimum_leaf_area_index(i) .gt. V_LARGE) then
      Rng(icell)%optimum_leaf_area_index(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%optimum_leaf_area_index(1) = Exceed%optimum_leaf_area_index(1) + 1
    end if
    if (Rng(icell)%optimum_leaf_area_index(i) .lt.   0.0  ) then
      Rng(icell)%optimum_leaf_area_index(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%optimum_leaf_area_index(2) = Exceed%optimum_leaf_area_index(2) + 1
    end if
  end do

  !  leaf_area_index
  do i = 1, FACETS
    if (Rng(icell)%leaf_area_index(i) .gt. V_LARGE) then
      Rng(icell)%leaf_area_index(i) = V_LARGE
      Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
      Exceed%leaf_area_index(1) = Exceed%leaf_area_index(1) + 1
    end if
    if (Rng(icell)%leaf_area_index(i) .lt.   0.0  ) then
      Rng(icell)%leaf_area_index(i) = 0.0
      Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
      Exceed%leaf_area_index(2) = Exceed%leaf_area_index(2) + 1
    end if
  end do

  !  water_function
  if (Rng(icell)%water_function .gt. V_LARGE) then
    Rng(icell)%water_function = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%water_function(1) = Exceed%water_function(1) + 1
  end if
  if (Rng(icell)%water_function .lt.   0.0  ) then
    Rng(icell)%water_function = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%water_function(2) = Exceed%water_function(2) + 1
  end if


  !  fire_severity
  if (Rng(icell)%fire_severity .gt. V_LARGE) then
    Rng(icell)%fire_severity = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%fire_severity(1) = Exceed%fire_severity(1) + 1
  end if
  if (Rng(icell)%fire_severity .lt.   0.0  ) then
    Rng(icell)%fire_severity = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%fire_severity(2) = Exceed%fire_severity(2) + 1
  end if

  ! burned_carbon
  if (Rng(icell)%burned_carbon .gt. V_LARGE) then
    Rng(icell)%burned_carbon = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%burned_carbon(1) = Exceed%burned_carbon(1) + 1
  end if
  if (Rng(icell)%burned_carbon .lt.   0.0  ) then
    Rng(icell)%burned_carbon = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%burned_carbon(2) = Exceed%burned_carbon(2) + 1
  end if

  !  burned_nitrogen
  if (Rng(icell)%burned_nitrogen .gt. V_LARGE) then
    Rng(icell)%burned_nitrogen = V_LARGE
    Rng(icell)%large_error_count = Rng(icell)%large_error_count + 1
    Exceed%burned_nitrogen(1) = Exceed%burned_nitrogen(1) + 1
  end if
  if (Rng(icell)%burned_nitrogen .lt.   0.0  ) then
    Rng(icell)%burned_nitrogen = 0.0
    Rng(icell)%neg_error_count = Rng(icell)%neg_error_count + 1
    Exceed%burned_nitrogen(2) = Exceed%burned_nitrogen(2) + 1
  end if

  iunit = Rng(icell)%range_type

  live_carbon = Rng(icell)%leaf_carbon(H_FACET) + Rng(icell)%leaf_carbon(S_FACET) + Rng(icell)%leaf_carbon(T_FACET) + &
                 Rng(icell)%fine_branch_carbon(S_FACET) + Rng(icell)%fine_branch_carbon(T_FACET)
  dead_carbon = Rng(icell)%dead_standing_carbon(H_FACET) + Rng(icell)%dead_standing_carbon(S_FACET) + &
                 Rng(icell)%dead_standing_carbon(T_FACET)
  if ( (live_carbon + dead_carbon) .gt. 0.0) then
     Rng(icell)%fraction_live_removed_grazing = ( Parms(iunit)%fraction_grazed / 12.0 ) * &
                                                ( live_carbon / ( live_carbon + dead_carbon ) )
     Rng(icell)%fraction_dead_removed_grazing = ( Parms(iunit)%fraction_grazed / 12.0 ) * &
                                                ( 1.0 - ( live_carbon / ( live_carbon + dead_carbon ) ) )
  else
     Rng(icell)%fraction_live_removed_grazing = 0.0
     Rng(icell)%fraction_dead_removed_grazing = 0.0
  end if

 end do
 
 if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'EACH_MN')

end subroutine


subroutine Zero_Accumulators
!**** Zero-out some accumulators each month.
!****
!**** R. Boone    Last modified: May 31, 2011
 use Parameter_Vars
 use Structures
 implicit none
 integer icell, ifacet

 do icell = 1, range_cells
   !Zeroing out the dead plant stores each month, as in Savanna (ZFLOW.F)
   !These have all been partitioned by now, either into forms of litter or passed to standing dead, which is not zeroed-out.
   do ifacet = 1, FACETS
     Rng(icell)%dead_fine_root_carbon(ifacet) = 0.0
     Rng(icell)%dead_fine_root_nitrogen(ifacet) = 0.0
     Rng(icell)%dead_fine_branch_carbon(ifacet) = 0.0
     Rng(icell)%dead_fine_branch_nitrogen(ifacet) = 0.0
     Rng(icell)%dead_seed_carbon(ifacet) = 0.0
     Rng(icell)%dead_seed_nitrogen(ifacet) = 0.0
     Rng(icell)%dead_leaf_carbon(ifacet) = 0.0
     Rng(icell)%dead_leaf_nitrogen(ifacet) = 0.0
     Rng(icell)%dead_coarse_branch_carbon(ifacet) = 0.0
     Rng(icell)%dead_coarse_branch_nitrogen(ifacet) = 0.0
     Rng(icell)%dead_coarse_root_carbon(ifacet) = 0.0
     Rng(icell)%dead_coarse_root_nitrogen(ifacet) = 0.0
   end do
 end do

 if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'ZERO_AC')

end subroutine


subroutine Read_Other
!**** Read in the maps that are associated with fire or fertilization, or other management.
!****
!**** R. Boone    Last modified: May 31, 2011
   use Parameter_Vars
   use Structures
   implicit none
   integer ix, iy
   integer in_year, in_month
   character(180) big_string
   character(20) in_string

   ! If fire is being modeled by maps ...
   if (Sim_Parm%fire_maps_used .ne. 0) then
     ! Open the fire history parameter file.  This file is stored in PARMS, regardless of where surfaces are.
     ! Determine the appropriate map to be using for the current month and year.
     ! This could be read more efficiently, but it should never be too lengthy, so a future modification
     ! NOTE: At a broader scale, this is also inefficient, in that the same map could be read month after month, but again, no matter right now.
     open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//'/Parms/'//Sim_Parm%fire_file_name, ACTION='READ', IOSTAT=ioerr)
     if (ioerr == 0) then
       big_string = ''
       in_year = 1
       do while (in_year .gt. 0)
         read(SHORT_USE_FILE,*) in_year, in_month, in_string
         if (in_year .le. year .and. in_year .gt. 0) then
           if ((in_year .lt. year .and. in_year .gt. 0) .or. (in_year .eq. year .and. in_month .le. month)) then
             in_string = trim(adjustl(in_string))
             write(big_string,*) trim(Sim_Parm%fire_path), trim(in_string(1:len_trim(in_string)))
           end if
         end if
       end do
     else
       write(*,*) 'Cannot read the fire history parameter file: ', &
                   app_path(1:len_trim(app_path))//'/Parms/'//Sim_Parm%fire_file_name
       stop
     end if

     ! Read in the fire map indicated
     if (Sim_Parm%echo_level .gt. 0) write(*,*) 'Reading in the fire map: ',adjustl(big_string)
     call Read_Map ('fire    ', big_string)
     ! Transfer the map
     do iy=1,y_dim
       do ix=1,x_dim
         if (Globe(ix,iy)%temporary .lt. 0) then
           Globe(ix,iy)%prop_burned = 0
         else
           Globe(ix,iy)%prop_burned = Globe(ix,iy)%temporary
         end if
       end do
     end do
     ! These maps are not echoed, due to their large size.  Specific maps may be printed as output, rather than echoed, as a check of their form.
   end if

   ! If fertilization is being modeled by maps ...
   if (Sim_Parm%fertilize_maps_used .ne. 0) then
     ! Open the fire history parameter file.  This file is stored in PARMS, regardless of where surfaces are.
     ! Determine the appropriate map to be using for the current month and year.
     ! This could be read more efficiently, but it should never be too lengthy, so a future modification
     ! NOTE: At a broader scale, this is also inefficient, in that the same map could be read month after month, but again, no matter right now.
     open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//'/Parms/'//Sim_Parm%fertilization_file_name, &
                          ACTION='READ', IOSTAT=ioerr)
     if (ioerr == 0) then
       big_string = ''
       in_year = 1
       do while (in_year .gt. 0)
         read(SHORT_USE_FILE,*) in_year, in_month, in_string
         if (in_year .le. year .and. in_year .gt. 0) then
           if ((in_year .lt. year .and. in_year .gt. 0) .or. (in_year .eq. year .and. in_month .le. month)) then
             in_string = trim(adjustl(in_string))
             write(big_string,*) trim(Sim_Parm%fertilization_path), trim(in_string(1:len_trim(in_string)))
           end if
         end if
       end do
     else
       write(*,*) 'Cannot read the fertilization history parameter file: ', &
                   app_path(1:len_trim(app_path))//'/Parms/'//Sim_Parm%fertilization_file_name
       stop
     end if

     ! Read in the fertilization map indicated
     if (Sim_Parm%echo_level .gt. 0) write(*,*) 'Reading in the fertilization map: ',adjustl(big_string)
     call Read_Map ('fert    ', big_string)
     ! Transfer the map
     do iy=1,y_dim
       do ix=1,x_dim
         if (Globe(ix,iy)%temporary .lt. 0) then
           Globe(ix,iy)%prop_fertilized = 0
         else
           Globe(ix,iy)%prop_fertilized = Globe(ix,iy)%temporary
         end if
       end do
     end do
     ! These maps are not echoed, due to their large size.  Specific maps may be printed as output, rather than echoed, as a check of their form.
   end if

end subroutine


subroutine Management (icell)
!**** Read in the maps that are associated with fire or fertilization, or other management.
!****
!**** R. Boone    Last modified: May 31, 2011
   use Parameter_Vars
   use Structures
   implicit none
   integer icell, iunit
   real proportion_cell_fertilized, c_added, n_added, organic_fert_lignin
   real harvest(1)

  iunit = Rng(icell)%range_type

  ! Simulate fertilization.  This will almost certainly be set through maps, but perhaps deposition or the like may be simulated
  ! using the probability option, so I will leave it.  It is there because of the parallels with fire.
  if (Sim_Parm%fertilize_maps_used .eq. .TRUE. .or. &
   (Sim_Parm%fertilize_maps_used .eq. .FALSE. .and. Parms(iunit)%frequency_of_fertilization .gt. 0.0 &
     .and. Parms(iunit)%fraction_fertilized .gt. 0.0)) then
    if (Sim_Parm%fertilize_maps_used .eq. .TRUE.) then
      ! Model fire, with their occurrence determined in maps.  The maps will store the proportion of each cell burned.
      ! No month is included here.  If someone wants to give detailed month-level fire maps, they may
      proportion_cell_fertilized = Globe(Rng(icell)%x, Rng(icell)%y)%prop_fertilized
    else
      ! Fertilized based on probabilies and percentages
      ! Fertilized is confined to one month, otherwise the method would be overly complex, requiring checks to judge which months are appropriate.
      if (Parms(iunit)%fertilize_month .eq. month) then
        ! The cell may burn
        call random_number(harvest)       ! Ensure that the seed is being set.
        if (Parms(iunit)%frequency_of_fertilization .gt. harvest(1)) then
          ! Some portion of the cell will burn ...
          proportion_cell_fertilized = Parms(iunit)%fraction_fertilized
        end if
      end if
    end if

    ! If some of the cell is to be fertilized, do that
    if (proportion_cell_fertilized .gt. 0.0) then
      Rng(icell)%mineral_nitrogen(SURFACE_INDEX) = Rng(icell)%mineral_nitrogen(SURFACE_INDEX) + &
         ( Parms(iunit)%fertilize_nitrogen_added * Parms(iunit)%fraction_fertilized)
      Rng(icell)%fertilized_nitrogen_added = Rng(icell)%fertilized_nitrogen_added + &
         ( Parms(iunit)%fertilize_nitrogen_added * Parms(iunit)%fraction_fertilized)
      ! Add organic matter if requested.  It gets partitioned in litter, following CENTURY
      ! NOTE that the lignin ratio is hardwired here.
      c_added = Parms(iunit)%fertilize_carbon_added * Parms(iunit)%fraction_fertilized
      n_added = Parms(iunit)%fertilize_nitrogen_added * Parms(iunit)%fraction_fertilized
      organic_fert_lignin = 0.20
      call Partition_Litter (icell, SOIL_INDEX, c_added, n_added, organic_fert_lignin, 'fert')
      Rng(icell)%fertilized_nitrogen_added = Rng(icell)%fertilized_nitrogen_added + n_added
      Rng(icell)%fertilized_carbon_added = Rng(icell)%fertilized_carbon_added + c_added
    end if
! else ... no fertilization
  end if

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'MANAGED')

end subroutine


subroutine check_for_nan (icell, calling)
!**** Check the contents of the entire Rng structure (for a single cell) for -NaN.  This should reduce the risk
!**** of -NaN appearing invisibly and propagating through several places in Rng before it becomes evident.   
!**** The error will still need to be tracked, but it should help to have a routine to focus upon.
!****
!**** R. Boone    Last modified: July 3, 2014
  use Parameter_Vars
  use Structures
  implicit none

  integer  icell, ifacet, ilyr, ipart
  character*7  calling


  ! X                                                                                                                                                                                
  if ( Rng(icell)%x .ne. Rng(icell)%x ) then                                                                                                                                         
    write(*,*) 'A -NaN error has occurred in rangeland variable x for cell ', icell, ' called by ', calling                                                                          
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable x due to -NaN.'                                                                                            
      Rng(icell)%x = 0.0                                                                                                                                                             
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Y                                                                                                                                                                                
  if ( Rng(icell)%y .ne. Rng(icell)%y ) then                                                                                                                                         
    write(*,*) 'A -NaN error has occurred in rangeland variable y for cell ', icell, ' called by ', calling                                                                          
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable y due to -NaN.'                                                                                            
      Rng(icell)%y = 0.0                                                                                                                                                             
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Range_Type                                                                                                                                                                       
  if ( Rng(icell)%range_type .ne. Rng(icell)%range_type ) then                                                                                                                       
    write(*,*) 'A -NaN error has occurred in rangeland variable range_type  for cell ', icell, ' called by ', calling                                                                
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable range_type due to -NaN.'                                                                                   
      Rng(icell)%range_type = 0.0                                                                                                                                                    
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Last_Month_Day_Length                                                                                                                                                            
  if ( Rng(icell)%last_month_day_length .ne. Rng(icell)%last_month_day_length ) then                                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable last_month_day_length for cell ', icell, ' called by ', calling                                                      
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable last_month_day_length due to -NaN.'                                                                        
      Rng(icell)%last_month_day_length = 0.0                                                                                                                                         
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Day_Length_Increasing                                                                                                                                                            
  if ( Rng(icell)%day_length_increasing .ne. Rng(icell)%day_length_increasing ) then                                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable day_length_increasing for cell ', icell, ' called by ', calling                                                      
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to FALSE for cell ', icell,' for variable day_length_increasing due to -NaN.'                                                                        
      Rng(icell)%day_length_increasing = .FALSE.                                                                                                                                         
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Day_Length                                                                                                                                                                       
  if ( Rng(icell)%day_length .ne. Rng(icell)%day_length ) then                                                                                                                       
    write(*,*) 'A -NaN error has occurred in rangeland variable day_length for cell ', icell, ' called by ', calling                                                                 
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable day_length due to -NaN.'                                                                                   
      Rng(icell)%day_length = 0.0                                                                                                                                                    
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Heat_Accumulation                                                                                                                                                                
  if ( Rng(icell)%heat_accumulation .ne. Rng(icell)%heat_accumulation ) then                                                                                                         
    write(*,*) 'A -NaN error has occurred in rangeland variable heat_accumulation for cell ', icell, ' called by ', calling                                                          
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable heat_accumulation due to -NaN.'                                                                            
      Rng(icell)%heat_accumulation = 0.0                                                                                                                                             
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Facet_Cover                                                                                                                                                                      
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%facet_cover(ifacet) .ne. Rng(icell)%facet_cover(ifacet) ) then                                                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable facet_cover for ', ifacet, ' for cell ', icell, ' called by ', &
                 calling                                             
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable facet_cover for ', ifacet, ' due to -NaN.'                                                               
        Rng(icell)%facet_cover(ifacet) = 0.0                                                                                                                                         
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Total_Population                                                                                                                                                                 
  do ilyr = 1, V_LYRS                                                                                                                                                                
    if ( Rng(icell)%total_population(ilyr) .ne. Rng(icell)%total_population(ilyr) ) then                                                                                             
      write(*,*) 'A -NaN error has occurred in rangeland variable total_population for layer ', ilyr, ' for cell ', icell, &       
                 ' called by ', calling                                    
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable total_population for layer ', ilyr, ' due to -NaN.'                                                      
        Rng(icell)%total_population(ilyr) = 0.0                                                                                                                                      
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Bare_Cover                                                                                                                                                                       
  if ( Rng(icell)%bare_cover .ne. Rng(icell)%bare_cover ) then                                                                                                                       
    write(*,*) 'A -NaN error has occurred in rangeland variable bare_cover for cell ', icell, ' called by ', calling                                                                 
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable bare_cover due to -NaN.'                                                                                   
      Rng(icell)%bare_cover = 0.0                                                                                                                                                    
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Prop_Annual_Decid                                                                                                                                                                
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%prop_annual_decid(ifacet) .ne. Rng(icell)%prop_annual_decid(ifacet) ) then                                                                                       
      write(*,*) 'A -NaN error has occurred in rangeland variable prop_annual_decid for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                 
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable prop_annual_decid for facet ', ifacet, &
                   ' due to -NaN.'                                                   
        Rng(icell)%prop_annual_decid(ifacet) = 0.0                                                                                                                                   
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Pot_Evap                                                                                                                                                                         
  if ( Rng(icell)%pot_evap .ne. Rng(icell)%pot_evap ) then                                                                                                                           
    write(*,*) 'A -NaN error has occurred in rangeland variable pot_evap for cell ', icell, ' called by ', calling                                                                   
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable pot_evap due to -NaN.'                                                                                     
      Rng(icell)%pot_evap = 0.0                                                                                                                                                      
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Evaporation                                                                                                                                                                      
  if ( Rng(icell)%evaporation .ne. Rng(icell)%evaporation ) then                                                                                                                     
    write(*,*) 'A -NaN error has occurred in rangeland variable evaporation for cell ', icell, ' called by ', calling                                                                
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable evaporation due to -NaN.'                                                                                  
      Rng(icell)%evaporation = 0.0                                                                                                                                                   
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Snow                                                                                                                                                                             
  if ( Rng(icell)%snow .ne. Rng(icell)%snow ) then                                                                                                                                   
    write(*,*) 'A -NaN error has occurred in rangeland variable snow for cell ', icell, ' called by ', calling                                                                       
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable snow due to -NaN.'                                                                                         
      Rng(icell)%snow = 0.0                                                                                                                                                          
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Snow_Liquid                                                                                                                                                                      
  if ( Rng(icell)%snow_liquid .ne. Rng(icell)%snow_liquid ) then                                                                                                                     
    write(*,*) 'A -NaN error has occurred in rangeland variable snow_liquid for cell ', icell, ' called by ', calling                                                                
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable snow_liquid due to -NaN.'                                                                                  
      Rng(icell)%snow_liquid = 0.0                                                                                                                                                   
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Old_Snow                                                                                                                                                                         
  if ( Rng(icell)%old_snow .ne. Rng(icell)%old_snow ) then                                                                                                                           
    write(*,*) 'A -NaN error has occurred in rangeland variable old_snow for cell ', icell, ' called by ', calling                                                                   
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable old_snow due to -NaN.'                                                                                     
      Rng(icell)%old_snow = 0.0                                                                                                                                                      
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Old_Snow_Liquid                                                                                                                                                                  
  if ( Rng(icell)%old_snow_liquid .ne. Rng(icell)%old_snow_liquid ) then                                                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable old_snow_liquid for cell ', icell, ' called by ', calling                                                            
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable old_snow_liquid due to -NaN.'                                                                              
      Rng(icell)%old_snow_liquid = 0.0                                                                                                                                               
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Melt                                                                                                                                                                             
  if ( Rng(icell)%melt .ne. Rng(icell)%melt ) then                                                                                                                                   
    write(*,*) 'A -NaN error has occurred in rangeland variable melt for cell ', icell, ' called by ', calling                                                                       
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable melt due to -NaN.'                                                                                         
      Rng(icell)%melt = 0.0                                                                                                                                                          
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Pet_Remaining                                                                                                                                                                    
  if ( Rng(icell)%pet_remaining .ne. Rng(icell)%pet_remaining ) then                                                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable pet_remaining for cell ', icell, ' called by ', calling                                                              
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable pet_remaining due to -NaN.'                                                                                
      Rng(icell)%pet_remaining = 0.0                                                                                                                                                 
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! PPT_Soil                                                                                                                                                                         
  if ( Rng(icell)%ppt_soil .ne. Rng(icell)%ppt_soil ) then                                                                                                                           
    write(*,*) 'A -NaN error has occurred in rangeland variable ppt_soil for cell ', icell, ' called by ', calling                                                                   
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable ppt_soil due to -NaN.'                                                                                     
      Rng(icell)%ppt_soil = 0.0                                                                                                                                                      
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Runoff                                                                                                                                                                           
  if ( Rng(icell)%runoff .ne. Rng(icell)%runoff ) then                                                                                                                               
    write(*,*) 'A -NaN error has occurred in rangeland variable runoff for cell ', icell, ' called by ', calling                                                                     
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable runoff due to -NaN.'                                                                                       
      Rng(icell)%runoff = 0.0                                                                                                                                                        
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Ratio_Water_PET                                                                                                                                                                  
  if ( Rng(icell)%ratio_water_pet .ne. Rng(icell)%ratio_water_pet ) then                                                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable ratio_water_pet for cell ', icell, ' called by ', calling                                                            
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable ratio_water_pet due to -NaN.'                                                                              
      Rng(icell)%ratio_water_pet = 0.0                                                                                                                                               
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
!  ! CO2_Value                                                                                                                                                                        
!  if ( Rng(icell)%co2_value .ne. Rng(icell)%co2_value ) then                                                                                                                       
!    write(*,*) 'A -NaN error has occurred in rangeland variable co2_value for cell ', icell, ' called by ', calling                                                                  
!    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
!      stop                                                                                                                                                                           
!    else                                                                                                                                                                             
!      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable co2_value due to -NaN.'                                                                                    
!      Rng(icell)%co2_value = 0.0                                                                                                                                                     
!    end if                                                                                                                                                                           
!  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! PET_Top_Soil                                                                                                                                                                     
  if ( Rng(icell)%pet_top_soil .ne. Rng(icell)%pet_top_soil ) then                                                                                                                   
    write(*,*) 'A -NaN error has occurred in rangeland variable pet_top_soil for cell ', icell, ' called by ', calling                                                               
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable pet_top_soil due to -NaN.'                                                                                 
      Rng(icell)%pet_top_soil = 0.0                                                                                                                                                  
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! N_Leached                                                                                                                                                                        
  do ilyr = 1, SOIL_LAYERS                                                                                                                                                           
    if ( Rng(icell)%n_leached(ilyr) .ne. Rng(icell)%n_leached(ilyr) ) then                                                                                                           
      write(*,*) 'A -NaN error has occurred in rangeland variable n_leached for layer ', ilyr, ' for cell ', icell, ' called by ', &
                  calling                                           
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable n_leached for layer ', ilyr, ' due to -NaN.'                                                             
        Rng(icell)%n_leached(ilyr) = 0.0                                                                                                                                             
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
!  ! Stream                                                                                                                                                                           
!  do ilyr = 1, 8                                                                                                                                                                     
!    if ( Rng(icell)%stream(ilyr) .ne. Rng(icell)% stream(ilyr) ) then                                                                                                                
!      write(*,*) 'A -NaN error has occurred in rangeland variable stream for layer ', ilyr, '  for cell ', icell, ' called by ', &
!                 calling                                             
!      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
!        stop                                                                                                                                                                         
!      else                                                                                                                                                                           
!        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable stream for layer ', ilyr, ' due to -NaN.'                                                                
!        Rng(icell)%stream(ilyr) = 0.0                                                                                                                                                
!      end if                                                                                                                                                                         
!    end if                                                                                                                                                                           
!  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! ASMOS                                                                                                                                                                            
  do ilyr = 1, SOIL_LAYERS                                                                                                                                                           
    if ( Rng(icell)%asmos(ilyr) .ne. Rng(icell)%asmos(ilyr) ) then                                                                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable asmos for layer ', ilyr, ' for cell ', icell, ' called by ', &
                 calling                                               
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable asmos for layer ', ilyr, ' due to -NaN.'                                                                 
        Rng(icell)%asmos(ilyr) = 0.0                                                                                                                                                 
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! AMOV                                                                                                                                                                             
  do ilyr = 1, SOIL_LAYERS                                                                                                                                                           
    if ( Rng(icell)%amov(ilyr) .ne. Rng(icell)%amov(ilyr) ) then                                                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable amov for layer ', ilyr, ' for cell ', icell, ' called by ', &
                  calling                                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable amov for layer ', ilyr, ' due to -NaN.'                                                                  
        Rng(icell)%amov(ilyr) = 0.0                                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Storm_Flow                                                                                                                                                                       
  if ( Rng(icell)%storm_flow .ne. Rng(icell)%storm_flow ) then                                                                                                                       
    write(*,*) 'A -NaN error has occurred in rangeland variable storm_flow for cell ', icell, ' called by ', calling                                                                 
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable storm_flow due to -NaN.'                                                                                   
      Rng(icell)%storm_flow = 0.0                                                                                                                                                    
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Holding_Tank                                                                                                                                                                     
  if ( Rng(icell)%holding_tank .ne. Rng(icell)%holding_tank ) then                                                                                                                   
    write(*,*) 'A -NaN error has occurred in rangeland variable holding_tank for cell ', icell, ' called by ', calling                                                               
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable holding_tank due to -NaN.'                                                                                 
      Rng(icell)%holding_tank = 0.0                                                                                                                                                  
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Transpiration                                                                                                                                                                    
  if ( Rng(icell)%transpiration .ne. Rng(icell)%transpiration ) then                                                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable transpiration for cell ', icell, ' called by ', calling                                                              
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable transpiration due to -NaN.'                                                                                
      Rng(icell)%transpiration = 0.0                                                                                                                                                 
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Relative_Water_Content                                                                                                                                                           
  do ilyr = 1, SOIL_LAYERS                                                                                                                                                           
    if ( Rng(icell)%relative_water_content(ilyr) .ne. Rng(icell)%relative_water_content(ilyr) ) then                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable relative_water_content for layer ', ilyr, &
                 ' for cell ', icell, ' called by ', calling                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable relative_water_content for layer ', ilyr, &
                   ' due to -NaN.'                                                
        Rng(icell)%relative_water_content(ilyr) = 0.0                                                                                                                                
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Water_Available                                                                                                                                                                  
  do ilyr = 1, 3                                                                                                                                                                     
    if ( Rng(icell)%water_available(ilyr) .ne. Rng(icell)%water_available(ilyr) ) then                                                                                               
      write(*,*) 'A -NaN error has occurred in rangeland variable water_available for layer ', ilyr, ' for cell ', &
                 icell, ' called by ', calling                                     
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable water_available for layer ', ilyr, &
                   ' due to -NaN.'                                                       
        Rng(icell)%water_available(ilyr) = 0.0                                                                                                                                       
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Annual_Evpotranspiration                                                                                                                                                         
  if ( Rng(icell)%annual_evapotranspiration .ne. Rng(icell)%annual_evapotranspiration ) then                                                                                         
    write(*,*) 'A -NaN error has occurred in rangeland variable annual_evapotranspiration for cell ', icell, ' called by ', calling                                                  
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable annual_evapotranspiration due to -NaN.'                                                                    
      Rng(icell)%annual_evapotranspiration = 0.0                                                                                                                                     
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Total_AGround_Live_Biomass                                                                                                                                                       
  if ( Rng(icell)%total_aground_live_biomass .ne. Rng(icell)%total_aground_live_biomass ) then                                                                                       
    write(*,*) 'A -NaN error has occurred in rangeland variable total_aground_live_biomass for cell ', icell, ' called by ', &
                calling                                                 
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable total_aground_live_biomass due to -NaN.'                                                                   
      Rng(icell)%total_aground_live_biomass = 0.0                                                                                                                                    
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Total_BGround_Live_Biomass                                                                                                                                                       
  if ( Rng(icell)%total_bground_live_biomass .ne. Rng(icell)%total_bground_live_biomass ) then                                                                                       
    write(*,*) 'A -NaN error has occurred in rangeland variable total_bground_live_biomass for cell ', icell, ' called by ', &
                calling                                                 
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable total_bground_live_biomass due to -NaN.'                                                                   
      Rng(icell)%total_bground_live_biomass = 0.0                                                                                                                                    
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Total_Litter_Carbon                                                                                                                                                              
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%total_litter_carbon(ilyr) .ne. Rng(icell)%total_litter_carbon(ilyr) ) then                                                                                       
      write(*,*) 'A -NaN error has occurred in rangeland variable total_litter_carbon for layer ', ilyr, ' for cell ', icell, &
                 ' called by ', calling                                 
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable total_litter_carbon for layer ', ilyr, &
                   ' due to -NaN.'                                                   
        Rng(icell)%total_litter_carbon(ilyr) = 0.0                                                                                                                                   
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Total_Litter_Nitrogen                                                                                                                                                            
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%total_litter_nitrogen(ilyr) .ne. Rng(icell)%total_litter_nitrogen(ilyr) ) then                                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable total_litter_nitrogen for layer ', ilyr, ' for cell ', icell, &
                 ' called by ', calling                               
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable total_litter_nitrogen for layer ', ilyr, &
                   ' due to -NaN.'                                                 
        Rng(icell)%total_litter_nitrogen(ilyr) = 0.0                                                                                                                                 
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
!  ! Non-Symiotic_Soil_N_Fixing                                                                                                                                                       
!  if ( Rng(icell)%non_symbiotic_soil_n_fixing .ne. Rng(icell)%non_symbiotic_soil_n_fixing ) then                                                                                     
!    write(*,*) 'A -NaN error has occurred in rangeland variable non_symbiotic_soil_n_fixing for cell ', icell, ' called by ', &
!                calling                                                
!    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
!      stop                                                                                                                                                                           
!    else                                                                                                                                                                             
!      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable non_symbiotic_soil_n_fixing due to -NaN.'                                                                  
!      Rng(icell)%non_symbiotic_soil_n_fixing = 0.0                                                                                                                                   
!    end if                                                                                                                                                                           
!  end if                                                                                                                                                                             
!                                                                                                                                                                                     
!  ! Atmospheric_N_Fixing                                                                                                                                                             
!  if ( Rng(icell)%atmosphere_n_fixing .ne. Rng(icell)%atmosphere_n_fixing ) then                                                                                                     
!    write(*,*) 'A -NaN error has occurred in rangeland variable atmosphere_n_fixing for cell ', icell, ' called by ', calling                                                        
!    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
!      stop                                                                                                                                                                           
!    else                                                                                                                                                                             
!      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable atmosphere_n_fixing due to -NaN.'                                                                          
!      Rng(icell)%atmosphere_n_fixing = 0.0                                                                                                                                           
!    end if                                                                                                                                                                           
!  end if                                                                                                                                                                             
!                                                                                                                                                                                     
!  ! Base_N_Deposition                                                                                                                                                                
!  if ( Rng(icell)%base_n_deposition .ne. Rng(icell)%base_n_deposition ) then                                                                                                         
!    write(*,*) 'A -NaN error has occurred in rangeland variable base_n_deposition for cell ', icell, ' called by ', calling                                                          
!    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
!      stop                                                                                                                                                                           
!    else                                                                                                                                                                             
!      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable base_n_deposition due to -NaN.'                                                                            
!      Rng(icell)%base_n_deposition = 0.0                                                                                                                                             
!    end if                                                                                                                                                                           
!  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Root_Shoot_Ratio                                                                                                                                                                 
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%root_shoot_ratio(ifacet) .ne. Rng(icell)%root_shoot_ratio(ifacet) ) then                                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable root_shoot_ratio for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                  
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable root_shoot_ratio for facet ', ifacet, &
                   ' due to -NaN.'                                                    
        Rng(icell)%root_shoot_ratio(ifacet) = 0.0                                                                                                                                    
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Tree_Basal_Area                                                                                                                                                                  
  if ( Rng(icell)%tree_basal_area .ne. Rng(icell)%tree_basal_area ) then                                                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable tree_basal_area for cell ', icell, ' called by ', calling                                                            
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable tree_basal_area due to -NaN.'                                                                              
      Rng(icell)%tree_basal_area = 0.0                                                                                                                                               
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Soil_Surface_Temperature                                                                                                                                                         
  if ( Rng(icell)%soil_surface_temperature .ne. Rng(icell)%soil_surface_temperature ) then                                                                                           
    write(*,*) 'A -NaN error has occurred in rangeland variable soil_surface_temperature for cell ', icell, ' called by ', calling                                                   
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable soil_surface_temperature due to -NaN.'                                                                     
      Rng(icell)%soil_surface_temperature = 0.0                                                                                                                                      
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Sand                                                                                                                                                                             
  do ilyr = 1, 4                                                                                                                                                                     
    if ( Rng(icell)%sand(ilyr) .ne. Rng(icell)%sand(ilyr) ) then                                                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable sand for layer ', ilyr, ' for cell ', icell, ' called by ', &
                  calling                                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable sand for layer ', ilyr, ' due to -NaN.'                                                                  
        Rng(icell)%sand(ilyr) = 0.0                                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Silt                                                                                                                                                                             
  do ilyr = 1, 4                                                                                                                                                                     
    if ( Rng(icell)%silt(ilyr) .ne. Rng(icell)%silt(ilyr) ) then                                                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable silt for layer ', ilyr, ' for cell ', icell, ' called by ', &
                  calling                                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable silt for layer ', ilyr, ' due to -NaN.'                                                                  
        Rng(icell)%silt(ilyr) = 0.0                                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Clay                                                                                                                                                                             
  do ilyr = 1, 4                                                                                                                                                                     
    if ( Rng(icell)%clay(ilyr) .ne. Rng(icell)%clay(ilyr) ) then                                                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable clay for layer ', ilyr, ' for cell ', icell, ' called by ', &
                  calling                                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable clay for layer ', ilyr, ' due to -NaN.'                                                                  
        Rng(icell)%clay(ilyr) = 0.0                                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
!  ! Gravel                                                                                                                                                                           
!  do ilyr = 1, 4                                                                                                                                                                     
!    if ( Rng(icell)%gravel(ilyr) .ne. Rng(icell)%gravel(ilyr) ) then                                                                                                                 
!      write(*,*) 'A -NaN error has occurred in rangeland variable gravel for layer ', ilyr, ' for cell ', icell, ' called by ', &
!                  calling                                              
!      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
!        stop                                                                                                                                                                         
!      else                                                                                                                                                                           
!        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable gravel for layer ', ilyr, ' due to -NaN.'                                                                
!        Rng(icell)%gravel(ilyr) = 0.0                                                                                                                                                
!      end if                                                                                                                                                                         
!    end if                                                                                                                                                                           
!  end do                                                                                                                                                                             
!                                                                                                                                                                                     
!  ! Bulk_Density                                                                                                                                                                     
!  do ilyr = 1, 4                                                                                                                                                                     
!    if ( Rng(icell)%bulk_density(ilyr) .ne. Rng(icell)%bulk_density(ilyr) ) then                                                                                                     
!      write(*,*) 'A -NaN error has occurred in rangeland variable bulk_density for layer ', ilyr, ' for cell ', icell, &
!                 ' called by ', calling                                        
!      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
!        stop                                                                                                                                                                         
!      else                                                                                                                                                                           
!        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable bulk_density for layer ', ilyr, ' due to -NaN.'                                                          
!        Rng(icell)%bulk_density(ilyr) = 0.0                                                                                                                                          
!      end if                                                                                                                                                                         
!    end if                                                                                                                                                                           
!  end do                                                                                                                                                                             
!                                                                                                                                                                                     
!  ! Organic_Carbon                                                                                                                                                                   
!  do ilyr = 1, 4                                                                                                                                                                     
!    if ( Rng(icell)%organic_carbon(ilyr) .ne. Rng(icell)%organic_carbon(ilyr) ) then                                                                                                 
!      write(*,*) 'A -NaN error has occurred in rangeland variable organic_carbon for layer ', ilyr, ' for cell ', icell, &
!                 ' called by ', calling                                      
!      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
!        stop                                                                                                                                                                         
!      else                                                                                                                                                                           
!        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable organic_carbon for layer ', ilyr, ' due to -NaN.'                                                        
!        Rng(icell)%organic_carbon(ilyr) = 0.0                                                                                                                                        
!      end if                                                                                                                                                                         
!    end if                                                                                                                                                                           
!  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Mineral_Nitrogen                                                                                                                                                                 
  do ilyr = 1, 4                                                                                                                                                                     
    if ( Rng(icell)%mineral_nitrogen(ilyr) .ne. Rng(icell)%mineral_nitrogen(ilyr) ) then                                                                                             
      write(*,*) 'A -NaN error has occurred in rangeland variable mineral_nitrogen for layer ', ilyr, ' for cell ', icell, &
                 ' called by ', calling                                    
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable mineral_nitrogen for layer ', ilyr, ' due to -NaN.'                                                      
        Rng(icell)%mineral_nitrogen(ilyr) = 0.0                                                                                                                                      
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Field_Capacity                                                                                                                                                                   
  do ilyr = 1, 4                                                                                                                                                                     
    if ( Rng(icell)%field_capacity(ilyr) .ne. Rng(icell)%field_capacity(ilyr) ) then                                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable field_capacity for layer ', ilyr, ' for cell ', icell, &
                 ' called by ', calling                                      
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable field_capacity for layer ', ilyr, ' due to -NaN.'                                                        
        Rng(icell)%field_capacity(ilyr) = 0.0                                                                                                                                        
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Wilting_Point                                                                                                                                                                    
  do ilyr = 1, 4                                                                                                                                                                     
    if ( Rng(icell)%wilting_point(ilyr) .ne. Rng(icell)%wilting_point(ilyr) ) then                                                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable wilting_point for layer ', ilyr, ' for cell ', icell, &
                 ' called by ', calling                                       
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable wilting_point for layer ', ilyr, ' due to -NaN.'                                                         
        Rng(icell)%wilting_point(ilyr) = 0.0                                                                                                                                         
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Soil_Total_Carbon                                                                                                                                                                
  if ( Rng(icell)%soil_total_carbon .ne. Rng(icell)%soil_total_carbon ) then                                                                                                         
    write(*,*) 'A -NaN error has occurred in rangeland variable soil_total_carbon for cell ', icell, ' called by ', calling                                                          
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable soil_total_carbon due to -NaN.'                                                                            
      Rng(icell)%soil_total_carbon = 0.0                                                                                                                                             
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Tree_Carbon                                                                                                                                                                      
  do ipart = 1, WOODY_PARTS                                                                                                                                                          
    if ( Rng(icell)%tree_carbon(ipart) .ne. Rng(icell)%tree_carbon(ipart) ) then                                                                                                    
      write(*,*) 'A -NaN error has occurred in rangeland variable tree_carbon for part ', ipart, ' for cell ', icell, &
                 ' called by ', calling                                         
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable tree_carbon for part ', ipart, ' due to -NaN.'                                                           
        Rng(icell)%tree_carbon(ipart) = 0.0                                                                                                                                          
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Tree_Nitrogen                                                                                                                                                                    
  do ipart = 1, WOODY_PARTS                                                                                                                                                          
    if ( Rng(icell)%tree_nitrogen(ipart) .ne. Rng(icell)%tree_nitrogen(ipart) ) then                                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable tree_nitrogen for part ', ipart, ' for cell ', icell, &
                 ' called by ', calling                                       
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable tree_nitrogen for part ', ipart, ' due to -NaN.'                                                         
        Rng(icell)%tree_nitrogen(ipart) = 0.0                                                                                                                                        
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Shrub_Carbon                                                                                                                                                                     
  do ipart = 1, WOODY_PARTS                                                                                                                                                          
    if ( Rng(icell)%shrub_carbon(ipart) .ne. Rng(icell)%shrub_carbon(ipart) ) then                                                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable shrub_carbon for part ', ipart, ' for cell ', icell, &
                 ' called by ', calling                                        
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable shrub_carbon for part ', ipart, ' due to -NaN.'                                                          
        Rng(icell)%shrub_carbon(ipart) = 0.0                                                                                                                                         
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Shrub_Nitrogen                                                                                                                                                                   
  do ipart = 1, WOODY_PARTS                                                                                                                                                          
    if ( Rng(icell)%shrub_nitrogen(ipart) .ne. Rng(icell)%shrub_nitrogen(ipart) ) then                                                                                               
      write(*,*) 'A -NaN error has occurred in rangeland variable shrub_nitrogen for part ', ipart, ' for cell ', icell, &
                 ' called by ', calling                                      
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable shrub_nitrogen for part ', ipart, ' due to -NaN.'                                                        
        Rng(icell)%shrub_nitrogen(ipart) = 0.0                                                                                                                                       
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Carbon_Nitrogen_Ratio                                                                                                                                                            
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%carbon_nitrogen_ratio(ilyr) .ne. Rng(icell)%carbon_nitrogen_ratio(ilyr) ) then                                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable carbon_nitrogen_ratio for layer ', ilyr, ' for cell ', icell, &
                 ' called by ', calling                               
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable carbon_nitrogen_ratio for layer ', ilyr, &
                   ' due to -NaN.'                                                 
        Rng(icell)%carbon_nitrogen_ratio(ilyr) = 0.0                                                                                                                                 
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fast_Soil_Carbon                                                                                                                                                                 
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%fast_soil_carbon(ilyr) .ne. Rng(icell)%fast_soil_carbon(ilyr) ) then                                                                                             
      write(*,*) 'A -NaN error has occurred in rangeland variable fast_soil_carbon for layer ', ilyr, ' for cell ', icell, &
                 ' called by ', calling                                    
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fast_soil_carbon for layer ', ilyr, ' due to -NaN.'                                                      
        Rng(icell)%fast_soil_carbon(ilyr) = 0.0                                                                                                                                      
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Intermediate_Soil_Carbon                                                                                                                                                         
  if ( Rng(icell)%intermediate_soil_carbon .ne. Rng(icell)%intermediate_soil_carbon ) then                                                                                           
    write(*,*) 'A -NaN error has occurred in rangeland variable intermediate_soil_carbon for cell ', icell, ' called by ', calling                                                   
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable intermediate_soil_carbon due to -NaN.'                                                                     
      Rng(icell)%intermediate_soil_carbon = 0.0                                                                                                                                      
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Passive_Soil_Carbon                                                                                                                                                              
  if ( Rng(icell)%passive_soil_carbon .ne. Rng(icell)%passive_soil_carbon ) then                                                                                                     
    write(*,*) 'A -NaN error has occurred in rangeland variable passive_soil_carbon for cell ', icell, ' called by ', calling                                                        
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable passive_soil_carbon due to -NaN.'                                                                          
      Rng(icell)%passive_soil_carbon = 0.0                                                                                                                                           
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fast_Soil_Nitrogen                                                                                                                                                               
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%fast_soil_nitrogen(ilyr) .ne. Rng(icell)%fast_soil_nitrogen(ilyr) ) then                                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable fast_soil_nitrogen for layer ', ilyr, ' for cell ', icell, &
                 ' called by ', calling                                  
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fast_soil_nitrogen for layer ', ilyr, &
                   ' due to -NaN.'                                                    
        Rng(icell)%fast_soil_nitrogen(ilyr) = 0.0                                                                                                                                    
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Intermediate_Soil_Nitrogen                                                                                                                                                       
  if ( Rng(icell)%intermediate_soil_nitrogen .ne. Rng(icell)%intermediate_soil_nitrogen ) then                                                                                       
    write(*,*) 'A -NaN error has occurred in rangeland variable intermediate_soil_nitrogen for cell ', icell, ' called by ', &
               calling                                                 
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable intermediate_soil_nitrogen due to -NaN.'                                                                   
      Rng(icell)%intermediate_soil_nitrogen = 0.0                                                                                                                                    
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Passive_Soil_Nitrogen                                                                                                                                                            
  if ( Rng(icell)%passive_soil_nitrogen .ne. Rng(icell)%passive_soil_nitrogen ) then                                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable passive_soil_nitrogen for cell ', icell, ' called by ', calling                                                      
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable passive_soil_nitrogen due to -NaN.'                                                                        
      Rng(icell)%passive_soil_nitrogen = 0.0                                                                                                                                         
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Potential_Production                                                                                                                                                             
  if ( Rng(icell)%potential_production .ne. Rng(icell)%potential_production ) then                                                                                                   
    write(*,*) 'A -NaN error has occurred in rangeland variable potential_production for cell ', icell, ' called by ', calling                                                       
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable potential_production due to -NaN.'                                                                         
      Rng(icell)%potential_production = 0.0                                                                                                                                          
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Belowground_Pot_Production                                                                                                                                                       
  do ilyr = 1, V_LYRS                                                                                                                                                                
    if ( Rng(icell)%belowground_pot_production(ilyr) .ne. Rng(icell)%belowground_pot_production(ilyr) ) then                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable belowground_potential_production for layer ', ilyr, &
                 ' for cell ', icell, ' called by ', calling                    
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable belowground_potential_production for layer ', &
                   ilyr, ' due to -NaN.'                                      
        Rng(icell)%belowground_pot_production(ilyr) = 0.0                                                                                                                            
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Aboveground_Pot_Production                                                                                                                                                       
  do ilyr = 1, V_LYRS                                                                                                                                                                
    if ( Rng(icell)%aboveground_pot_production(ilyr) .ne. Rng(icell)%aboveground_pot_production(ilyr) ) then                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable aboveground_potential_production for layer ', ilyr, &
                 ' for cell ', icell, ' called by ', calling                    
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable aboveground_potential_production for layer ', &
                   ilyr, ' due to -NaN.'                                      
        Rng(icell)%aboveground_pot_production(ilyr) = 0.0                                                                                                                            
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do  
                                                                                                                                                                                     
  ! Total_Pot_Production                                                                                                                                                             
  do ilyr = 1, V_LYRS                                                                                                                                                                
    if ( Rng(icell)%total_pot_production(ilyr) .ne. Rng(icell)%total_pot_production(ilyr) ) then                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable total_potential_production for layer ', ilyr, ' for cell ', &
                 icell, ' called by ', calling                          
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable total_potential_production for layer ', ilyr, &
                   ' due to -NaN.'                                            
        Rng(icell)%total_pot_production(ilyr) = 0.0                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! CO2_Effect_on_Production                                                                                                                                                         
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%co2_effect_on_production(ifacet) .ne. Rng(icell)%co2_effect_on_production(ifacet) ) then                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable co2_effect_on_production for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                          
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable co2_effect_on_production for facet ', ifacet, &
                   ' due to -NaN.'                                            
        Rng(icell)%co2_effect_on_production(ifacet) = 0.0                                                                                                                            
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Total_Pot_Prod_Limited_by_N                                                                                                                                                      
  do ilyr = 1, V_LYRS                                                                                                                                                                
    if ( Rng(icell)%total_pot_prod_limited_by_n(ilyr) .ne. Rng(icell)%total_pot_prod_limited_by_n(ilyr) ) then                                                                       
      write(*,*) 'A -NaN error has occurred in rangeland variable total_potential_production_limited_by_n for layer ', &                                                             
                 ilyr, ' for cell ', icell, ' called by ', calling                                                                                                                   
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell, &
                   ' for variable total_potential_production_limited_by_n for layer ', ilyr, ' due to -NaN.'                                                                                                                                               
        Rng(icell)%total_pot_prod_limited_by_n(ilyr) = 0.0                                                                                                                           
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             

  ! Monthly net primary production
  if ( Rng(icell)%monthly_net_primary_production .ne. Rng(icell)%monthly_net_primary_production ) then                                                                       
    write(*,*) 'A -NaN error has occurred in rangeland variable monthly_net_primary_production ', &                                                             
               ' for cell ', icell, ' called by ', calling                                                                                                                   
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
      stop                                                                                                                                                                         
    else                                                                                                                                                                           
      write(*,*) 'Warning, resetting to zero for cell ', icell, &
                 ' for variable monthly_net_primary_production due to -NaN.'                                                                                                                                               
      Rng(icell)%monthly_net_primary_production = 0.0                                                                                                                           
    end if                                                                                                                                                                         
  end if                                                                                                                                                                           
                                                                                                                                                                                     
  ! Fraction_Live_Removed_Grazing                                                                                                                                                    
  if ( Rng(icell)%fraction_live_removed_grazing .ne. Rng(icell)%fraction_live_removed_grazing ) then                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable fraction_live_removed_grazing for cell ', icell, ' called by ', &
               calling                                              
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fraction_live_removed_grazing due to -NaN.'                                                                
      Rng(icell)%fraction_live_removed_grazing = 0.0                                                                                                                                 
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fraction_Dead_Removed_Grazing                                                                                                                                                    
  if ( Rng(icell)%fraction_dead_removed_grazing .ne. Rng(icell)%fraction_dead_removed_grazing ) then                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable fraction_dead_removed_grazing for cell ', icell, ' called by ', &
               calling                                              
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fraction_dead_removed_grazing due to -NaN.'                                                                
      Rng(icell)%fraction_dead_removed_grazing = 0.0                                                                                                                                 
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Temp_Effect_on_Decomp                                                                                                                                                            
  if ( Rng(icell)%temp_effect_on_decomp .ne. Rng(icell)%temp_effect_on_decomp ) then                                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable temp_effect_on_decomp for cell ', icell, ' called by ', calling                                                      
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable temp_effect_on_decomp due to -NaN.'                                                                        
      Rng(icell)%temp_effect_on_decomp = 0.0                                                                                                                                         
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Water_Effect_on_Decomp                                                                                                                                                           
  if ( Rng(icell)%water_effect_on_decomp .ne. Rng(icell)%water_effect_on_decomp ) then                                                                                               
    write(*,*) 'A -NaN error has occurred in rangeland variable water_effect_on_decomp for cell ', icell, ' called by ', calling                                                     
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable water_effect_on_decomp due to -NaN.'                                                                       
      Rng(icell)%water_effect_on_decomp = 0.0                                                                                                                                        
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Anerobic_Effect_on_Decomp                                                                                                                                                        
  if ( Rng(icell)%anerobic_effect_on_decomp .ne. Rng(icell)%anerobic_effect_on_decomp ) then                                                                                         
    write(*,*) 'A -NaN error has occurred in rangeland variable anerobic_effect_on_decomp for cell ', icell, ' called by ', calling                                                  
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable anerobic_effect_on_decomp due to -NaN.'                                                                    
      Rng(icell)%anerobic_effect_on_decomp = 0.0                                                                                                                                     
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! All_Effects_on_Decomp                                                                                                                                                            
  if ( Rng(icell)%all_effects_on_decomp .ne. Rng(icell)%all_effects_on_decomp ) then                                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable all_effects_on_decomp for cell ', icell, ' called by ', calling                                                      
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable all_effects_on_decomp due to -NaN.'                                                                        
      Rng(icell)%all_effects_on_decomp = 0.0                                                                                                                                         
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Fine_Root_Carbon                                                                                                                                                            
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_fine_root_carbon(ifacet) .ne. Rng(icell)%dead_fine_root_carbon(ifacet) ) then                                                                               
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_fine_root_carbon for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                             
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_fine_root_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                               
        Rng(icell)%dead_fine_root_carbon(ifacet) = 0.0                                                                                                                               
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Fine_Root_Nitrogen                                                                                                                                                          
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_fine_root_nitrogen(ifacet) .ne. Rng(icell)%dead_fine_root_nitrogen(ifacet) ) then                                                                           
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_fine_root_nitrogen for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                           
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_fine_root_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                             
        Rng(icell)%dead_fine_root_nitrogen(ifacet) = 0.0                                                                                                                             
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Standing_Carbon                                                                                                                                                             
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_standing_carbon(ifacet) .ne. Rng(icell)%dead_standing_carbon(ifacet) ) then                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_standing_carbon for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_standing_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                                
        Rng(icell)%dead_standing_carbon(ifacet) = 0.0                                                                                                                                
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Standing_Nitrogen                                                                                                                                                           
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_standing_nitrogen(ifacet) .ne. Rng(icell)%dead_standing_nitrogen(ifacet) ) then                                                                             
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_standing_nitrogen for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                            
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_standing_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                              
        Rng(icell)%dead_standing_nitrogen(ifacet) = 0.0                                                                                                                              
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Seed_Carbon                                                                                                                                                                 
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_seed_carbon(ifacet) .ne. Rng(icell)%dead_seed_carbon(ifacet) ) then                                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_seed_carbon for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                  
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_seed_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                                    
        Rng(icell)%dead_seed_carbon(ifacet) = 0.0                                                                                                                                    
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Seed_Nitrogen                                                                                                                                                               
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_seed_nitrogen(ifacet) .ne. Rng(icell)%dead_seed_nitrogen(ifacet) ) then                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_seed_nitrogen for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_seed_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                                  
        Rng(icell)%dead_seed_nitrogen(ifacet) = 0.0                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Leaf_Carbon                                                                                                                                                                 
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_leaf_carbon(ifacet) .ne. Rng(icell)%dead_leaf_carbon(ifacet) ) then                                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_leaf_carbon for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                  
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_leaf_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                                    
        Rng(icell)%dead_leaf_carbon(ifacet) = 0.0                                                                                                                                    
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Leaf_Nitrogen                                                                                                                                                               
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_leaf_nitrogen(ifacet) .ne. Rng(icell)%dead_leaf_nitrogen(ifacet) ) then                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_leaf_nitrogen for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_leaf_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                                  
        Rng(icell)%dead_leaf_nitrogen(ifacet) = 0.0                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Fine_Branch_Carbon                                                                                                                                                          
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_fine_branch_carbon(ifacet) .ne. Rng(icell)%dead_fine_branch_carbon(ifacet) ) then                                                                           
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_fine_branch_carbon for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                           
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_fine_branch_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                             
        Rng(icell)%dead_fine_branch_carbon(ifacet) = 0.0                                                                                                                             
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Total_Fine_Branch_Carbon                                                                                                                                                    
  if ( Rng(icell)%dead_total_fine_branch_carbon .ne. Rng(icell)%dead_total_fine_branch_carbon ) then                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable dead_total_fine_branch carbon for cell ', icell, ' called by ', &
               calling                                              
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_total_fine_branch carbon due to -NaN.'                                                                
      Rng(icell)%dead_total_fine_branch_carbon = 0.0                                                                                                                                 
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Fine_Branch_Nitrogen                                                                                                                                                        
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_fine_branch_nitrogen(ifacet) .ne. Rng(icell)%dead_fine_branch_nitrogen(ifacet) ) then                                                                       
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_fine_branch_nitrogen for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                         
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_fine_branch_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                           
        Rng(icell)%dead_fine_branch_nitrogen(ifacet) = 0.0                                                                                                                           
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Total_Fine_Branch_Nitrogen                                                                                                                                                  
  if ( Rng(icell)%dead_total_fine_branch_nitrogen .ne. Rng(icell)%dead_total_fine_branch_nitrogen ) then                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable dead_total_fine_branch_nitrogen for cell ', icell, &
               ' called by ', calling                                            
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_total_fine_branch_nitrogen due to -NaN.'                                                              
      Rng(icell)%dead_total_fine_branch_nitrogen = 0.0                                                                                                                               
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Coarse_Root_Carbon                                                                                                                                                          
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_coarse_root_carbon(ifacet) .ne. Rng(icell)%dead_coarse_root_carbon(ifacet) ) then                                                                           
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_coarse_root_carbon for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                           
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_coarse_root_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                             
        Rng(icell)%dead_coarse_root_carbon(ifacet) = 0.0                                                                                                                             
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Total_Coarse_Root_Carbon                                                                                                                                                    
  if ( Rng(icell)%dead_total_coarse_root_carbon .ne. Rng(icell)%dead_total_coarse_root_carbon ) then                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable dead_total_coarse_root_carbon for cell ', icell, &
               ' called by ', calling                                              
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_total_coarse_root_carbon due to -NaN.'                                                                
      Rng(icell)%dead_total_coarse_root_carbon = 0.0                                                                                                                                 
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Coarse_Root_Nitrogen                                                                                                                                                        
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_coarse_root_nitrogen(ifacet) .ne. Rng(icell)%dead_coarse_root_nitrogen(ifacet) ) then                                                                       
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_coarse_root_nitrogen for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                         
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_coarse_root_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                           
        Rng(icell)%dead_coarse_root_nitrogen(ifacet) = 0.0                                                                                                                           
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Total_Coarse_Root_Nitrogen                                                                                                                                                  
  if ( Rng(icell)%dead_total_coarse_root_nitrogen .ne. Rng(icell)%dead_total_coarse_root_nitrogen ) then                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable dead_total_coarse_root_nitrogen for cell ', icell, ' called by ', &
               calling                                            
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_total_coarse_root_nitrogen due to -NaN.'                                                              
      Rng(icell)%dead_total_coarse_root_nitrogen = 0.0                                                                                                                               
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Coarse_Branch_Carbon                                                                                                                                                        
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_coarse_branch_carbon(ifacet) .ne. Rng(icell)%dead_coarse_branch_carbon(ifacet) ) then                                                                       
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_coarse_branch_carbon for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                         
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_coarse_branch_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                           
        Rng(icell)%dead_coarse_branch_carbon(ifacet) = 0.0                                                                                                                           
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Total_Coarse_Branch_Carbon                                                                                                                                                  
  if ( Rng(icell)%dead_total_coarse_branch_carbon .ne. Rng(icell)%dead_total_coarse_branch_carbon ) then                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable dead_total_coarse_branch_carbon for cell ', icell, ' called by ', &
               calling                                            
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_total_coarse_branch_carbon due to -NaN.'                                                              
      Rng(icell)%dead_total_coarse_branch_carbon = 0.0                                                                                                                               
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Coarse_Branch_Nitrogen                                                                                                                                                      
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%dead_coarse_branch_nitrogen(ifacet) .ne. Rng(icell)%dead_coarse_branch_nitrogen(ifacet) ) then                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable dead_coarse_branch_nitrogen for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                       
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_coarse_branch_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                         
        Rng(icell)%dead_coarse_branch_nitrogen(ifacet) = 0.0                                                                                                                         
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Dead_Total_Coarse_Branch_Nitrogen                                                                                                                                                
  if ( Rng(icell)%dead_total_coarse_branch_nitrogen .ne. Rng(icell)%dead_total_coarse_branch_nitrogen ) then                                                                         
    write(*,*) 'A -NaN error has occurred in rangeland variable dead_total_coarse_branch_nitrogen for cell ', icell, &
               ' called by ', calling                                          
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable dead_total_coarse_branch_nitrogen due to -NaN.'                                                            
      Rng(icell)%dead_total_coarse_branch_nitrogen = 0.0                                                                                                                             
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Lignin_Fine_Root                                                                                                                                                                 
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%lignin_fine_root(ifacet) .ne. Rng(icell)%lignin_fine_root(ifacet) ) then                                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable lignin_fine_root for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                  
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable lignin_fine_root for facet ', ifacet, &
                   ' due to -NaN.'                                                    
        Rng(icell)%lignin_fine_root(ifacet) = 0.0                                                                                                                                    
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Lignin_Coarse_Root                                                                                                                                                               
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%lignin_coarse_root(ifacet) .ne. Rng(icell)%lignin_coarse_root(ifacet) ) then                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable lignin_coarse_root for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable lignin_coarse_root for facet ', ifacet, &
                   ' due to -NaN.'                                                  
        Rng(icell)%lignin_coarse_root(ifacet) = 0.0                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Lignin_Fine_Branch                                                                                                                                                               
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%lignin_fine_branch(ifacet) .ne. Rng(icell)%lignin_fine_branch(ifacet) ) then                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable lignin_fine_branch for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable lignin_fine_branch for facet ', ifacet, &
                   ' due to -NaN.'                                                  
        Rng(icell)%lignin_fine_branch(ifacet) = 0.0                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Lignin_Coarse_Branch                                                                                                                                                             
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%lignin_coarse_branch(ifacet) .ne. Rng(icell)%lignin_coarse_branch(ifacet) ) then                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable lignin_coarse_branch for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable lignin_coarse_branch for facet ', ifacet, &
                   ' due to -NaN.'                                                
        Rng(icell)%lignin_coarse_branch(ifacet) = 0.0                                                                                                                                
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Lignin_Leaf                                                                                                                                                                      
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%lignin_leaf(ifacet) .ne. Rng(icell)%lignin_leaf(ifacet) ) then                                                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable lignin_leaf for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                       
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable lignin_leaf for facet ', ifacet, ' due to -NaN.'                                                         
        Rng(icell)%lignin_leaf(ifacet) = 0.0                                                                                                                                         
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
 end do                                                                                                                                                                              
                                                                                                                                                                                     
  ! Plant_Lignin_Fraction                                                                                                                                                            
  do ifacet = 1, FACETS                                                                                                                                                              
    do ilyr = 1, SOIL_INDEX                                                                                                                                                          
      if ( Rng(icell)%plant_lignin_fraction(ifacet,ilyr) .ne. Rng(icell)%plant_lignin_fraction(ifacet,ilyr) ) then                                                                   
        write(*,*) 'A -NaN error has occurred in rangeland variable plant_lignin_fraction for facet ', ifacet, ' and layer ', &
                   ilyr, ' for cell ', icell, ' called by ', calling      
        if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                      
          stop                                                                                                                                                                       
        else                                                                                                                                                                         
          write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable plant_lignin_fraction for facet ', ifacet, &
                     ' and layer ', ilyr, ' due to -NaN.'                        
          Rng(icell)%plant_lignin_fraction(ifacet,ilyr) = 0.0                                                                                                                        
        end if                                                                                                                                                                       
      end if                                                                                                                                                                         
    end do                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Litter_Structural_Carbon                                                                                                                                                         
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%litter_structural_carbon(ilyr) .ne. Rng(icell)%litter_structural_carbon(ilyr) ) then                                                                             
      write(*,*) 'A -NaN error has occurred in rangeland variable litter_structural_carbon for layer ', ilyr, ' for cell ', &
                 icell, ' called by ', calling                            
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable litter_structural_carbon for layer ', ilyr, &
                   ' due to -NaN.'                                              
        Rng(icell)%litter_structural_carbon(ilyr) = 0.0                                                                                                                              
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do    
                                                                                                                                                                                     
  ! Litter_Metabolic_Carbon                                                                                                                                                          
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%litter_metabolic_carbon(ilyr) .ne. Rng(icell)%litter_metabolic_carbon(ilyr) ) then                                                                               
      write(*,*) 'A -NaN error has occurred in rangeland variable litter_metabolic_carbon for layer ', ilyr, ' for cell ', icell, &
                 ' called by ', calling                             
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable litter_metabolic_carbon for layer ', ilyr, &
                   ' due to -NaN.'                                               
        Rng(icell)%litter_metabolic_carbon(ilyr) = 0.0                                                                                                                               
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Litter_Structural_Nitrogen                                                                                                                                                       
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%litter_structural_nitrogen(ilyr) .ne. Rng(icell)%litter_structural_nitrogen(ilyr) ) then                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable litter_structural_nitrogen for layer ', ilyr, ' for cell ', &
                 icell, ' called by ', calling                          
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable litter_structural_nitrogen for layer ', ilyr, &
                   ' due to -NaN.'                                            
        Rng(icell)%litter_structural_nitrogen(ilyr) = 0.0                                                                                                                            
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Litter_Metabolic_Nitrogen                                                                                                                                                        
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%litter_metabolic_nitrogen(ilyr) .ne. Rng(icell)%litter_metabolic_nitrogen(ilyr) ) then                                                                           
      write(*,*) 'A -NaN error has occurred in rangeland variable litter_metabolic_nitrogen for layer ', ilyr, ' for cell ', &
                 icell, ' called by ', calling                           
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable litter_metabolic_nitrogen for layer ', ilyr, &
                   ' due to -NaN.'                                             
        Rng(icell)%litter_metabolic_nitrogen(ilyr) = 0.0                                                                                                                             
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! TNetMin                                                                                                                                                                          
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%tnetmin(ilyr) .ne. Rng(icell)%tnetmin(ilyr) ) then                                                                                                               
      write(*,*) 'A -NaN error has occurred in rangeland variable tnetmin for layer ', ilyr, ' for cell ', icell, ' called by ', &
                 calling                                             
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable tnetmin for layer ', ilyr, ' due to -NaN.'                                                               
        Rng(icell)%tnetmin(ilyr) = 0.0                                                                                                                                               
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! TMinUp                                                                                                                                                                           
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%tminup(ilyr) .ne. Rng(icell)%tminup(ilyr) ) then                                                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable tminup for layer ', ilyr, ' for cell ', icell, ' called by ', &
                 calling                                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable tminup for layer ', ilyr, ' due to -NaN.'                                                                
        Rng(icell)%tminup(ilyr) = 0.0                                                                                                                                                
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! GrossMin                                                                                                                                                                         
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%grossmin(ilyr) .ne. Rng(icell)%grossmin(ilyr) ) then                                                                                                             
      write(*,*) 'A -NaN error has occurred in rangeland variable grossmin for layer ', ilyr, ' for cell ', icell, ' called by ', &
                 calling                                            
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable grossmin for layer ', ilyr, ' due to -NaN.'                                                              
        Rng(icell)%grossmin(ilyr) = 0.0                                                                                                                                              
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! VolitN                                                                                                                                                                           
  do ilyr = 1, SOIL_INDEX                                                                                                                                                            
    if ( Rng(icell)%volitn(ilyr) .ne. Rng(icell)%volitn(ilyr) ) then                                                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable volitn for layer ', ilyr, ' for cell ', icell, ' called by ', &
                 calling                                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable volitn for layer ', ilyr, ' due to -NaN.'                                                                
        Rng(icell)%volitn(ilyr) = 0.0                                                                                                                                                
      end if                                                                                                                                                                         
    end if
  end do                                                                                                                                                                           
                                                                                                                                                                                     
  ! FixNit                                                                                                                                                                           
  if ( Rng(icell)%fixnit .ne. Rng(icell)%fixnit ) then                                                                                                                               
    write(*,*) 'A -NaN error has occurred in rangeland variable fixnit for cell ', icell, ' called by ', calling                                                                     
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fixnit due to -NaN.'                                                                                       
      Rng(icell)%fixnit = 0.0                                                                                                                                                        
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! RunoffN                                                                                                                                                                          
  if ( Rng(icell)%runoffn .ne. Rng(icell)%runoffn ) then                                                                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable runoffn for cell ', icell, ' called by ', calling                                                                    
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable runoffn due to -NaN.'                                                                                      
      Rng(icell)%runoffn = 0.0                                                                                                                                                       
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! E_Up                                                                                                                                                                             
  do ifacet = 1, FACETS                                                                                                                                                              
    do ipart = 1, WOODY_PARTS                                                                                                                                                        
      if ( Rng(icell)%e_up(ifacet,ipart) .ne. Rng(icell)%e_up(ifacet,ipart) ) then                                                                                                   
        write(*,*) 'A -NaN error has occurred in rangeland variable e_up for facet ', ifacet, ' for part ', ipart, ' for cell ', &
                   icell, ' called by ', calling                       
        if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                      
          stop                                                                                                                                                                       
        else                                                                                                                                                                         
          write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable e_up for facet ', ifacet, ' for part ', ipart, &
                     ' due to -NaN.'                                         
          Rng(icell)%e_up(ifacet,ipart) = 0.0                                                                                                                                        
        end if                                                                                                                                                                       
      end if                                                                                                                                                                         
    end do                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Volatized_N                                                                                                                                                                      
  if ( Rng(icell)%volatized_n .ne. Rng(icell)%volatized_n ) then                                                                                                                     
    write(*,*) 'A -NaN error has occurred in rangeland variable volatized_n for cell ', icell, ' called by ', calling                                                                
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable volatized_n due to -NaN.'                                                                                  
      Rng(icell)%volatized_n = 0.0                                                                                                                                                   
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Maintain_Respiration                                                                                                                                                             
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%maintain_respiration(ifacet) .ne. Rng(icell)%maintain_respiration(ifacet) ) then                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable maintain_respiration for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable maintain_respiration for facet ', ifacet, &
                   ' due to -NaN.'                                                
        Rng(icell)%maintain_respiration(ifacet) = 0.0                                                                                                                                
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Phenology                                                                                                                                                                        
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%phenology(ifacet) .ne. Rng(icell)%phenology(ifacet) ) then                                                                                                       
      write(*,*) 'A -NaN error has occurred in rangeland variable phenology for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                         
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable phenology for facet ', ifacet, ' due to -NaN.'                                                           
        Rng(icell)%phenology(ifacet) = 0.0                                                                                                                                           
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fine_Root_Carbon                                                                                                                                                                 
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%fine_root_carbon(ifacet) .ne. Rng(icell)%fine_root_carbon(ifacet) ) then                                                                                         
      write(*,*) 'A -NaN error has occurred in rangeland variable fine_root_carbon for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                  
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fine_root_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                                    
        Rng(icell)%fine_root_carbon(ifacet) = 0.0                                                                                                                                    
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fine_Root_Nitrogen                                                                                                                                                               
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%fine_root_nitrogen(ifacet) .ne. Rng(icell)%fine_root_nitrogen(ifacet) ) then                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable fine_root_nitrogen for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fine_root_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                                  
        Rng(icell)%fine_root_nitrogen(ifacet) = 0.0                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Seed_Carbon                                                                                                                                                                      
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%seed_carbon(ifacet) .ne. Rng(icell)%seed_carbon(ifacet) ) then                                                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable seed_carbon for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                       
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable seed_carbon for facet ', ifacet, ' due to -NaN.'                                                         
        Rng(icell)%seed_carbon(ifacet) = 0.0                                                                                                                                         
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Seed_Nitrogen                                                                                                                                                                    
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%seed_nitrogen(ifacet) .ne. Rng(icell)%seed_nitrogen(ifacet) ) then                                                                                               
      write(*,*) 'A -NaN error has occurred in rangeland variable seed_nitrogen for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                     
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable seed_nitrogen for facet ', ifacet, ' due to -NaN.'                                                       
        Rng(icell)%seed_nitrogen(ifacet) = 0.0                                                                                                                                       
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Leaf_Carbon                                                                                                                                                                      
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%leaf_carbon(ifacet) .ne. Rng(icell)%leaf_carbon(ifacet) ) then                                                                                                   
      write(*,*) 'A -NaN error has occurred in rangeland variable leaf_carbon for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                       
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable leaf_carbon for facet ', ifacet, ' due to -NaN.'                                                         
        Rng(icell)%leaf_carbon(ifacet) = 0.0                                                                                                                                         
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Leaf_Nitrogen                                                                                                                                                                    
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%leaf_nitrogen(ifacet) .ne. Rng(icell)%leaf_nitrogen(ifacet) ) then                                                                                               
      write(*,*) 'A -NaN error has occurred in rangeland variable leaf_nitrogen for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                     
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable leaf_nitrogen for facet ', ifacet, ' due to -NaN.'                                                       
        Rng(icell)%leaf_nitrogen(ifacet) = 0.0                                                                                                                                       
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fine_Branch_Carbon                                                                                                                                                               
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%fine_branch_carbon(ifacet) .ne. Rng(icell)%fine_branch_carbon(ifacet) ) then                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable fine_branch_carbon for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fine_branch_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                                  
        Rng(icell)%fine_branch_carbon(ifacet) = 0.0                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fine_Branch_Nitrogen                                                                                                                                                             
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%fine_branch_nitrogen(ifacet) .ne. Rng(icell)%fine_branch_nitrogen(ifacet) ) then                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable fine_branch_nitrogen for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fine_branch_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                                
        Rng(icell)%fine_branch_nitrogen(ifacet) = 0.0                                                                                                                                
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Coarse_Root_Carbon                                                                                                                                                               
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%coarse_root_carbon(ifacet) .ne. Rng(icell)%coarse_root_carbon(ifacet) ) then                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable coarse_root_carbon for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable coarse_root_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                                  
        Rng(icell)%coarse_root_carbon(ifacet) = 0.0                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Coarse_Root_Nitrogen                                                                                                                                                             
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%coarse_root_nitrogen(ifacet) .ne. Rng(icell)%coarse_root_nitrogen(ifacet) ) then                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable coarse_root_nitrogen for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable coarse_root_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                                
        Rng(icell)%coarse_root_nitrogen(ifacet) = 0.0                                                                                                                                
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Coarse_Branch_Carbon                                                                                                                                                             
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%coarse_branch_carbon(ifacet) .ne. Rng(icell)%coarse_branch_carbon(ifacet) ) then                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable coarse_branch_carbon for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable coarse_branch_carbon for facet ', ifacet, &
                   ' due to -NaN.'                                                
        Rng(icell)%coarse_branch_carbon(ifacet) = 0.0                                                                                                                                
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Coarse_Branch_Nitrogen                                                                                                                                                           
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%coarse_branch_nitrogen(ifacet) .ne. Rng(icell)%coarse_branch_nitrogen(ifacet) ) then                                                                             
      write(*,*) 'A -NaN error has occurred in rangeland variable coarse_branch_nitrogen for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                            
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable coarse_branch_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                              
        Rng(icell)%coarse_branch_nitrogen(ifacet) = 0.0                                                                                                                              
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Stored_Nitrogen                                                                                                                                                                  
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%stored_nitrogen(ifacet) .ne. Rng(icell)%stored_nitrogen(ifacet) ) then                                                                                           
      write(*,*) 'A -NaN error has occurred in rangeland variable stored_nitrogen for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                                   
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable stored_nitrogen for facet ', ifacet, &
                   ' due to -NaN.'                                                     
        Rng(icell)%stored_nitrogen(ifacet) = 0.0                                                                                                                                     
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Plant_Nitrogen_Fixed                                                                                                                                                             
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%plant_nitrogen_fixed(ifacet) .ne. Rng(icell)%plant_nitrogen_fixed(ifacet) ) then                                                                                 
      write(*,*) 'A -NaN error has occurred in rangeland variable plant_nitrogen_fixed for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                              
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable plant_nitrogen_fixed for facet ', ifacet, &
                   ' due to -NaN.'                                                
        Rng(icell)%plant_nitrogen_fixed(ifacet) = 0.0                                                                                                                                
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Nitrogen_Fixed                                                                                                                                                                   
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%nitrogen_fixed(ifacet) .ne. Rng(icell)%nitrogen_fixed(ifacet) ) then                                                                                             
      write(*,*) 'A -NaN error has occurred in rangeland variable nitrogen_fixed for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                    
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable nitrogen_fixed for facet ', ifacet, ' due to -NaN.'                                                      
        Rng(icell)%nitrogen_fixed(ifacet) = 0.0                                                                                                                                      
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Respiration_Flows                                                                                                                                                                
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%respiration_flows(ifacet) .ne. Rng(icell)%respiration_flows(ifacet) ) then                                                                                       
      write(*,*) 'A -NaN error has occurred in rangeland variable respiration_flows for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                 
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable respiration_flows for facet ', ifacet, &
                   ' due to -NaN.'                                                   
        Rng(icell)%respiration_flows(ifacet) = 0.0                                                                                                                                   
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Respiration_Annual                                                                                                                                                               
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%respiration_annual(ifacet) .ne. Rng(icell)%respiration_annual(ifacet) ) then                                                                                     
      write(*,*) 'A -NaN error has occurred in rangeland variable respiration_annual for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable respiration_annual for facet ', ifacet, &
                   ' due to -NaN.'                                                  
        Rng(icell)%respiration_annual(ifacet) = 0.0                                                                                                                                  
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Carbon_Source_Sink                                                                                                                                                               
  if ( Rng(icell)%carbon_source_sink .ne. Rng(icell)%carbon_source_sink ) then                                                                                                       
    write(*,*) 'A -NaN error has occurred in rangeland variable carbon_source_sink for cell ', icell, ' called by ', calling                                                         
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable carbon_source_sink due to -NaN.'                                                                           
      Rng(icell)%carbon_source_sink = 0.0                                                                                                                                            
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Nitrogen_Source_Sink                                                                                                                                                             
  if ( Rng(icell)%nitrogen_source_sink .ne. Rng(icell)%nitrogen_source_sink ) then                                                                                                   
    write(*,*) 'A -NaN error has occurred in rangeland variable nitrogen_source_sink for cell ', icell, ' called by ', calling                                                       
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable nitrogen_source_sink due to -NaN.'                                                                         
      Rng(icell)%nitrogen_source_sink = 0.0                                                                                                                                          
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Carbon_Allocation                                                                                                                                                                
  do ifacet = 1, FACETS                                                                                                                                                              
    do ipart = 1, WOODY_PARTS                                                                                                                                                        
      if ( Rng(icell)%carbon_allocation(ifacet,ipart) .ne. Rng(icell)%carbon_allocation(ifacet,ipart) ) then                                                                         
        write(*,*) 'A -NaN error has occurred in rangeland variable carbon_allocation for facet ', ifacet, ' for part ', ipart, &
                   ' for cell ', icell, ' called by ', calling          
        if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                      
          stop                                                                                                                                                                       
        else                                                                                                                                                                         
          write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable carbon_allocation for facet ', ifacet, &
                     ' for part ', ipart, ' due to -NaN.'                            
          Rng(icell)%carbon_allocation(ifacet,ipart) = 0.0                                                                                                                           
        end if                                                                                                                                                                       
      end if                                                                                                                                                                         
    end do                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Optimum_Leaf_Area_Index                                                                                                                                                          
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%optimum_leaf_area_index(ifacet) .ne. Rng(icell)%optimum_leaf_area_index(ifacet) ) then                                                                           
      write(*,*) 'A -NaN error has occurred in rangeland variable optimum_leaf_area_index for facet ', ifacet, ' for cell ', &
                 icell, ' called by ', calling                           
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable optimum_leaf_area_index for facet ', ifacet, &
                   ' due to -NaN.'                                             
        Rng(icell)%optimum_leaf_area_index(ifacet) = 0.0                                                                                                                             
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Leaf_Area_Index                                                                                                                                                                  
  do ifacet = 1, FACETS                                                                                                                                                              
    if ( Rng(icell)%leaf_area_index(ifacet) .ne. Rng(icell)%leaf_area_index(ifacet) ) then                                                                                           
      write(*,*) 'A -NaN error has occurred in rangeland variable leaf_area_index for facet ', ifacet, ' for cell ', icell, &
                 ' called by ', calling                                   
      if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                        
        stop                                                                                                                                                                         
      else                                                                                                                                                                           
        write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable leaf_area_index for facet ', ifacet, ' due to -NaN.'                                                     
        Rng(icell)%leaf_area_index(ifacet) = 0.0                                                                                                                                     
      end if                                                                                                                                                                         
    end if                                                                                                                                                                           
  end do                                                                                                                                                                             
                                                                                                                                                                                     
  ! Water_Function                                                                                                                                                                   
  if ( Rng(icell)%water_function .ne. Rng(icell)%water_function ) then                                                                                                               
    write(*,*) 'A -NaN error has occurred in rangeland variable water_function for cell ', icell, ' called by ', calling                                                             
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable water_function due to -NaN.'                                                                               
      Rng(icell)%water_function = 0.0                                                                                                                                                
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fire_Severity                                                                                                                                                                    
  if ( Rng(icell)%fire_severity .ne. Rng(icell)%fire_severity ) then                                                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable fire_severity for cell ', icell, ' called by ', calling                                                              
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fire_severity due to -NaN.'                                                                                
      Rng(icell)%fire_severity = 0.0                                                                                                                                                 
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Burned_Carbon                                                                                                                                                                    
  if ( Rng(icell)%burned_carbon .ne. Rng(icell)%burned_carbon ) then                                                                                                                 
    write(*,*) 'A -NaN error has occurred in rangeland variable burned_carbon for cell ', icell, ' called by ', calling                                                              
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable burned_carbon due to -NaN.'                                                                                
      Rng(icell)%burned_carbon = 0.0                                                                                                                                                 
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Burned_Nitrogen                                                                                                                                                                  
  if ( Rng(icell)%burned_nitrogen .ne. Rng(icell)%burned_nitrogen ) then                                                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable burned_nitrogen for cell ', icell, ' called by ', calling                                                            
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable burned_nitrogen due to -NaN.'                                                                              
      Rng(icell)%burned_nitrogen = 0.0                                                                                                                                               
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fertilized_Nitrogen_Added                                                                                                                                                        
  if ( Rng(icell)%fertilized_nitrogen_added .ne. Rng(icell)%fertilized_nitrogen_added ) then                                                                                         
    write(*,*) 'A -NaN error has occurred in rangeland variable fertilized_nitrogen_added for cell ', icell, ' called by ', calling                                                  
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fertlizied_nitrogen_added due to -NaN.'                                                                    
      Rng(icell)%fertilized_nitrogen_added = 0.0                                                                                                                                     
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Fertilized_Carbon_Added                                                                                                                                                          
  if ( Rng(icell)%fertilized_carbon_added .ne. Rng(icell)%fertilized_carbon_added ) then                                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable fertilized_carbon_added for cell ', icell, ' called by ', calling                                                    
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable fertilized_carbon_added due to -NaN.'                                                                      
      Rng(icell)%fertilized_carbon_added = 0.0                                                                                                                                       
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Large_Error_Count                                                                                                                                                                
  if ( Rng(icell)%large_error_count .ne. Rng(icell)%large_error_count ) then                                                                                                         
    write(*,*) 'A -NaN error has occurred in rangeland variable large_error_count for cell ', icell, ' called by ', calling                                                          
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable large_error_count due to -NaN.'                                                                            
      Rng(icell)%large_error_count = 0.0                                                                                                                                             
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             
                                                                                                                                                                                     
  ! Neg_Error_Count                                                                                                                                                                  
  if ( Rng(icell)%neg_error_count .ne. Rng(icell)%neg_error_count ) then                                                                                                             
    write(*,*) 'A -NaN error has occurred in rangeland variable neg_error_count for cell ', icell, ' called by ', calling                                                            
    if ( stop_on_nan_flag .eq. .TRUE.) then                                                                                                                                          
      stop                                                                                                                                                                           
    else                                                                                                                                                                             
      write(*,*) 'Warning, resetting to zero for cell ', icell,' for variable neg_error_count due to -NaN.'                                                                              
      Rng(icell)%neg_error_count = 0.0                                                                                                                                               
    end if                                                                                                                                                                           
  end if                                                                                                                                                                             


end subroutine