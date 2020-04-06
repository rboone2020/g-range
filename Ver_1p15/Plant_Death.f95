subroutine Plant_Part_Death (icell)
!**** Plant material dies for various reasons.  This routine includes a call to Woody_Plant_Part_Death.
!****
!**** R. Boone   Last modified: Oct 5, 2013
  use Parameter_Vars
  use Structures
  implicit none
  
  real deck5, rtdh, death_rate, death_carbon, death_nitrogen, frac_lignin
  real fuel, green, perc_green, effect_of_green, data_val(4)
  real prop_shoots_burned, prop_standing_dead_burned, prop_litter_burned, prop_burned_carbon_ash
  real prop_burned_nitrogen_ash, burned_ash_c, burned_ash_n
  real proportion_cell_burned              ! Note scale dependence
  integer icell, iunit, ifacet
  real harvest(1)                          ! Used in Fortran 95 to get a random value.
  real linear
  
  deck5 = 5.0                              ! Value appears standard in 100 files.  Not documented anywhere I can find.
  iunit = Rng(icell)%range_type

  ! Death of herb fine roots  (DROOT.F)
  if ( (Rng(icell)%fine_root_carbon(H_FACET) .gt. 0.0) .and. Rng(icell)%soil_surface_temperature .gt. 0.0) then
    ! In the following, water_available(1) is for growth, 2 is for survival, 3 for the top two layers.  
    rtdh = 1.0 - Rng(icell)%water_available(1) / (deck5 + Rng(icell)%water_available(1))                         ! Deck5 is set, so no division by 0.  Century uses the first layer ... avh20(1) like here.
    death_rate = Parms(iunit)%max_herb_root_death_rate * rtdh
    if (death_rate .gt. 0.95) then
      death_rate = 0.95
    end if
    if (death_rate .gt. 0.0) then
      death_carbon = death_rate * Rng(icell)%fine_root_carbon(H_FACET)
    else
      death_carbon = 0.0
    end if
    ! Moving away from CENTURY here, which goes into labeled materials, etc.
    if (death_rate .gt. 0.0 .and. Rng(icell)%fine_root_nitrogen(H_FACET) .gt. 0.0) then
      death_nitrogen = death_rate * Rng(icell)%fine_root_nitrogen(H_FACET) 
    else
      death_nitrogen = 0.0
    end if
    ! Do flows for plant parts
    ! Dead fine root carbon and nitrogen are short-term data holders, reset at the end of each month.  The 
    ! material is passed directly to litter.
    Rng(icell)%dead_fine_root_carbon(H_FACET) = Rng(icell)%dead_fine_root_carbon(H_FACET) + death_carbon
    Rng(icell)%fine_root_carbon(H_FACET) = Rng(icell)%fine_root_carbon(H_FACET) - death_carbon
    Rng(icell)%dead_fine_root_nitrogen(H_FACET) = Rng(icell)%dead_fine_root_nitrogen(H_FACET) + death_nitrogen
    Rng(icell)%fine_root_nitrogen(H_FACET) = Rng(icell)%fine_root_nitrogen(H_FACET) - death_nitrogen
    ! Do flows to litter, keeping track of structural and metabolic components
    ! Dead fine root carbon is already partitioned to litter in the main decomposition program.   The same is true for seeds.
    ! A portion of stored carbon for maintence respiration is loss associated with death of plant part.                 ! FLOWS may be misnamed, or the wrong indicator.    
    Rng(icell)%respiration_flows(H_FACET) = Rng(icell)%respiration_flows(H_FACET) - &
                                   ( Rng(icell)%respiration_flows(H_FACET) * death_rate ) 
  end if

  ! Death of leaves and shoots  (DSHOOT.F)   
  ! ***** INCLUDE STANDING DEAD HERE, RATHER THAN TRANSFERS TO BELOW-GROUND *****
  ! (Century uses aboveground live carbon)  
  if (Rng(icell)%leaf_carbon(H_FACET) .gt. 0.00001) then
    ! Century increases to the maximum the death rate during month of senescencee.  I am going to experiment using phenology 
    ! Century uses a series of four additional parameters for shoot death.  They describe losses due to 1) water stress, 2) phenology, and 3) shading as indicated in 4.  
    if (Rng(icell)%phenology(H_FACET) .ge. 3.95) then                          ! Comparison to 4.0 exactly may be causing an issue.
      ! Weighted so that annuals move entirely to standing dead, but are likely just a portion of the total herbs
      ! Should include death rate 
      ! Edited to include water function in death rate for non-annual plants - 02/13/2013
      death_rate = Parms(iunit)%shoot_death_rate(2) * ( 1.0 - Rng(icell)%prop_annual_decid(H_FACET) ) * &
                   Rng(icell)%water_function
    else
      death_rate = Parms(iunit)%shoot_death_rate(1) * Rng(icell)%water_function
    end if
    if (month .eq. Parms(iunit)%month_to_remove_annuals) then
      death_rate = death_rate + Rng(icell)%prop_annual_decid(H_FACET)     ! A one-time loss, so not corrected by month
    end if
    if (Rng(icell)%leaf_carbon(H_FACET) .gt. Parms(iunit)%shoot_death_rate(4)) then           ! Shoot death rate 4 stores g / m^2, a threshold for shading affecting shoot death, stored in 3
      death_rate = death_rate + Parms(iunit)%shoot_death_rate(3) 
    end if
    death_rate = min(1.0, death_rate)
    death_rate = max(0.0, death_rate)
    death_carbon = death_rate * Rng(icell)%leaf_carbon(H_FACET)
    ! Moving away from CENTURY here, which goes into labeled materials, etc.
    death_nitrogen = death_rate * Rng(icell)%leaf_nitrogen(H_FACET) 
    ! Do flows for plant parts  
    ! NOTE:  dead_leaf_carbon, dead_leaf_nitrogen are accumulators only, and cleared at the end of the month.  They
    !        are not used in modeling.  Leaves go from living to standing dead, and standing dead is the main operator.
    Rng(icell)%dead_leaf_carbon(H_FACET) = Rng(icell)%dead_leaf_carbon(H_FACET) + death_carbon
    Rng(icell)%dead_standing_carbon(H_FACET) = Rng(icell)%dead_standing_carbon(H_FACET) + death_carbon
    Rng(icell)%leaf_carbon(H_FACET) = Rng(icell)%leaf_carbon(H_FACET) - death_carbon
    Rng(icell)%dead_standing_nitrogen(H_FACET) = Rng(icell)%dead_standing_nitrogen(H_FACET) + death_nitrogen
    Rng(icell)%dead_leaf_nitrogen(H_FACET) = Rng(icell)%dead_leaf_nitrogen(H_FACET) + death_nitrogen
    Rng(icell)%leaf_nitrogen(H_FACET) = Rng(icell)%leaf_nitrogen(H_FACET) - death_nitrogen
    ! Do flows to litter, keeping track of structural and metabolic components
    ! Not using Partition_Litter here.  The leaves and shoots from herbs go to standing dead biomass.  That is handled elsewhere
    ! A portion of stored carbon for maintence respiration is loss associated with death of plant part.    ! FLOWS may be misnamed, or the wrong indicator.
    Rng(icell)%respiration_flows(H_FACET) = Rng(icell)%respiration_flows(H_FACET) - &
                                          ( Rng(icell)%respiration_flows(H_FACET) * death_rate ) 
    ! Death of seeds
    if (Rng(icell)%seed_carbon(H_FACET) .gt. 0.00001) then
      death_carbon = Rng(icell)%seed_carbon(H_FACET) * ( Parms(iunit)%fraction_seeds_not_germinated(H_FACET) / 12.0 ) 
    else
      death_carbon = 0.0
    end if
    if (Rng(icell)%seed_nitrogen(H_FACET) .gt. 0.00001) then
      death_nitrogen = Rng(icell)%seed_nitrogen(H_FACET) * ( Parms(iunit)%fraction_seeds_not_germinated(H_FACET) / 12.0 ) 
    else
      death_nitrogen = 0.0
    end if
    Rng(icell)%dead_seed_carbon(H_FACET) = Rng(icell)%dead_seed_carbon(H_FACET) + death_carbon 
    Rng(icell)%seed_carbon(H_FACET) = Rng(icell)%seed_carbon(H_FACET) - death_carbon
    Rng(icell)%dead_seed_nitrogen(H_FACET) = Rng(icell)%dead_seed_nitrogen(H_FACET) + death_nitrogen
    Rng(icell)%seed_nitrogen(H_FACET) = Rng(icell)%seed_nitrogen(H_FACET) - death_nitrogen
    ! Seeds won't play a role in respiration.  These seeds are destined for decomposition only.
    ! Do flows to litter for seeds, keeping track of structural and metabolic components
    ! Seeds are partitioned to litter in the main decomposition module, so removed here.
  end if          

  ! Leaf death, Fine branch death, Coarse stem death, Fine root death, Coarse root death
  call Woody_Plant_Part_Death (icell)   

  do ifacet = 1, FACETS
    ! Simulate fall of standing dead to litter, for all facets
    if (Rng(icell)%dead_standing_carbon(ifacet) .gt. 0.00001) then
      death_carbon = Rng(icell)%dead_standing_carbon(ifacet) * Parms(iunit)%fall_rate_of_standing_dead(ifacet)
    else
      death_carbon = 0.0
    end if
    if (Rng(icell)%dead_standing_nitrogen(ifacet) .gt. 0.00001) then
      death_nitrogen = Rng(icell)%dead_standing_nitrogen(ifacet) * Parms(iunit)%fall_rate_of_standing_dead(ifacet)
    else
      death_nitrogen = 0.0
    end if
    ! Do flows for plant parts
    Rng(icell)%dead_standing_carbon(ifacet) = Rng(icell)%dead_standing_carbon(ifacet) - death_carbon
    Rng(icell)%dead_standing_nitrogen(ifacet) = Rng(icell)%dead_standing_nitrogen(ifacet) - death_nitrogen
    ! Do flows to litter, keeping track of structural and metabolic components
    call Partition_Litter (icell, SURFACE_INDEX, death_carbon, death_nitrogen, &
        Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX), 'sdead')
  end do

  ! Simulate fire for all the facets
  ! Decide whether fire is being modeled.   If maps are being used, then move ahead (too many possibilities to judge if fire is not occurring in the maps).
  ! If maps are not being used, then check to see that the frequency of fire and the fraction burned are both greater than 0.  
  ! Otherwise, there is no fire.
  if (Sim_Parm%fire_maps_used .ne. 0 .or. &
   (Sim_Parm%fire_maps_used .eq. 0 .and. Parms(iunit)%frequency_of_fire .gt. 0.0 &
     .and. Parms(iunit)%fraction_burned .gt. 0.0)) then
    if (Sim_Parm%fire_maps_used .ne. 0) then
      ! Model fire, with their occurrence determined in maps.  The maps will store the proportion of each cell burned.
      ! No month is included here.  If someone wants to give detailed month-level fire maps, they may
      proportion_cell_burned = Globe(Rng(icell)%x, Rng(icell)%y)%prop_burned
    else
      ! Fire based on probabilies and percentages
      ! Fire is confined to one month, otherwise the method would be overly complex, requiring checks to judge which months are appropriate.
      if (Parms(iunit)%burn_month .eq. month) then
        ! The cell may burn
        call random_number(harvest)       ! Ensure that the seed is being set.
        if (Parms(iunit)%frequency_of_fire .gt. harvest(1)) then
          ! Some portion of the cell will burn ...
          proportion_cell_burned = Parms(iunit)%fraction_burned
        end if
      end if
    end if

    ! If some of the cell is to burn, do that
    if (proportion_cell_burned .gt. 0.0009) then
      ! Calculate the intensity of the fire, using the method in SAVANNA
      fuel = 0.0
      green = 0.0
      do ifacet = 1, FACETS
        fuel = fuel + Rng(icell)%leaf_carbon(ifacet) + Rng(icell)%dead_standing_carbon(ifacet) + &
               Rng(icell)%fine_branch_carbon(ifacet) + Rng(icell)%coarse_branch_carbon(ifacet)
        green = green + Rng(icell)%leaf_carbon(ifacet)

        perc_green = green / ( fuel + 0.000001 )
        data_val(1) = Parms(iunit)%green_vs_intensity(1, 1)
        data_val(2) = Parms(iunit)%green_vs_intensity(1, 2)
        data_val(3) = Parms(iunit)%green_vs_intensity(2, 1)
        data_val(4) = Parms(iunit)%green_vs_intensity(2, 2)
        effect_of_green = linear(perc_green, data_val, 2)
        data_val(1) = Parms(iunit)%fuel_vs_intensity(1)
        data_val(2) = 0.0 
        data_val(3) = Parms(iunit)%fuel_vs_intensity(2)
        data_val(4) = 1.0 
        Rng(icell)%fire_severity = linear(fuel, data_val, 2)
        ! Now calculate the proportion burned of different plant parts, and proportion ash
        ! data_val(2) = 0.0 and data_val(4) = 1.0 still, and in all that follows
          data_val(1) = Parms(iunit)%fraction_shoots_burned(ifacet, 1)
          data_val(3) = Parms(iunit)%fraction_shoots_burned(ifacet, 2)
        ! Proportion of cell burned is captured here, in these proportion burned entries that are used below.
        prop_shoots_burned = linear(Rng(icell)%fire_severity, data_val, 2) * proportion_cell_burned
          data_val(1) = Parms(iunit)%fraction_standing_dead_burned(ifacet, 1)
          data_val(3) = Parms(iunit)%fraction_standing_dead_burned(ifacet, 2)
        prop_standing_dead_burned = linear(Rng(icell)%fire_severity, data_val, 2) * proportion_cell_burned
          data_val(1) = Parms(iunit)%fraction_litter_burned(ifacet, 1)
          data_val(3) = Parms(iunit)%fraction_litter_burned(ifacet, 2)
        prop_litter_burned = linear(Rng(icell)%fire_severity, data_val, 2) * proportion_cell_burned
        prop_burned_carbon_ash = Parms(iunit)%fraction_burned_carbon_as_ash
        prop_burned_nitrogen_ash = Parms(iunit)%fraction_burned_nitrogen_as_ash
     
        ! Burn the materials
        ! LEAVES
        Rng(icell)%leaf_carbon(ifacet) = Rng(icell)%leaf_carbon(ifacet) - ( Rng(icell)%leaf_carbon(ifacet) * prop_shoots_burned )
        Rng(icell)%burned_carbon = Rng(icell)%burned_carbon + ( Rng(icell)%leaf_carbon(ifacet) * prop_shoots_burned )
        burned_ash_c = ( Rng(icell)%leaf_carbon(ifacet) * prop_shoots_burned ) * prop_burned_carbon_ash
        Rng(icell)%leaf_nitrogen(ifacet) = Rng(icell)%leaf_nitrogen(ifacet) - &
                                            ( Rng(icell)%leaf_nitrogen(ifacet) * prop_shoots_burned )
        Rng(icell)%burned_nitrogen = Rng(icell)%burned_nitrogen + ( Rng(icell)%leaf_nitrogen(ifacet) * prop_shoots_burned )
        burned_ash_n = ( Rng(icell)%leaf_nitrogen(ifacet) * prop_shoots_burned ) * prop_burned_nitrogen_ash
        ! Do flows to litter, keeping track of structural and metabolic components
        ! Assume the fraction lignin does not change with combustion ... EDIT AS NEEDED
        frac_lignin = Rng(icell)%lignin_leaf(ifacet) / Rng(icell)%leaf_carbon(ifacet)
        frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
        frac_lignin = min(0.50, frac_lignin)     
        call Partition_Litter (icell, SURFACE_INDEX, burned_ash_c, burned_ash_n, frac_lignin, 'fire_leaves')  
        ! A portion of stored carbon for maintence respiration is loss associated with death of plant part.            
        ! This is done only once, otherwise there would be multiple occurrences of the transfer
        Rng(icell)%respiration_flows(ifacet) = Rng(icell)%respiration_flows(ifacet) - &
                                               ( Rng(icell)%respiration_flows(ifacet) * prop_shoots_burned ) 
        ! SEEDS
        Rng(icell)%seed_carbon(ifacet) = Rng(icell)%seed_carbon(ifacet) - ( Rng(icell)%seed_carbon(ifacet) * prop_shoots_burned )
        Rng(icell)%burned_carbon = Rng(icell)%burned_carbon + ( Rng(icell)%seed_carbon(ifacet) * prop_shoots_burned )
        burned_ash_c = ( Rng(icell)%seed_carbon(ifacet) * prop_shoots_burned ) * prop_burned_carbon_ash
        Rng(icell)%seed_nitrogen(ifacet) = Rng(icell)%seed_nitrogen(ifacet) - &
                                           ( Rng(icell)%seed_nitrogen(ifacet) * prop_shoots_burned )
        Rng(icell)%burned_nitrogen = Rng(icell)%burned_nitrogen + ( Rng(icell)%seed_nitrogen(ifacet) * prop_shoots_burned )
        burned_ash_n = ( Rng(icell)%seed_nitrogen(ifacet) * prop_shoots_burned ) * prop_burned_nitrogen_ash
        ! Do flows to litter, keeping track of structural and metabolic components
        ! Assume the fraction lignin does not change with combustion ... EDIT AS NEEDED
        ! Using leaf lignin for seed lignin, as an approximate
        frac_lignin = Rng(icell)%lignin_leaf(ifacet) / Rng(icell)%leaf_carbon(ifacet)
        frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
        frac_lignin = min(0.50, frac_lignin)     
        call Partition_Litter (icell, SURFACE_INDEX, burned_ash_c, burned_ash_n, frac_lignin,'fire_seeds')  
        ! FINE BRANCHES
        Rng(icell)%fine_branch_carbon(ifacet) = Rng(icell)%fine_branch_carbon(ifacet) - &
                                                ( Rng(icell)%fine_branch_carbon(ifacet) * prop_shoots_burned )
        Rng(icell)%burned_carbon = Rng(icell)%burned_carbon + &
                                   ( Rng(icell)%fine_branch_carbon(ifacet) * prop_shoots_burned )
        burned_ash_c = ( Rng(icell)%fine_branch_carbon(ifacet) * prop_shoots_burned ) * prop_burned_carbon_ash
        Rng(icell)%fine_branch_nitrogen(ifacet) = Rng(icell)%fine_branch_nitrogen(ifacet) - &
                                                  ( Rng(icell)%fine_branch_nitrogen(ifacet) * prop_shoots_burned )
        Rng(icell)%burned_nitrogen = Rng(icell)%burned_nitrogen + &
                                     ( Rng(icell)%fine_branch_nitrogen(ifacet) * prop_shoots_burned )
        burned_ash_n = ( Rng(icell)%fine_branch_nitrogen(ifacet) * prop_shoots_burned ) * prop_burned_nitrogen_ash
        ! Do flows to litter, keeping track of structural and metabolic components.  BURNED MATERIALS are not going to standing dead
        ! Assume the fraction lignin does not change with combustion ... EDIT AS NEEDED
        frac_lignin = Rng(icell)%lignin_fine_branch(ifacet) / Rng(icell)%fine_branch_carbon(ifacet)
        frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
        frac_lignin = min(0.50, frac_lignin)        
        call Partition_Litter (icell, SURFACE_INDEX, burned_ash_c, burned_ash_n, frac_lignin, 'fire_fire')  
        ! COARSE BRANCHES
        Rng(icell)%coarse_branch_carbon(ifacet) = Rng(icell)%coarse_branch_carbon(ifacet) - &
                                                  ( Rng(icell)%coarse_branch_carbon(ifacet) * prop_shoots_burned )
        Rng(icell)%burned_carbon = Rng(icell)%burned_carbon + &
                                   ( Rng(icell)%coarse_branch_carbon(ifacet) * prop_shoots_burned )
        burned_ash_c = ( Rng(icell)%coarse_branch_carbon(ifacet) * prop_shoots_burned ) * prop_burned_carbon_ash
        Rng(icell)%coarse_branch_nitrogen(ifacet) = Rng(icell)%coarse_branch_nitrogen(ifacet) - &
                                                    ( Rng(icell)%coarse_branch_nitrogen(ifacet) * prop_shoots_burned )
        Rng(icell)%burned_nitrogen = Rng(icell)%burned_nitrogen + &
                                     ( Rng(icell)%coarse_branch_nitrogen(ifacet) * prop_shoots_burned )
        burned_ash_n = ( Rng(icell)%coarse_branch_nitrogen(ifacet) * prop_shoots_burned ) * prop_burned_nitrogen_ash
        ! Do flows to litter, keeping track of structural and metabolic components.  BURNED MATERIALS are not going to standing dead
        ! Assume the fraction lignin does not change with combustion ... EDIT AS NEEDED
        frac_lignin = Rng(icell)%lignin_coarse_branch(ifacet) / Rng(icell)%coarse_branch_carbon(ifacet)
        frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
        frac_lignin = min(0.50, frac_lignin)     
        call Partition_Litter (icell, SURFACE_INDEX, burned_ash_c, burned_ash_n, frac_lignin, 'fire2_fire2')  
        ! STANDING DEAD
        Rng(icell)%dead_standing_carbon(ifacet) = Rng(icell)%dead_standing_carbon(ifacet) - &
                                                  ( Rng(icell)%dead_standing_carbon(ifacet) * prop_standing_dead_burned )
        Rng(icell)%burned_carbon = Rng(icell)%burned_carbon + &
                                   ( Rng(icell)%dead_standing_carbon(ifacet) * prop_standing_dead_burned )
        burned_ash_c = ( Rng(icell)%dead_standing_carbon(ifacet) * prop_standing_dead_burned ) * prop_burned_carbon_ash
        Rng(icell)%dead_standing_nitrogen(ifacet) = Rng(icell)%dead_standing_nitrogen(ifacet) - &
                                                    ( Rng(icell)%dead_standing_nitrogen(ifacet) * prop_standing_dead_burned )
        Rng(icell)%burned_nitrogen = Rng(icell)%burned_nitrogen + &
                                     ( Rng(icell)%dead_standing_nitrogen(ifacet) * prop_standing_dead_burned )
        burned_ash_n = ( Rng(icell)%dead_standing_nitrogen(ifacet) * prop_standing_dead_burned ) * prop_burned_nitrogen_ash
        ! Do flows to litter, keeping track of structural and metabolic components.  BURNED MATERIALS are not going to standing dead
        ! Assume the fraction lignin does not change with combustion ... EDIT AS NEEDED
        ! Using leaf lignin for standing dead lignin, as an approximate
        frac_lignin = Rng(icell)%lignin_leaf(ifacet) / Rng(icell)%leaf_carbon(ifacet)
        frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
        frac_lignin = min(0.50, frac_lignin)        
        call Partition_Litter (icell, SURFACE_INDEX, burned_ash_c, burned_ash_n, frac_lignin, 'fire3_fire3')  
      end do        ! End of facet loop
      ! LITTER - STRUCTURAL CARBON
      Rng(icell)%litter_structural_carbon(SURFACE_INDEX) = Rng(icell)%litter_structural_carbon(SURFACE_INDEX) - &
                                                    ( Rng(icell)%litter_structural_carbon(SURFACE_INDEX) * prop_litter_burned )
      Rng(icell)%burned_carbon = Rng(icell)%burned_carbon + &
                                     ( Rng(icell)%litter_structural_carbon(SURFACE_INDEX) * prop_litter_burned )
      burned_ash_c = ( Rng(icell)%litter_structural_carbon(SURFACE_INDEX) * prop_litter_burned ) * prop_burned_carbon_ash
      Rng(icell)%litter_structural_nitrogen(SURFACE_INDEX) = Rng(icell)%litter_structural_nitrogen(SURFACE_INDEX) - &
                                     ( Rng(icell)%litter_structural_nitrogen(SURFACE_INDEX) * prop_litter_burned )
      Rng(icell)%burned_nitrogen = Rng(icell)%burned_nitrogen + &
                                  ( Rng(icell)%litter_structural_nitrogen(SURFACE_INDEX) * prop_litter_burned )
      burned_ash_n = ( Rng(icell)%litter_structural_nitrogen(SURFACE_INDEX) * prop_litter_burned ) * prop_burned_nitrogen_ash
      ! Flows to litter not required, given litter is burning
      ! LITTER - METABOLIC CARBON
      Rng(icell)%litter_metabolic_carbon(SURFACE_INDEX) = Rng(icell)%litter_metabolic_carbon(SURFACE_INDEX) - &
                                   ( Rng(icell)%litter_metabolic_carbon(SURFACE_INDEX) * prop_litter_burned )
      Rng(icell)%burned_carbon = Rng(icell)%burned_carbon + &
                                   ( Rng(icell)%litter_metabolic_carbon(SURFACE_INDEX) * prop_litter_burned )
      burned_ash_c = ( Rng(icell)%litter_metabolic_carbon(SURFACE_INDEX) + prop_litter_burned ) * prop_burned_carbon_ash
      Rng(icell)%litter_metabolic_nitrogen(SURFACE_INDEX) = Rng(icell)%litter_metabolic_nitrogen(SURFACE_INDEX) - &
                                   ( Rng(icell)%litter_metabolic_nitrogen(SURFACE_INDEX) * prop_litter_burned )
      Rng(icell)%burned_nitrogen = Rng(icell)%burned_nitrogen + &
                                   ( Rng(icell)%litter_metabolic_nitrogen(SURFACE_INDEX) * prop_litter_burned )
      burned_ash_n = ( Rng(icell)%litter_structural_nitrogen(SURFACE_INDEX) * prop_litter_burned ) * &
                                   prop_burned_nitrogen_ash
      ! Flows to litter not required, given litter is burning    
    end if
! else ... no fire
  end if

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'P_P_DTH')

end subroutine


subroutine Woody_Plant_Part_Death (icell)
!**** Woody plant material dies for various reasons.
!****
!**** R. Boone   Last modified: Feb 13, 2013
  use Parameter_Vars
  use Structures
  implicit none
  
  real temp_average, death_carbon, c_prop, death_nitrogen, to_storage
  integer  icell, iunit, ifacet
  logical  in_winter, l1, l2, l3, l4
  
  ! Note construct not usually used in G-Range.  This saves the toggled variable "in_winter" between months.
  save in_winter
  
  iunit = Rng(icell)%range_type
  temp_average = ( Globe(Rng(icell)%x, Rng(icell)%y)%max_temp + Globe(Rng(icell)%x, Rng(icell)%y)%min_temp ) / 2.0
  
  do ifacet = S_FACET, T_FACET
    ! Death of leaves
    if (Rng(icell)%leaf_carbon(ifacet) .gt. 0.0001) then
      ! Must account for deciduous leaves.  Century uses a death rate for non-deciduous months and another for deciduous months.  I am going to do the same, but the deciduous will be additive.
      ! Extremely complex if-then being replaced with a series of tests ...
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      l4 = .FALSE.
      if (temp_average .lt. Parms(iunit)%temperature_leaf_out_and_fall(2)) l1 = .TRUE.
      if (Rng(icell)%day_length_increasing .eq. .FALSE.) l2 = .TRUE.
      ! The following cannot use .eq. months as CENTURY does.  In CENTURY, a SAVE stores the state of "drop leaves" between months.  
      ! I will edit it to more closely match CENTURY, so that IN_WINTER *is* saved
      if (in_winter .eq. .FALSE. .and. Globe(Rng(icell)%x, Rng(icell)%y)%latitude .gt. 0. .and. month .eq. 12) l3 = .TRUE.
      if (in_winter .eq. .FALSE. .and. Globe(Rng(icell)%x, Rng(icell)%y)%latitude .le. 0. .and. month .eq. 6) l4 = .TRUE.
      if ((l1 .eq. .TRUE. .and. l2 .eq. .TRUE.) .or. l3 .eq. .TRUE. .or. l4 .eq. .TRUE.) then
        in_winter = .TRUE.
      end if
      ! Do the base line leaf death rate
      death_carbon = Rng(icell)%leaf_carbon(ifacet) * Parms(iunit)%leaf_death_rate(ifacet)                ! Century has an option to adjust death rate associated with element concentration, but I will skip that.
      ! And incorporate deciduous death.  That's the proportion that are deciduous times the carbon times the death rate)
      if ((in_winter .eq. .TRUE.) .and. (Rng(icell)%prop_annual_decid(ifacet) .gt. 0.001)) then
        death_carbon = death_carbon + ( Rng(icell)%prop_annual_decid(ifacet) * &
                        ( Rng(icell)%leaf_carbon(ifacet) * Parms(iunit)%death_rate_of_deciduous_leaves ) )
      end if
      if (Parms(iunit)%drought_deciduous(ifacet) .gt. 0.0001) then
        death_carbon = death_carbon + ( ( Rng(icell)%leaf_carbon(ifacet) * (1. - Rng(icell)%water_function) * &
                       Parms(iunit)%death_rate_of_deciduous_leaves ) * Parms(iunit)%drought_deciduous(ifacet) )
      end if
      Rng(icell)%dead_standing_carbon(ifacet) = Rng(icell)%dead_standing_carbon(ifacet) + death_carbon
      ! Make sure the following is only a placeholder.  Death_Carbon should not accumulate twice.   If it is happening, only use standing dead carbon.
      Rng(icell)%dead_leaf_carbon(ifacet) = Rng(icell)%dead_leaf_carbon(ifacet) + death_carbon
      ! Calculate the proportion of carbon that has been killed
      ! A different pathway than in Century.  I may misinterpret, but perhaps it is just simplier given the structure I have used.
      if (Rng(icell)%leaf_carbon(ifacet) .gt. 0.0001) then
        c_prop = death_carbon / Rng(icell)%leaf_carbon(ifacet)
      else
        c_prop = 0.0
      end if
      Rng(icell)%leaf_carbon(ifacet) = Rng(icell)%leaf_carbon(ifacet) - death_carbon
      to_storage = Rng(icell)%leaf_nitrogen(ifacet) * c_prop * Parms(iunit)%fraction_woody_leaf_n_translocated      
      Rng(icell)%stored_nitrogen(ifacet) = Rng(icell)%stored_nitrogen(ifacet) + to_storage
      death_nitrogen = ( Rng(icell)%leaf_nitrogen(ifacet) * c_prop ) - to_storage
      Rng(icell)%leaf_nitrogen(ifacet) = Rng(icell)%leaf_nitrogen(ifacet) - death_nitrogen
      Rng(icell)%dead_standing_nitrogen(ifacet) = Rng(icell)%dead_standing_nitrogen(ifacet) + death_nitrogen
      Rng(icell)%dead_leaf_nitrogen(ifacet) = Rng(icell)%dead_leaf_nitrogen(ifacet) + death_nitrogen

      ! A portion of stored carbon for maintence respiration is loss associated with death of plant part.    ! FLOWS may be misnamed, or the wrong indicator.  
      Rng(icell)%respiration_flows(H_FACET) = Rng(icell)%respiration_flows(H_FACET) - &
                                          ( Rng(icell)%respiration_flows(H_FACET) * c_prop ) 

      ! Do flows to litter, keeping track of structural and metabolic components
      ! DON'T PASS TO LITTER UNTIL STANDING DEAD LEAVES FALL.  CENTURY doesn't use standing dead for woody plants, but leaves can take some time to fall, so I will use them here.
    end if

    ! Death of seeds
    if (Rng(icell)%seed_carbon(ifacet) .gt. 0.0) then
      ! Death of seeds
      death_carbon = Rng(icell)%seed_carbon(ifacet) * ( Parms(iunit)%fraction_seeds_not_germinated(ifacet) / 12.0 ) 
    else
      death_carbon = 0.0
    end if
    if (Rng(icell)%seed_nitrogen(ifacet) .gt. 0.0) then
      death_nitrogen = Rng(icell)%seed_nitrogen(ifacet) * ( Parms(iunit)%fraction_seeds_not_germinated(ifacet) / 12.0 )      
    else
      death_nitrogen = 0.0
    end if
    Rng(icell)%dead_seed_carbon(ifacet) = Rng(icell)%dead_seed_carbon(ifacet) + death_carbon 
    Rng(icell)%seed_carbon(ifacet) = Rng(icell)%seed_carbon(ifacet) - death_carbon
    Rng(icell)%dead_seed_nitrogen(ifacet) = Rng(icell)%dead_seed_nitrogen(ifacet) + death_nitrogen
    Rng(icell)%seed_nitrogen(ifacet) = Rng(icell)%seed_nitrogen(ifacet) - death_nitrogen
    ! Seeds won't play a role in respiration.  These seeds are destined for decomposition only.
    ! Seeds don't move to standing dead.  They are viable until they drop from the plant and a portion becomes litter.
    ! Do flows to litter for seeds, keeping track of structural and metabolic components
    ! Seed litter is already partitioned in the main decomposition routine.   Commented out here.

    ! Death of fine branches 
    if (Rng(icell)%fine_branch_carbon(ifacet) .gt. 0.0) then
      death_carbon = Rng(icell)%fine_branch_carbon(ifacet) * Parms(iunit)%fine_branch_death_rate(ifacet)
    else
      death_carbon = 0.0
    end if
    Rng(icell)%fine_branch_carbon(ifacet) = Rng(icell)%fine_branch_carbon(ifacet) - death_carbon
    Rng(icell)%dead_fine_branch_carbon(ifacet) = Rng(icell)%dead_fine_branch_carbon(ifacet) + death_carbon     
    ! Fine branches should remain in standing dead until fall
    Rng(icell)%dead_standing_carbon(ifacet) = Rng(icell)%dead_standing_carbon(ifacet) + death_carbon
    if (Rng(icell)%fine_branch_nitrogen(ifacet) .gt. 0.0001) then
      death_nitrogen = Rng(icell)%fine_branch_nitrogen(ifacet) * Parms(iunit)%fine_branch_death_rate(ifacet)
    else
      death_nitrogen = 0.0
    end if
    Rng(icell)%fine_branch_nitrogen(ifacet) = Rng(icell)%fine_branch_nitrogen(ifacet) - death_nitrogen
    Rng(icell)%dead_fine_branch_nitrogen(ifacet) = Rng(icell)%dead_fine_branch_nitrogen(ifacet) + death_nitrogen
    Rng(icell)%dead_standing_nitrogen(ifacet) = Rng(icell)%dead_standing_nitrogen(ifacet) + death_nitrogen
    ! WAIT until standing dead falls to litter before partitioning ... FINE BRANCHES JOIN STANDING DEAD

    ! Death of fine roots
    if (Rng(icell)%fine_root_carbon(ifacet) .gt. 0.0) then
      death_carbon = Rng(icell)%fine_root_carbon(ifacet) * Parms(iunit)%fine_root_death_rate(ifacet)
    else 
      death_carbon = 0.0
    end if
    Rng(icell)%fine_root_carbon(ifacet) = Rng(icell)%fine_root_carbon(ifacet) - death_carbon
    Rng(icell)%dead_fine_root_carbon(ifacet) = Rng(icell)%dead_fine_root_carbon(ifacet) + death_carbon
    if (Rng(icell)%fine_root_nitrogen(ifacet) .gt. 0.00001) then
      death_nitrogen = Rng(icell)%fine_root_nitrogen(ifacet) * Parms(iunit)%fine_root_death_rate(ifacet)
    else
      death_nitrogen = 0.0
    end if
    Rng(icell)%fine_root_nitrogen(ifacet) = Rng(icell)%fine_root_nitrogen(ifacet) - death_nitrogen
    Rng(icell)%dead_fine_root_nitrogen(ifacet) = Rng(icell)%dead_fine_root_nitrogen(ifacet) + death_nitrogen
    ! Do flows to litter, keeping track of structural and metabolic components
    ! Fine roots are already partitioned in the main decomposition module.  

    ! Death of coarse branches
    if (Rng(icell)%coarse_branch_carbon(ifacet) .gt. 0.00001) then
      death_carbon = Rng(icell)%coarse_branch_carbon(ifacet) * Parms(iunit)%coarse_branch_death_rate(ifacet) 
    else
      death_carbon = 0.0
    end if
    Rng(icell)%coarse_branch_carbon(ifacet) = Rng(icell)%coarse_branch_carbon(ifacet) - death_carbon
    Rng(icell)%dead_coarse_branch_carbon(ifacet) = Rng(icell)%dead_coarse_branch_carbon(ifacet) + death_carbon
    if (Rng(icell)%coarse_branch_nitrogen(ifacet) .gt. 0.00001) then
      death_nitrogen = Rng(icell)%coarse_branch_nitrogen(ifacet) * Parms(iunit)%coarse_branch_death_rate(ifacet)
    else
      death_nitrogen = 0.0
    end if
    Rng(icell)%coarse_branch_nitrogen(ifacet) = Rng(icell)%coarse_branch_nitrogen(ifacet) - death_nitrogen
    Rng(icell)%dead_coarse_branch_nitrogen(ifacet) = Rng(icell)%dead_coarse_branch_nitrogen(ifacet) + death_nitrogen

    ! Death of coarse root
    if (Rng(icell)%coarse_root_carbon(ifacet) .gt. 0.00001) then
      death_carbon = Rng(icell)%coarse_root_carbon(ifacet) * Parms(iunit)%coarse_root_death_rate(ifacet) 
    else
      death_carbon = 0.0
    end if
    Rng(icell)%coarse_root_carbon(ifacet) = Rng(icell)%coarse_root_carbon(ifacet) - death_carbon
    Rng(icell)%dead_coarse_root_carbon(ifacet) = Rng(icell)%dead_coarse_root_carbon(ifacet) + death_carbon
    if (Rng(icell)%coarse_root_nitrogen(ifacet) .gt. 0.00001) then
      death_nitrogen = Rng(icell)%coarse_root_nitrogen(ifacet) * Parms(iunit)%coarse_root_death_rate(ifacet)
    else
      death_nitrogen = 0.0
    end if
    Rng(icell)%coarse_root_nitrogen(ifacet) = Rng(icell)%coarse_root_nitrogen(ifacet) - death_nitrogen
    Rng(icell)%dead_coarse_root_nitrogen(ifacet) = Rng(icell)%dead_coarse_root_nitrogen(ifacet) + death_nitrogen
  end do

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'W_P_DTH')

end subroutine