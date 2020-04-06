subroutine Whole_Plant_Death (icell)
!**** Simulate whole plant deaths
!****
!**** To represent populations, I will use an area 1 km x 1 km, and explicit tallying of space-filling plants.  There is the
!**** option to do this for entire cells modeled, but their areas will vary across the globe.  Moreover, none of the results
!**** will be reported on a per-cell basis, so the results from the 1 km^2 area will be suitable for reporting covers and
!**** concentrations, etc.  The 1 km^2 area was selected because it is convenient and sufficiently large enough to minimize
!**** rounding effects.  For example, if the effective root area of a tree is 8 x 8 m, 15,625 trees could fit within the
!**** 1 km^2 area.  Many thousands of herbs could fit, etc.  I speak of root area rather than volume because the soil depths
!**** and layers define volume (i.e., herbs have access to the first two layers, shrubs up to layer three, trees all four 
!**** layers).  Incidentally, an herb that occupies 0.2 x 0.2 m would fill the entire 1 km^2 area with 25,000,000 individuals,
!**** so there is no need for specially defined storage spaces (double and the like)
!****
!****
!****
!**** R. Boone   Last modified: Feb 13, 2013
  use Parameter_Vars
  use Structures
  implicit none
  
  real temperature, death_rate, data_val(4), temp_rate
  real linear, harvest(1)
  real proportion_cell_burned
  real proportion_plants_killed
  integer icell, ifacet, ilayer, iunit, i
  integer died(6)    ! temporary for debugging
  
  iunit = Rng(icell)%range_type
  temperature = ( Globe(Rng(icell)%x, Rng(icell)%y)%min_temp + Globe(Rng(icell)%x, Rng(icell)%y)%max_temp ) / 2.0
    
  ! Note that this is whole plant death and removal.  This is not related to standing dead on scenecent perennials, etc.  But to count flows here
  ! or not?  There is a risk of double-counting flows of nutrients and stocks.  NO ... nutrient flows are on a cell-basis, not captured here.
    
  ! Death occurs mostly when plants are not dormant.  Here I will use temperature as a surrogate for that.  I could use phenology, or a dormancy flag. 
  if (temperature .gt. 0.0) then
     do ifacet = 1, FACETS   
       death_rate = Parms(iunit)%nominal_plant_death_rate(ifacet)
       ! Probably pass correctly, but to decrease risk of error ...
       ! Those at the facet level are handled outside of the layers level
       do i = 1,4
         data_val(i) = Parms(iunit)%water_effect_on_death_rate(ifacet,i)
       end do
       death_rate = death_rate + linear(Rng(icell)%ratio_water_pet, data_val, 2) 
       do i = 1,4
         data_val(i) = Parms(iunit)%grazing_effect_on_death_rate(ifacet,i)
       end do
       death_rate = death_rate + linear(Rng(icell)%fraction_live_removed_grazing, data_val, 2) 
       temp_rate = death_rate               ! Save the death rate up to this point, so that recalculating based on LAI below won't keep adding to death rate for the different layers
       ! Incorproate annuals.  Do so after the season has ended and standing dead has fallen to litter, etc. The best time to account for 
       ! annual death may be the beginning of the following year, when phenology is reset to 0.  
       if (month .eq. Parms(iunit)%month_to_remove_annuals .and. ifacet .eq. H_FACET) then
         death_rate = death_rate + Rng(icell)%prop_annual_decid(H_FACET)
       end if
       
       ! Plants die regardless of their placement, whether in the understory of another plant, or defining a facet.   The rate is the same, except for LAI effects.  
       ! Kill the plants ...
       ! Shading effect on death rate is at the facet level, but the leaf area index is at the layer level, so store 
       ! death rate in temp_rate and use that in the following functions.
       do i = 1,4
         data_val(i) = Parms(iunit)%shading_effect_on_death_rate(ifacet,i)
       end do
       select case (ifacet)
         case (H_FACET)
           ! Do the three herbaceous layers, open, under shrubs, and under trees
           do ilayer = H_LYR, H_T_LYR
             death_rate = temp_rate + linear(Rng(icell)%leaf_area_index(H_FACET), data_val, 2) 
             death_rate = min(1.0, death_rate)
             death_rate = max(0.0, death_rate)
             died(ilayer) = Rng(icell)%total_population(ilayer) * death_rate
             Rng(icell)%total_population(ilayer) = Rng(icell)%total_population(ilayer) - &
                                                 ( Rng(icell)%total_population(ilayer) * death_rate )
           end do
         case (S_FACET)
           ! Do the two shrub layers, open, and under trees
           do ilayer = S_LYR, S_T_LYR
             death_rate = temp_rate + linear(Rng(icell)%leaf_area_index(S_FACET), data_val, 2) 
             death_rate = min(1.0, death_rate)
             death_rate = max(0.0, death_rate)             
             died(ilayer) = Rng(icell)%total_population(ilayer) * death_rate
             Rng(icell)%total_population(ilayer) = Rng(icell)%total_population(ilayer) - &
                                                 ( Rng(icell)%total_population(ilayer) * death_rate )
           end do
         case (T_FACET)
           ! Do the tree layer, open
           ilayer = T_LYR
           death_rate = temp_rate + linear(Rng(icell)%leaf_area_index(T_FACET), data_val, 2) 
           death_rate = min(1.0, death_rate)
           death_rate = max(0.0, death_rate)     
           died(ilayer) = Rng(icell)%total_population(ilayer) * death_rate
           Rng(icell)%total_population(ilayer) = Rng(icell)%total_population(ilayer) - &
                                                 ( Rng(icell)%total_population(ilayer) * death_rate )
       end select 
     end do

     ! Now calculate facet covers.  This is much streamlined with the 6-layer approach than the 3-facet only approach.
     ! ilayers 1, 4, and 6 are the overstory layers defining facets directly.     
     ! TREES
     Rng(icell)%facet_cover(T_FACET) = ( Rng(icell)%total_population(T_LYR) * Parms(iunit)%indiv_plant_area(T_FACET) ) / REF_AREA
     ! SHRUBS
     Rng(icell)%facet_cover(S_FACET) = ( Rng(icell)%total_population(S_LYR) * Parms(iunit)%indiv_plant_area(S_FACET) ) / REF_AREA
     ! HERBS
     Rng(icell)%facet_cover(H_FACET) = ( Rng(icell)%total_population(H_LYR) * Parms(iunit)%indiv_plant_area(H_FACET) ) / REF_AREA
     
     ! And update bare ground proportion.
     Rng(icell)%bare_cover = (1.0 - ( Rng(icell)%facet_cover(T_FACET) + Rng(icell)%facet_cover(S_FACET) + &
                             Rng(icell)%facet_cover(H_FACET) ) )  
995 format (A40, I5, 6(I15) )
  
   end if
1000 format(A30, I10, I10, I10, F10.7)

  ! Fire will be handled separately, because fire may occur when temperatures are below freezing.  A little duplicative, but no matter.
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
        
    ! If some of the cell is to burn, do that, removing whole dead plants
    if (proportion_cell_burned .gt. 0.0) then
      do ifacet = 1, FACETS
        ! Use fire severity, already calculated in modeling plant part death.
        ! Calculate the proportion of plants to be killed.
        data_val(2) = 0.0 
        data_val(4) = 1.0 
        data_val(1) = Parms(iunit)%fraction_plants_burned_dead(ifacet, 1)
        data_val(3) = Parms(iunit)%fraction_plants_burned_dead(ifacet, 2)
        proportion_plants_killed = linear(Rng(icell)%fire_severity, data_val, 2) * proportion_cell_burned

        select case (ifacet)
          case (H_FACET)
            ilayer = H_LYR
            ! Porportion of plants to be killed by fire ... 
            Rng(icell)%total_population(ilayer) = Rng(icell)%total_population(ilayer) - &
                                                 ( Rng(icell)%total_population(ilayer) * proportion_plants_killed )
            ilayer = H_S_LYR
            ! Porportion of plants to be killed by fire ... 
            Rng(icell)%total_population(ilayer) = Rng(icell)%total_population(ilayer) - &
                                                 ( Rng(icell)%total_population(ilayer) * proportion_plants_killed )
            ilayer = H_T_LYR
            ! Porportion of plants to be killed by fire ... 
            Rng(icell)%total_population(ilayer) = Rng(icell)%total_population(ilayer) - &
                                                 ( Rng(icell)%total_population(ilayer) * proportion_plants_killed )
          case (S_FACET)
            ilayer = S_LYR
            ! Porportion of plants to be killed by fire ... 
            Rng(icell)%total_population(ilayer) = Rng(icell)%total_population(ilayer) - &
                                                 ( Rng(icell)%total_population(ilayer) * proportion_plants_killed )
            ilayer = S_T_LYR
            ! Porportion of plants to be killed by fire ... 
            Rng(icell)%total_population(ilayer) = Rng(icell)%total_population(ilayer) - &
                                                 ( Rng(icell)%total_population(ilayer) * proportion_plants_killed )
          case (T_FACET)
            ilayer = T_LYR
            ! Porportion of plants to be killed by fire ... 
            Rng(icell)%total_population(ilayer) = Rng(icell)%total_population(ilayer) - &
                                                 ( Rng(icell)%total_population(ilayer) * proportion_plants_killed )
        end select
      end do
      ! Now recalculate facet covers.  This is much streamlined with the 6-layer approach than the 3-facet only approach.
      ! ilayers 1, 4, and 6 are the overstory layers defining facets directly.     
      ! TREES
      Rng(icell)%facet_cover(T_FACET) = ( Rng(icell)%total_population(T_LYR) * Parms(iunit)%indiv_plant_area(T_FACET) ) / REF_AREA
      ! SHRUBS
      Rng(icell)%facet_cover(S_FACET) = ( Rng(icell)%total_population(S_LYR) * Parms(iunit)%indiv_plant_area(S_FACET) ) / REF_AREA
      ! HERBS
      Rng(icell)%facet_cover(H_FACET) = ( Rng(icell)%total_population(H_LYR) * Parms(iunit)%indiv_plant_area(H_FACET) ) / REF_AREA
     
      ! And update bare ground proportion.
      Rng(icell)%bare_cover = (1.0 - ( Rng(icell)%facet_cover(T_FACET) + Rng(icell)%facet_cover(S_FACET) + &
                              Rng(icell)%facet_cover(H_FACET) ) )  
    end if

  end if

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'WHL_DTH')
  
end subroutine


subroutine Plant_Reproduction (icell)
!**** Simulate whole plant reproduction
!****
!****
!****
!**** R. Boone   Last modified: April 27, 2013
  use Parameter_Vars
  use Structures
  implicit none
  
  real temperature, data_val(4)
  real relative_establishment(V_LYRS)
  real linear, max_population
  integer icell, iunit, i
  integer pot_established(V_LYRS), tot_pop(FACETS)
!  real harvest(1)               ! Odd part of ABSOFT FORTRAN 95.   DEBUG  Just supporting a random echoing out.
  
  iunit = Rng(icell)%range_type
  temperature = ( Globe(Rng(icell)%x, Rng(icell)%y)%min_temp + Globe(Rng(icell)%x, Rng(icell)%y)%max_temp ) / 2.0
    
  ! Establishment occurs only when plants are not dormant.  Here I will use temperature as a surrogate for that.  I could use phenology, or a dormancy flag. 
  if (temperature .gt. 0.0) then
    ! Get the empty space, to determine establishment
    ! Note: As long as a plant is alive, it maintains control of that area.  This is based on herbaceous size classes
    ! only, since even though the plants reserve their full size, they can become established on a smaller area.  Really it is more
    ! about the logic, so that any given herb-size patch has a chance to be occupeid by an herb, shrub, or tree ... or nothing.
    ! Note that deaths have already been accounted for, and bare cover updated.  
    ! The only means of shifts in facet cover is through differential allotment of bare cover to the three facets.
    ! So I am going to calculate those values first, then do the three special cases after the main facets are attended to.
    
    relative_establishment(H_LYR)   = Parms(iunit)%relative_seed_production(H_FACET)
    relative_establishment(H_S_LYR) = Parms(iunit)%relative_seed_production(H_FACET)
    relative_establishment(H_T_LYR) = Parms(iunit)%relative_seed_production(H_FACET)
    relative_establishment(S_LYR)   = Parms(iunit)%relative_seed_production(S_FACET)
    relative_establishment(S_T_LYR) = Parms(iunit)%relative_seed_production(S_FACET)
    relative_establishment(T_LYR)   = Parms(iunit)%relative_seed_production(T_FACET)    

!    write(ECHO_FILE,'(A40, I5, F9.3, 6(F10.3))') ' REL_EST pre water FACETS: ',icell, &
!                             Rng(icell)%ratio_water_pet, relative_establishment    
    ! Do water limitations
    do i = 1,4  ;  data_val(i) = Parms(iunit)%water_effect_on_establish(H_FACET,i)  ;  end do
     if (Rng(icell)%ratio_water_pet .ge. 0.0 .and. Rng(icell)%ratio_water_pet .lt. V_LARGE) then
       ! Do nothing
     else
       Rng(icell)%ratio_water_pet = 0.0
     end if
     relative_establishment(H_LYR)   = relative_establishment(H_LYR)   * linear(Rng(icell)%ratio_water_pet, data_val, 2)
     relative_establishment(H_S_LYR) = relative_establishment(H_S_LYR) * linear(Rng(icell)%ratio_water_pet, data_val, 2)
     relative_establishment(H_T_LYR) = relative_establishment(H_T_LYR) * linear(Rng(icell)%ratio_water_pet, data_val, 2)
    do i = 1,4  ;  data_val(i) = Parms(iunit)%water_effect_on_establish(S_FACET,i)  ;  end do
     relative_establishment(S_LYR)   = relative_establishment(S_LYR)   * linear(Rng(icell)%ratio_water_pet, data_val, 2)
     relative_establishment(S_T_LYR) = relative_establishment(S_T_LYR) * linear(Rng(icell)%ratio_water_pet, data_val, 2)    
    do i = 1,4  ;  data_val(i) = Parms(iunit)%water_effect_on_establish(T_FACET,i)  ;  end do
     relative_establishment(T_LYR)   = relative_establishment(T_LYR)   * linear(Rng(icell)%ratio_water_pet, data_val, 2)

    ! Do litter limitation
    do i = 1,4  ;  data_val(i) = Parms(iunit)%litter_effect_on_establish(H_FACET,i)  ;  end do
     if (Rng(icell)%total_litter_carbon(SURFACE_INDEX) .ge. 0.0 .and. &
                                  Rng(icell)%total_litter_carbon(SURFACE_INDEX) .lt. V_LARGE) then
       ! Do nothing
     else
       Rng(icell)%total_litter_carbon(SURFACE_INDEX) = 0.0
     end if    
     relative_establishment(H_LYR)   = relative_establishment(H_LYR) * &
                                       linear((Rng(icell)%total_litter_carbon(SURFACE_INDEX) * 2.5), data_val, 2)
     relative_establishment(H_S_LYR) = relative_establishment(H_S_LYR) * &
                                       linear((Rng(icell)%total_litter_carbon(SURFACE_INDEX) * 2.5), data_val, 2)
     relative_establishment(H_T_LYR) = relative_establishment(H_T_LYR) * &
                                       linear((Rng(icell)%total_litter_carbon(SURFACE_INDEX) * 2.5), data_val, 2)
    do i = 1,4  ;  data_val(i) = Parms(iunit)%litter_effect_on_establish(S_FACET,i)  ;  end do
     relative_establishment(S_LYR)   = relative_establishment(S_LYR) * &
                                       linear((Rng(icell)%total_litter_carbon(SURFACE_INDEX) * 2.5), data_val, 2)
     relative_establishment(S_T_LYR) = relative_establishment(S_T_LYR) * &
                                       linear((Rng(icell)%total_litter_carbon(SURFACE_INDEX) * 2.5), data_val, 2)
    do i = 1,4  ;  data_val(i) = Parms(iunit)%litter_effect_on_establish(T_FACET,i)  ;  end do
     relative_establishment(T_LYR)   = relative_establishment(T_LYR) * &
                                       linear((Rng(icell)%total_litter_carbon(SURFACE_INDEX) * 2.5), data_val, 2)

    ! Do limitation due to herbaceous root biomass limiations.  Only empty patches are candidates for establishment, but that
    !    is taking the logic too rigorously.  Root biomass will play a role, as represented in the non-spatial portion of the model.
    do i = 1,4  ;  data_val(i) = Parms(iunit)%herb_root_effect_on_establish(H_FACET,i)  ;  end do
     if (Rng(icell)%fine_root_carbon(H_FACET) .ge. 0.0 .and. Rng(icell)%fine_root_carbon(H_FACET) .lt. V_LARGE) then
       ! Do nothing
     else
       Rng(icell)%fine_root_carbon(H_FACET) = 0.0
     end if
     relative_establishment(H_LYR)   = relative_establishment(H_LYR) * &
                                       linear((Rng(icell)%fine_root_carbon(H_FACET) * 2.5), data_val, 2)
     relative_establishment(H_S_LYR) = relative_establishment(H_S_LYR) * &
                                       linear((Rng(icell)%fine_root_carbon(H_FACET) * 2.5), data_val, 2)
     relative_establishment(H_T_LYR) = relative_establishment(H_T_LYR) * &
                                       linear((Rng(icell)%fine_root_carbon(H_FACET) * 2.5), data_val, 2)
    do i = 1,4  ;  data_val(i) = Parms(iunit)%herb_root_effect_on_establish(S_FACET,i)  ;  end do
     relative_establishment(S_LYR)   = relative_establishment(S_LYR) * &
                                       linear((Rng(icell)%fine_root_carbon(H_FACET) * 2.5), data_val, 2)
     relative_establishment(S_T_LYR) = relative_establishment(S_T_LYR) * &
                                       linear((Rng(icell)%fine_root_carbon(H_FACET) * 2.5), data_val, 2)
    do i = 1,4  ;  data_val(i) = Parms(iunit)%herb_root_effect_on_establish(T_FACET,i)  ;  end do
     relative_establishment(T_LYR)   = relative_establishment(T_LYR) * &
                                       linear((Rng(icell)%fine_root_carbon(H_FACET) * 2.5), data_val, 2)
                                       
!    write(ECHO_FILE,998) ' REL_EST pre wood FACETS: ',icell,Rng(icell)%facet_cover(S_FACET), &
!                                                          Rng(icell)%facet_cover(T_FACET), relative_establishment
    ! Do woody cover limitation ... this is applicable to most of the entries, but not all.
    do i = 1,4  ;  data_val(i) = Parms(iunit)%woody_cover_effect_on_establish(H_FACET,i)  ;  end do
     relative_establishment(H_S_LYR) = relative_establishment(H_S_LYR) * &
                                       linear(Rng(icell)%facet_cover(S_FACET), data_val, 2)
     relative_establishment(H_T_LYR) = relative_establishment(H_T_LYR) * &
                                       linear(Rng(icell)%facet_cover(T_FACET), data_val, 2)
    do i = 1,4  ;  data_val(i) = Parms(iunit)%woody_cover_effect_on_establish(S_FACET,i)  ;  end do
     relative_establishment(S_LYR)   = relative_establishment(S_LYR) * &
                                       linear(Rng(icell)%facet_cover(S_FACET), data_val, 2)
     relative_establishment(S_T_LYR) = relative_establishment(S_T_LYR) * &
                                       linear(Rng(icell)%facet_cover(T_FACET), data_val, 2)
    do i = 1,4  ;  data_val(i) = Parms(iunit)%woody_cover_effect_on_establish(T_FACET,i)  ;  end do
     relative_establishment(T_LYR)   = relative_establishment(T_LYR) * &
                                       linear(Rng(icell)%facet_cover(T_FACET), data_val, 2)

    ! Calculate total seed production, to allow competition between facets for the empty patches.
    ! This needs to be based on populations, rather than facet area, because understory plants may produce seed as well.
    ! Baseline relative seed production, incorporating total population size for the facet.
    tot_pop(H_FACET) = Rng(icell)%total_population(H_LYR) + Rng(icell)%total_population(H_S_LYR) + &
                       Rng(icell)%total_population(H_T_LYR) 
    tot_pop(S_FACET) = Rng(icell)%total_population(S_LYR) + Rng(icell)%total_population(S_T_LYR)
    tot_pop(T_FACET) = Rng(icell)%total_population(T_LYR)

    ! Potential plants established
    pot_established(H_LYR)   = tot_pop(H_FACET) * relative_establishment(H_LYR)
    pot_established(H_S_LYR) = tot_pop(H_FACET) * relative_establishment(H_S_LYR)
    pot_established(H_T_LYR) = tot_pop(H_FACET) * relative_establishment(H_T_LYR)
    pot_established(S_LYR)   = tot_pop(S_FACET) * relative_establishment(S_LYR)
    pot_established(S_T_LYR) = tot_pop(S_FACET) * relative_establishment(S_T_LYR)
    pot_established(T_LYR)   = tot_pop(T_FACET) * relative_establishment(T_LYR)
995 format (A40, I5, 6(I15) )
    ! Now plant the plants that will fit. The ordering will be critical here.  
    ! First, trees will take precedence, shading out the other facets, and recognizing that woody plant limitations were considered above.
    ! Second, shrubs will take precedence.  Third, herbs will be planted.
    ! Plants are being put in place full-size.
    ! Recall that death is already simulated at this point - facets will not shrink, but they may grow.  The order in which they grow is relavant.
    ! TREE
    Rng(icell)%total_population(T_LYR) = Rng(icell)%total_population(T_LYR) + pot_established(T_LYR)
    if (Rng(icell)%total_population(T_LYR) .gt. Parms(iunit)%pot_population(T_FACET)) then
      Rng(icell)%total_population(T_LYR) = Parms(iunit)%pot_population(T_FACET)
    end if
    Rng(icell)%facet_cover(T_FACET) = ( Rng(icell)%total_population(T_LYR) * Parms(iunit)%indiv_plant_area(T_FACET) ) / REF_AREA
    ! SHRUB
    ! First, adjust pot_population(S_FACET) to incorporate the tree facet.
    max_population = Parms(iunit)%pot_population(S_FACET) * (1.0 - Rng(icell)%facet_cover(T_FACET))                 ! The pot_pop is at 100%.  So if trees are 30%, then pot_pop * 0.70 would be the maximum shrubs could occupy.  Again, an ordered approach, with trees, then shrubs, then herbs.
    Rng(icell)%total_population(S_LYR) = Rng(icell)%total_population(S_LYR) + pot_established(S_LYR)
    if (Rng(icell)%total_population(S_LYR) .gt. max_population) then
      Rng(icell)%total_population(S_LYR) = max_population
    end if
    Rng(icell)%facet_cover(S_FACET) = ( Rng(icell)%total_population(S_LYR) * Parms(iunit)%indiv_plant_area(S_FACET) ) / REF_AREA
    ! HERB
    ! First, adjust pot_population(H_FACET) to incorporate the tree and shrub facets. 
    max_population = Parms(iunit)%pot_population(H_FACET) * (1.0 - ( Rng(icell)%facet_cover(T_FACET) + &
                     Rng(icell)%facet_cover(S_FACET) ) )                                                             ! The pot_pop is at 100%.  So if trees are 30%, then pot_pop * 0.70 would be the maximum shrubs could occupy.  Again, an ordered approach, with trees, then shrubs, then herbs.
    Rng(icell)%total_population(H_LYR) = Rng(icell)%total_population(H_LYR) + pot_established(H_LYR)
    if (Rng(icell)%total_population(H_LYR) .gt. max_population) then
      Rng(icell)%total_population(H_LYR) = max_population
    end if    
    Rng(icell)%facet_cover(H_FACET) = ( Rng(icell)%total_population(H_LYR) * Parms(iunit)%indiv_plant_area(H_FACET) ) / REF_AREA
    ! HERBS UNDER SHRUBS
    max_population = ( Parms(iunit)%pot_population(H_FACET) * &
                       Rng(icell)%facet_cover(S_FACET) * Parms(iunit)%indiv_plant_area(S_FACET) ) / &
                       Parms(iunit)%indiv_plant_area(H_FACET)                                                          ! Parameter, so no division by 0 likely.  The same is true below.
    Rng(icell)%total_population(H_S_LYR) = Rng(icell)%total_population(H_S_LYR) + pot_established(H_S_LYR)
    if (Rng(icell)%total_population(H_S_LYR) .gt. max_population) then
      Rng(icell)%total_population(H_S_LYR) = max_population
    end if
    ! SHRUBS UNDER TREES
    max_population = ( Parms(iunit)%pot_population(S_FACET) * &
                       Rng(icell)%facet_cover(T_FACET) * Parms(iunit)%indiv_plant_area(T_FACET) ) / &
                       Parms(iunit)%indiv_plant_area(S_FACET)
    Rng(icell)%total_population(S_T_LYR) = Rng(icell)%total_population(S_T_LYR) + pot_established(S_T_LYR)
    if (Rng(icell)%total_population(S_T_LYR) .gt. max_population) then
      Rng(icell)%total_population(S_T_LYR) = max_population
    end if
    ! HERBS UNDER TREES
    max_population = ( Parms(iunit)%pot_population(H_FACET) * &
                       Rng(icell)%facet_cover(T_FACET) * Parms(iunit)%indiv_plant_area(T_FACET) ) / &
                       Parms(iunit)%indiv_plant_area(H_FACET)
    Rng(icell)%total_population(H_T_LYR) = Rng(icell)%total_population(H_T_LYR) + pot_established(H_T_LYR)
    if (Rng(icell)%total_population(H_T_LYR) .gt. max_population) then
      Rng(icell)%total_population(H_T_LYR) = max_population
    end if

    ! And update bare ground proportion.
    Rng(icell)%bare_cover = (1.0 - ( Rng(icell)%facet_cover(T_FACET) + Rng(icell)%facet_cover(S_FACET) + &
                            Rng(icell)%facet_cover(H_FACET) ) )  
999 format (A40, I5, F9.3, 6(I15) )
998 format (A40, I5, 2(F9.3), 6(F10.7) )
997 format (A40, I5, 6(I14))  

    if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'PL_REPO')
  
  end if
  
end subroutine