subroutine Potential_Production (icell)
! **** Calculate potential production and the attributes that go into potential production.
! **** This routine also includes a calculation of surface soil temperatures.
! **** POTPROD uses a crop system and forest system approach.   
! **** (In CENTURY, cancvr is passed at each call, but I am using the structured approach)
! **** 
! **** R. Boone   Last modified:  March 31, 2011
  use Parameter_Vars
  use Structures
  implicit none

  real    biomass, woody_biomass, tmin, tmax
  real    temp_min_melt, temp_max_melt
  real    temp_min_leaf, temp_min_wood, temp_max_leaf, temp_max_wood
  real    surface_litter_biomass, standing_dead_biomass, avg_wood_biomass
  real    avg_live_biomass
  real    maximum_soil_surface_temperature, minimum_soil_surface_temperature
  integer icell, iunit, ifacet

  iunit = Rng(icell)%range_type
  tmin = Globe(Rng(icell)%x, Rng(icell)%y)%min_temp
  tmax = Globe(Rng(icell)%x, Rng(icell)%y)%max_temp
    
! Calculate temperature ... in CENTURY, this uses three approaches (CURSYS), one for FOREST, on for SAVANNA, and one for GRASSLAND.  
!                           Forest isn't represented here.  I am seeking to avoid the SAVANNA/GRASSLAND split.  
! Live biomass
  surface_litter_biomass = ( Rng(icell)%litter_structural_carbon(SURFACE_INDEX) + &
                             Rng(icell)%litter_metabolic_carbon(SURFACE_INDEX) ) * 2.25                ! Averaged over types

  do ifacet = 1, FACETS
    select case (ifacet)
      case (H_FACET)
        avg_live_biomass = ( Rng(icell)%leaf_carbon(H_FACET) + Rng(icell)%seed_carbon(H_FACET) ) * 2.5
        avg_wood_biomass = 0.0
        standing_dead_biomass = Rng(icell)%dead_standing_carbon(H_FACET) * 2.5
      case (S_FACET)
        avg_live_biomass = ( Rng(icell)%leaf_carbon(S_FACET) + Rng(icell)%seed_carbon(S_FACET) ) * 2.5
        avg_wood_biomass = ( Rng(icell)%shrub_carbon(COARSE_BRANCH_INDEX) + Rng(icell)%shrub_carbon(FINE_BRANCH_INDEX) ) * 2.0
        standing_dead_biomass = Rng(icell)%dead_standing_carbon(S_FACET) * 2.5
      case (T_FACET)      
        avg_live_biomass = ( Rng(icell)%leaf_carbon(T_FACET) + Rng(icell)%seed_carbon(T_FACET) ) * 2.5                             
        avg_wood_biomass = ( Rng(icell)%tree_carbon(COARSE_BRANCH_INDEX) + Rng(icell)%tree_carbon(FINE_BRANCH_INDEX) ) * 2.0
        standing_dead_biomass = Rng(icell)%dead_standing_carbon(T_FACET) * 2.5
    end select
    ! ** Calculating soil surface temperatures
    ! Century makes a call to a routine to calculate soil surface temperature.  It is only used here, and this is brief, so I will merge that into here.
    ! Total biomass
    biomass = avg_live_biomass + standing_dead_biomass + &
              ( surface_litter_biomass * Parms(iunit)%litter_effect_on_soil_temp )
    biomass = min( biomass, Parms(iunit)%maximum_biomass_soil_temp )
    woody_biomass = min( avg_wood_biomass, 5000.0 )                             ! Number hardwired in CENTURY
    
    ! Maximum temperature with leaf shading
    temp_max_leaf = tmax + (25.4 / ( 1. + 18. * exp( -0.20 * tmax ))) * &                          ! 1 + avoids division by 0.
                    ( exp( Parms(iunit)%biomass_effect_on_max_soil_temp * biomass ) - 0.13 )
    ! Minimum temperature with leaf shading
    temp_min_leaf = tmin + ( Parms(iunit)%biomass_effect_on_min_soil_temp * biomass ) - 1.78 
    
    ! Maximum temperature with wood shading
    temp_max_wood = tmax + ( 25.4 / ( 1. + 18. * exp( -0.20 * tmax ))) * &                         ! 1 + avoids division by 0.
                    ( exp( Parms(iunit)%biomass_effect_on_max_soil_temp * 0.1 * woody_biomass ) - 0.13 )
    ! Minimum temperature with wood shading
    temp_min_wood = tmin + ( Parms(iunit)%biomass_effect_on_min_soil_temp * 0.1 * woody_biomass ) - 1.78
    
    maximum_soil_surface_temperature = min( temp_max_leaf, temp_max_wood )
    minimum_soil_surface_temperature = max( temp_min_leaf, temp_min_wood )
    
    ! Let soil surface temperature be affected by day length
    if (Rng(icell)%day_length .lt. 12.0) then
      temp_min_melt = (( 12.0 - Rng(icell)%day_length) * 3.0 + 12.0 ) / 24.0
    else
      temp_min_melt = (( 12.0 - Rng(icell)%day_length) * 1.2 + 12.0 ) / 24.0
    end if
    temp_min_melt = min(0.95, temp_min_melt)
    temp_min_melt = max(0.05, temp_min_melt)
    temp_max_melt = 1.0 - temp_min_melt
    Rng(icell)%soil_surface_temperature = temp_max_melt * maximum_soil_surface_temperature + &
                                                  temp_min_melt * minimum_soil_surface_temperature
  end do  
  ! End calculating soil surface temperatures

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'POT_PRD')
    
  ! Calculate potential production, using POTCRP as a guide, but not including cropping types.  An edit included POTTREE as well, but mostly already present in the model.
  call Potential(icell)
                                                  
end subroutine


subroutine Potential(icell)
!**** Calculate potential production given the information calculated in Potential_Production and other data.
!**** Potential draws from POTCRP, but without the specifics to cropping systems.
!**** WC (water_content) is calcualted elsewhere in CENTURY, but handiest here.
!****
!**** This has been modified to calculate potential for all six vegetation layers.
!****
!**** R. Boone      Last updated:  July 17, 2014
    use Parameter_Vars
    use Structures
    implicit none

    real shortwave
    real aisc, h2ogef, water_content, frac, bioc, bioprd, biof, fracrc
    real shading_modifier, woody_cover
    real intercept, slope
    real a, b, c, d
    real temp1, temp2, temp3, ratlc
    real prop_live_per_layer(V_LYRS), total_cover, w_cover
    
    integer icell, iunit, ifacet, ilyr

    iunit = Rng(icell)%range_type     ! Identifies the index to the Parm array structure, storing landscape units.

    ! From Century, the value for potential plant production is now calculated from the equation of a line whose intercept changes depending on water
    ! content based on soil type. 
    if ( Rng(icell)%pot_evap .ge. 0.01 ) then
      h2ogef = ( Rng(icell)%water_available(1) + Globe(Rng(icell)%x, Rng(icell)%y)%precip ) / Rng(icell)%pot_evap         ! Irrigation was not included, unlike CENTURY
    else
      h2ogef = 0.01
    endif

    water_content = Rng(icell)%field_capacity(1) - Rng(icell)%wilting_point(1)
    
    ! Doing PPRDWC in this subroutine, rather than another call to a function 
    ! Regression points were confirmed as 0.0, 1.0, and 0.8 (using the old method in Century to reduce parameters)
    intercept = Parms(iunit)%ppt_regression_points(1) + (Parms(iunit)%ppt_regression_points(2) * water_content)
    if (Parms(iunit)%ppt_regression_points(3) .ne. intercept) then
      slope = 1.0 / (Parms(iunit)%ppt_regression_points(3) - intercept)                                               
    else
      slope = 1.0
    end if
    ! Do the correction, altering h2ogef based on these corrections.  h2ogef is both x (in) and pprdwc (out) in the PPRDWC function in CENTURY
    h2ogef = 1.0 + slope * (h2ogef - Parms(iunit)%ppt_regression_points(3))
    if (h2ogef .gt. 1.0) then
      h2ogef = 1.0
    elseif (h2ogef .lt. 0.01) then
      h2ogef = 0.01
    endif    

    ! Calculate how much live aboveground biomass is in each vegetation layer, and use that as a guide to distribute production
    ! Get the total cover, to allow ignoring bare ground.
    do ilyr = 1, V_LYRS
      prop_live_per_layer(ilyr) = 0.0                              ! Filling the arrays in case something gets skipped below.
    end do
   
    total_cover =   Rng(icell)%facet_cover(H_FACET) + Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(T_FACET)
    if (total_cover .gt. 0.000001) then
      ! Using an approach that provides production estimates for each facet independently. 
      ! Also accounting for not looking at bare ground.
      prop_live_per_layer(H_LYR) = 1.0
      prop_live_per_layer(H_S_LYR) = Rng(icell)%facet_cover(S_FACET) / total_cover
      prop_live_per_layer(H_T_LYR) = Rng(icell)%facet_cover(T_FACET) / total_cover
      w_cover = Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(T_FACET)
      if (w_cover .gt. 0.000001) then
        prop_live_per_layer(S_LYR) = 1.0
        prop_live_per_layer(S_T_LYR) = Rng(icell)%facet_cover(T_FACET) / w_cover
        prop_live_per_layer(T_LYR) = 1.0
      end if

      do ilyr = 1, V_LYRS
        ! Compute shading modifier.  First, set woody cover 
        select case (ilyr)
          case (H_LYR)
            ifacet = H_FACET
            woody_cover = 0.0
            aisc = 0.0
          case (H_S_LYR)                                              ! NOTE that in the following, spatial cover is substituting for what would normally be a density measure of leaves at any one place.  Use LAI instead?  No, canopy cover is used in Century.  Across the entire landscape cell, this measure is appropriate.
            ifacet = H_FACET
            woody_cover = Rng(icell)%facet_cover(S_FACET)                                            ! Facet cover should never go to 0, as it is counter to the definition.  Perhaps include a catch-all in MISC_MATERIAL that makes sure a few shrubs are present, a few trees are present, in any cell.
            if (Rng(icell)%shrub_carbon(LEAF_INDEX) .lt. 0.00001) then
              aisc = 0.0
            else
              aisc = 5.0 * exp(-0.0035 * (Rng(icell)%shrub_carbon(LEAF_INDEX) * 2.5) / woody_cover+0.00000001)    ! Shading by seeds ignored.
            end if
          case (H_T_LYR)
            ifacet = H_FACET
            woody_cover = Rng(icell)%facet_cover(T_FACET)
            if (Rng(icell)%tree_carbon(LEAF_INDEX) .lt. 0.00001) then
              aisc = 0.0
            else
              aisc = 5.0 * exp(-0.0035 * (Rng(icell)%tree_carbon(LEAF_INDEX) * 2.5) / woody_cover+0.00000001)
            end if
          case (S_LYR)
            ifacet = S_FACET
            woody_cover = 0.0                                                                        ! Assumes that shrubs aren't shaded, although they do shade themselves.  But really, all these types shade themselves.  Perhaps adjust.
            aisc = 0.0
          case (S_T_LYR)
            ifacet = S_FACET
            woody_cover = Rng(icell)%facet_cover(T_FACET)
            if (Rng(icell)%tree_carbon(LEAF_INDEX) .lt. 0.00001) then
              aisc = 0.0
            else
              aisc = 5.0 * exp(-0.0035 * (Rng(icell)%tree_carbon(LEAF_INDEX) * 2.5) / woody_cover+0.00000001)           
            end if
          case (T_LYR)
            ifacet = T_FACET
            woody_cover = 0.0                                                                        ! Assumes trees don't shade themselves.
            aisc = 0.0
        end select
        shading_modifier = ( 1.0 - woody_cover )+( woody_cover * ( aisc / (aisc + 1.) ) )
        ! I suspect the modifier should be less than 1.   Check for this?
      
        ! Estimate plant production
        if (Rng(icell)%soil_surface_temperature .gt. 0.0) then
          ! Calculate temperature effect on growth
          ! Account for removal of litter effects on soil temperature as it drives plant production
          ! Century recalculates min, max, and average soil surface temperatures.  I have those already, so using the existing.
         
          ! Century uses a function call, but I will collapse it to here.  X is ctemp or average soil surface temperature.  a,b,c,d are PPDF values
          a = Parms(iunit)%temperature_production(1)
          b = Parms(iunit)%temperature_production(2)
          c = Parms(iunit)%temperature_production(3)
          d = Parms(iunit)%temperature_production(4)
        
          frac = ( b - Rng(icell)%soil_surface_temperature ) / ( b - a )                                     ! Based mostly on parameters, so division by 0 unlikely.
          Rng(icell)%potential_production = 0.0
          ! The following appears appropriate for both herbs and woodies, based on Century documentation.
          if (frac .gt. 0.0) then
            Rng(icell)%potential_production = exp(c/d * (1.0 - frac**d)) * (frac**c)
          endif

          ! Calculate the potential effect of standing dead on plant growth, the effect of physical obstruction of litter and standing dead
          bioc = Rng(icell)%dead_standing_carbon(ifacet) + 0.1 * Rng(icell)%litter_structural_carbon(1)
          if (bioc .le. 0.0) then
            bioc = 0.01
          end if
          if (bioc .gt. Parms(iunit)%maximum_biomass_soil_temp) then
            bioc = Parms(iunit)%maximum_biomass_soil_temp
          end if
          bioprd = 1. - ( bioc / (Parms(iunit)%standing_dead_production_halved + bioc))                      ! Parameters, so division by 0 unlikely.
  
          ! Calculate the effect of the ratio of live biomass to dead biomass on the reduction of potential growth rate.  The intercept of this equation (highest negative effect of dead plant biomass) is equal to bioprd when the ratio is zero.
          temp1 = (1. - bioprd)
          temp2 = temp1 * 0.75
          temp3 = temp1 * 0.25
          ratlc = Rng(icell)%leaf_carbon(H_FACET) / bioc                                                     ! Logic above prevents 0.
          if (ratlc .le. 1.0) then
            biof = bioprd + ( temp2 * ratlc )
          end if
          if (ratlc .gt. 1.0 .and. ratlc .le. 2.0) then
            biof = ( bioprd + temp2 ) + temp3 * ( ratlc-1. )
          end if
          if (ratlc .gt. 2.0) then
            biof = 1.0
          end if

          Rng(icell)%total_pot_production(ilyr) = shortwave(icell) * Parms(iunit)%radiation_production_coefficient * &
             Rng(icell)%potential_production * h2ogef * biof * shading_modifier * Rng(icell)%co2_effect_on_production(ifacet) * &
             prop_live_per_layer(ilyr)
!if (Rng(icell)%total_pot_production(2) .gt. 10000.0) then             
!  write(*,*) 'ZZ TOT_POT_PROD: ', icell, ilyr, month, Rng(icell)%total_pot_production(2), total_biomass, total_cover, & 
!                                  prop_live_per_layer(ilyr)
!  write(*,*) 'ZZ          S_W: ', icell, ilyr, month, shortwave(icell)
!  write(*,*) 'ZZ       h2ogef: ', icell, ilyr, month, h2ogef
!  write(*,*) 'ZZ         biof: ', icell, ilyr, month, biof
!  write(*,*) 'ZZ          S_M: ', icell, ilyr, month, shading_modifier
!  write(*,*) 'ZZ         PLPL: ', icell, ilyr, month, prop_live_per_layer(ilyr)
!  write(*,*) 'ZZ      RAD_P_C: ', icell, ilyr, month, Parms(iunit)%radiation_production_coefficient
!  write(*,*) 'ZZ          P_P: ', icell, ilyr, month, Rng(icell)%potential_production
!  write(*,*) 'ZZ          CO2: ', icell, ilyr, month, Rng(icell)%co2_effect_on_production(ifacet)
!end if
      
          ! Dynamic carbon allocation to compute root/shoot ratio
          if (Rng(icell)%total_pot_production(ilyr) .gt. 0.0) then
            ! call Crop_Dynamic_Carbon(Rng(icell)%root_shoot_ratio, fracrc)      I WON'T BE INCLUDING DYNAMIC CARBON ALLOCATION BETWEEN SHOOTS AND ROOTS FOR NOW.  TOO COMPLEX, MUST SIMPLIFY.
            fracrc = Parms(iunit)%fraction_carbon_to_roots(ifacet)               ! Gross simplifaction, but required for completion.
            ! Change root shoot ratio based on effects of co2
            ! The following can't be the same as effect on production.  Distruptive.  I will turn this off for now.  Specific to crops and trees, incidentally.
            ! Rng(icell)%root_shoot_ratio(ifacet) = Rng(icell)%root_shoot_ratio(ifacet) * Rng(icell)%co2_effect_on_production(ifacet)
          
            ! Allocate production
            Rng(icell)%belowground_pot_production(ilyr) = Rng(icell)%total_pot_production(ilyr) * fracrc
            Rng(icell)%aboveground_pot_production(ilyr) = Rng(icell)%total_pot_production(ilyr) - &
                                                        Rng(icell)%belowground_pot_production(ilyr)
                 
            ! Restrict production due to that taken by grazers
            call Grazing_Restrictions (icell, ifacet, ilyr)
        
            ! Update accumulators and compute potential C production
            ! Skipping this for now.  They seem to be detailed accumulators of C over the entire run, which may be of interest, but not now.        
          else
            ! No production this month ... total potential production is 0.0
            Rng(icell)%total_pot_production(ilyr) = 0.0
            Rng(icell)%aboveground_pot_production(ilyr) = 0.0
            Rng(icell)%belowground_pot_production(ilyr) = 0.0
          endif
        else
          ! No production this month ... too cold
          Rng(icell)%total_pot_production(ilyr) = 0.0
          Rng(icell)%aboveground_pot_production(ilyr) = 0.0
          Rng(icell)%belowground_pot_production(ilyr) = 0.0
        end if
      end do                 ! End of ilyr loop
    else
      ! No production this month ... there is either zero cover or zero biomass
      do ilyr = 1, V_LYRS
        prop_live_per_layer(ilyr) = 0.0                      ! Nothing alive, so no production.  Or no facet cover, so no production
        Rng(icell)%total_pot_production(ilyr) = 0.0                               ! Note this is biomass
        Rng(icell)%aboveground_pot_production(ilyr) = 0.0                         ! Note this is biomass
        Rng(icell)%belowground_pot_production(ilyr) = 0.0                         ! Note this is biomass.  Divide by 2.5 or 2 for carbon.
      end do
    end if  
    
    if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'POTENTL')

end subroutine


subroutine Grazing_Restrictions (icell, ifacet, ilyr)
!**** Grazing restrictions modifies production due to grazing pressures, which are calculated elsewhere.  
!**** Grazing_effects uses scores from 0 to 6 to describe the functional responses, right from CENTURY
!****
!**** R. Boone                 Last modified:  February 16, 2011
   use Parameter_Vars
   use Structures
   implicit none
   
   real agrd, bgrd, fracrmv, rtsht, bop, graze_mult
   integer icell, ifacet, ilyr, iunit
  
   iunit = Rng(icell)%range_type
  
   ! Only to make things more compressed
   agrd = Rng(icell)%aboveground_pot_production(ilyr)
   bgrd = Rng(icell)%belowground_pot_production(ilyr)
   fracrmv = Rng(icell)%fraction_live_removed_grazing
   if (agrd .gt. 0.0) then            
     rtsht = bgrd / agrd   
   else
     rtsht = 0.0
   end if
   graze_mult = Parms(iunit)%grazing_effect_multiplier

   select case (Parms(iunit)%grazing_effect)            ! Grazing effects 0 through 6 come right from CENTURY
!    case (0)   ! Grazing has no direct effect on production.  Captured in 'case default'
     case (1)   ! Linear impact of grazing on aboveground potential production
       agrd = (1-(2.21*fracrmv))*agrd
       if (agrd .lt. 0.02) then
         agrd = 0.02
       endif
       bgrd = rtsht * agrd
     case (2)   ! Quadratic impact of grazing on aboveground potential production and root:shoot ratio
       agrd = (1+(2.6*fracrmv-(5.83*(fracrmv**2))))*agrd
       if (agrd .lt. 0.02) then
         agrd = 0.02
       endif
       bop = rtsht + 3.05*fracrmv - 11.78*(fracrmv**2)
       if (bop .le. 0.01) then  
         bop = 0.01
       endif
       bgrd = agrd * bop            
     case (3)   ! Quadratic impact of grazing of grazing on root:shoot ratio
       bop = rtsht + 3.05*fracrmv - 11.78*(fracrmv**2)
       if (bop .le. 0.01) then  
         bop = 0.01
       endif
       bgrd = agrd * bop
     case (4)   ! Linear impact of grazing on root:shoot ratio
       bop = 1 - (fracrmv * graze_mult)
       bgrd = agrd * bop
     case (5)   ! Quadratic impact of grazing on aboveground potential production and linear impact on root:shoot ratio
       agrd = (1 + 2.6*fracrmv - (5.83*(fracrmv**2)))*agrd
       if (agrd .lt. 0.02) then
         agrd = 0.02
       endif
       bop = 1 - (fracrmv * graze_mult)
       bgrd = agrd * bop
     case (6)   ! Linear impact of grazing on aboveground potential production and root:shoot ratio
       agrd = (1 + 2.21*fracrmv) * agrd
       if (agrd .lt. 0.02) then
         agrd = 0.02
       endif
       bop = 1 - (fracrmv * graze_mult)
       bgrd = agrd * bop
     case default
       ! Do nothing.  grazing_effect = 0, and so values will be read and re-assigned without modification    
       ! This routine modifies the effects of grazing on above- and belowground potential production.  It removes no forage.
   end select
   Rng(icell)%aboveground_pot_production(ilyr) = agrd
   Rng(icell)%belowground_pot_production(ilyr) = bgrd

   Rng(icell)%total_pot_production(ilyr) = agrd + bgrd
   if (agrd .gt. 0.0) then
     Rng(icell)%root_shoot_ratio(ifacet) = bgrd/agrd                 
   else
     Rng(icell)%root_shoot_ratio(ifacet) = 1.0
   end if

   if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'GRZ_RST') 
         
end subroutine



subroutine Herb_Growth (icell)
!**** Growth builds from potential production estimates to calculate actual growth
!****
!**** R. Boone                 Last modified:  July 18, 2014
!****                                          Edited to avoid taking the mean of herb layers for mcprd calculation
!****                                          Edited to stop growth when phenology is 4.0
   use Parameter_Vars
   use Structures
   implicit none

   real tolerance, uptake(3), accum, mcprd(2), rimpct, cfrac(WOODY_PARTS), available_nitrogen                       ! Three cells in CFRAC won't be used, but dimensioning to 2 will cause it to fail, I suspect.
   real agfrac, bgfrac, resp_flow_shoots, resp_flow_roots, resp_temp_effect, cmrspflux(2), euf(2), amt, fsol, calcup
   real avg_total_pot_prod_carbon, avg_aground_pot_prod_carbon, avg_total_prod_limited_n
   integer icell, iunit, ilayer
   
   tolerance = 1.0E-30
   iunit = Rng(icell)%range_type

   uptake(N_STORE) = 0.0
   uptake(N_SOIL) = 0.0
   uptake(N_FIX) = 0.0
   accum = 0.0
   mcprd(SURFACE_INDEX) = 0.0
   mcprd(SOIL_INDEX) = 0.0
   Rng(icell)%maintain_respiration(H_FACET) = 0.0        
   Rng(icell)%respiration_flows(H_FACET) = 0.0
  
   ! Century includes flags set in the schedular to turn growth on or off.  I don't want to do that.  
   ! I wish to use degree-days, or topsoil available water to pet ratio, or temperature limits
   ! *** USE PHENOLOGY HERE?  USE DORMANCY INSTEAD?   DAY LENGTH (which is in RNG) ... no, not for herbs.
   if (Rng(icell)%ratio_water_pet > 0.0 .and. Globe(Rng(icell)%x, Rng(icell)%y)%min_temp .gt. 3.0 .and. &
       Rng(icell)%phenology(H_FACET) .lt. 3.999 ) then
     ! Calculate effect of root biomass on available nutrients          
     if ((Parms(iunit)%root_intercept_on_nutrients * Rng(icell)%fine_root_carbon(H_FACET) * 2.5) .gt. 33) then
       rimpct = 1.0
     else
       rimpct = (1.0 - Parms(iunit)%root_intercept_on_nutrients * &
               exp(-Parms(iunit)%root_effect_on_nutrients * Rng(icell)%fine_root_carbon(H_FACET) * 2.5))
     end if

     ! Calculate carbon fraction above and belowground
     if (Rng(icell)%total_pot_production(H_LYR) + Rng(icell)%total_pot_production(H_S_LYR) + &
         Rng(icell)%total_pot_production(H_T_LYR) .gt. 0.0) then
       cfrac(ABOVE) = Rng(icell)%aboveground_pot_production(H_FACET) / ( Rng(icell)%total_pot_production(H_LYR) + &
                      Rng(icell)%total_pot_production(H_S_LYR) + Rng(icell)%total_pot_production(H_T_LYR) )              ! Using ABOVE and BELOW in CFRAC when it is dimensioned as woody parts, but that is ok.
     else
       cfrac(ABOVE) = 0.0
     end if       
     cfrac(BELOW) = 1.0 - cfrac(ABOVE)

     available_nitrogen = 0.0
     do ilayer=1, SOIL_LAYERS                     ! Nutrients will be used to calculate mineral availablity for all layers, since they are only 15 cm each, 4, across the globe.  Recent CENTURY uses a parameter (CLAYPG) here.     
       available_nitrogen = available_nitrogen + Rng(icell)%mineral_nitrogen(ilayer)
     end do
     
     ! Determine actual production, restricted based on carbon to nitrogen ratios.  Note CFRAC & UPTAKE are arrays.
     call Restrict_Production (icell, H_FACET, 2, available_nitrogen, rimpct, cfrac, uptake)

     ! If growth occurs ...                       (Still stored in potential production ... move to actual production?)
     ! Get average potential production
     avg_total_pot_prod_carbon = ( ( Rng(icell)%total_pot_production(H_LYR) +  Rng(icell)%total_pot_production(H_S_LYR) + &
                                   Rng(icell)%total_pot_production(H_T_LYR) ) / 3.0 ) * 0.4 
     avg_aground_pot_prod_carbon = ( ( Rng(icell)%aboveground_pot_production(H_LYR) + &
            Rng(icell)%aboveground_pot_production(H_S_LYR) + Rng(icell)%aboveground_pot_production(H_T_LYR) ) / 3.0 ) * 0.4 
     if (avg_total_pot_prod_carbon .gt. 0.0) then         ! Wouldn't this be production limited by nitrogen?
       ! Compute nitrogen fixation which actually occurs and add to accumulator
       Rng(icell)%nitrogen_fixed(H_FACET) = Rng(icell)%nitrogen_fixed(H_FACET) + Rng(icell)%plant_nitrogen_fixed(H_FACET)
       ! Accumulators skipped for now.  Will be added as needed (Century includes so many it is bound to confuse)   EUPACC  SNFXAC  NFIXAC  TCNPRO
       
       ! Maintenance respiration calculations
       ! Growth of shoots
       if (avg_total_pot_prod_carbon .gt. 0.0) then
         agfrac = avg_aground_pot_prod_carbon / avg_total_pot_prod_carbon
       else
         agfrac = 0.0
       end if
       mcprd(ABOVE) = ( Rng(icell)%total_pot_prod_limited_by_n(H_LYR) + Rng(icell)%total_pot_prod_limited_by_n(H_S_LYR) + &
                        Rng(icell)%total_pot_prod_limited_by_n(H_T_LYR) ) * agfrac
       resp_flow_shoots = mcprd(ABOVE) * Parms(iunit)%fraction_npp_to_respiration(H_FACET)
       Rng(icell)%respiration_flows(H_FACET) = Rng(icell)%respiration_flows(H_FACET) + resp_flow_shoots
       Rng(icell)%leaf_carbon(H_FACET) = Rng(icell)%leaf_carbon(H_FACET) + &
                           ( ( 1.0 -  Parms(iunit)%fraction_aground_npp_to_seeds(H_FACET) ) * mcprd(ABOVE) )        ! Translated from CSHED call and complexity of CSRSNK and BGLCIS.  Not sure if interpretted correctly.
       Rng(icell)%seed_carbon(H_FACET) = Rng(icell)%seed_carbon(H_FACET) + &
                           ( Parms(iunit)%fraction_aground_npp_to_seeds(H_FACET) * mcprd(ABOVE) )                 ! Translated from CSHED call and complexity of CSRSNK and BGLCIS.  Not sure if interpretted correctly.
       Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink - mcprd(ABOVE)
       ! Growth of roots
       bgfrac = 1.0 - agfrac
       mcprd(BELOW) = ( ( Rng(icell)%total_pot_prod_limited_by_n(H_LYR) + Rng(icell)%total_pot_prod_limited_by_n(H_S_LYR) + &
                        Rng(icell)%total_pot_prod_limited_by_n(H_T_LYR) ) / 3.0 ) * bgfrac
       resp_flow_roots = mcprd(BELOW) * Parms(iunit)%fraction_npp_to_respiration(H_FACET)
       Rng(icell)%respiration_flows(H_FACET) = Rng(icell)%respiration_flows(H_FACET) + resp_flow_roots
       Rng(icell)%fine_root_carbon(H_FACET) = Rng(icell)%fine_root_carbon(H_FACET) + mcprd(BELOW)                  ! Translated from CSHED call and complexity of CSRSNK and BGLCIS.  Not sure if interpretted correctly.
       Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink - mcprd(BELOW)
       ! Store maintenance respiration to storage pool
       Rng(icell)%maintain_respiration(H_FACET) = Rng(icell)%maintain_respiration(H_FACET) + &
                                                     Rng(icell)%respiration_flows(H_FACET)
       Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink - Rng(icell)%respiration_flows(H_FACET)

       ! Maintenance respiration fluxes reduce maintenance respiration storage pool
       resp_temp_effect = 0.1 * exp(0.07 * Globe(Rng(icell)%x,Rng(icell)%y)%temperature_average)
       resp_temp_effect = min(1.0, resp_temp_effect)
       resp_temp_effect = max(0.0, resp_temp_effect)
       cmrspflux(ABOVE) = Parms(iunit)%herb_max_fraction_npp_to_respiration(ABOVE) * resp_temp_effect * &
                          Rng(icell)%leaf_carbon(H_FACET)
       cmrspflux(BELOW) = Parms(iunit)%herb_max_fraction_npp_to_respiration(BELOW) * resp_temp_effect * &
                          Rng(icell)%fine_root_carbon(H_FACET)
           
       Rng(icell)%respiration_annual(H_FACET) = Rng(icell)%respiration_annual(H_FACET) + cmrspflux(ABOVE) + cmrspflux(BELOW)
       Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink + cmrspflux(ABOVE)
       Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink + cmrspflux(BELOW)

       ! Get average potential production
       avg_total_prod_limited_n = ( Rng(icell)%total_pot_prod_limited_by_n(H_LYR) + &
             Rng(icell)%total_pot_prod_limited_by_n(H_S_LYR) + Rng(icell)%total_pot_production(H_T_LYR) ) * 2.5 
       ! Actual uptake
       if (avg_total_prod_limited_n .gt. 0.0) then
         euf(ABOVE) = Rng(icell)%e_up(H_FACET,ABOVE) / avg_total_prod_limited_n
         euf(BELOW) = Rng(icell)%e_up(H_FACET,BELOW) / avg_total_prod_limited_n
       else
         euf(ABOVE) = 0.0
         euf(BELOW) = 0.0
       end if
         
       ! Takeup nutrients from internal storage pool, and don't allow that if storage stored nitrogen (CRPSTG) is negative
       if (Rng(icell)%stored_nitrogen(H_FACET) .gt. 0.0) then
         amt = uptake(N_STORE) * euf(ABOVE)
         Rng(icell)%stored_nitrogen(H_FACET) = Rng(icell)%stored_nitrogen(H_FACET) - amt
         Rng(icell)%leaf_nitrogen(H_FACET) = Rng(icell)%leaf_nitrogen(H_FACET) + &
                           ( ( 1.0 -  Parms(iunit)%fraction_aground_npp_to_seeds(H_FACET) ) * amt )  
         Rng(icell)%seed_nitrogen(H_FACET) = Rng(icell)%seed_nitrogen(H_FACET) + &
                           ( Parms(iunit)%fraction_aground_npp_to_seeds(H_FACET) * amt )                 ! Translated from CSHED call and complexity of CSRSNK and BGLCIS.  Not sure if interpretted correctly.
         amt = uptake(N_STORE) * euf(BELOW)
         Rng(icell)%stored_nitrogen(H_FACET) = Rng(icell)%stored_nitrogen(H_FACET) - amt
         Rng(icell)%fine_root_nitrogen(H_FACET) = Rng(icell)%fine_root_nitrogen(H_FACET) + amt
       end if
         
       ! Takeup nutrients from the soil.  
       do ilayer = 1, 2                                                       ! Herbs are taking nutrients from the top two layers
         if (Rng(icell)%mineral_nitrogen(ilayer) .gt. tolerance) then
           fsol = 1.0
           if (available_nitrogen .gt. 0.00001) then
             calcup = uptake(N_SOIL) * Rng(icell)%mineral_nitrogen(ilayer) * fsol / available_nitrogen
           else
             calcup = 0.0
           end if
           amt = uptake(N_SOIL) * euf(ABOVE)
           Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) - amt
           Rng(icell)%leaf_nitrogen(H_FACET) = Rng(icell)%leaf_nitrogen(H_FACET) + &
                           ( ( 1.0 -  Parms(iunit)%fraction_aground_npp_to_seeds(H_FACET) ) * amt )  
           Rng(icell)%seed_nitrogen(H_FACET) = Rng(icell)%seed_nitrogen(H_FACET) + &
                           ( Parms(iunit)%fraction_aground_npp_to_seeds(H_FACET) * amt )                 ! Translated from CSHED call and complexity of CSRSNK and BGLCIS.  Not sure if interpretted correctly.
           amt = uptake(N_SOIL) * euf(BELOW)
           Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) - amt
           Rng(icell)%fine_root_nitrogen(H_FACET) = Rng(icell)%fine_root_nitrogen(H_FACET) + amt
         end if
       end do
       ! Takeup nutrients from nitrogen fixation
       if (Rng(icell)%plant_nitrogen_fixed(H_FACET) .gt. 0.0) then
         amt = uptake(N_FIX) * euf(ABOVE)
         Rng(icell)%nitrogen_source_sink = Rng(icell)%nitrogen_source_sink - amt
         Rng(icell)%leaf_nitrogen(H_FACET) = Rng(icell)%leaf_nitrogen(H_FACET) + &
                           ( ( 1.0 -  Parms(iunit)%fraction_aground_npp_to_seeds(H_FACET) ) * amt )  
         Rng(icell)%seed_nitrogen(H_FACET) = Rng(icell)%seed_nitrogen(H_FACET) + &
                           ( Parms(iunit)%fraction_aground_npp_to_seeds(H_FACET) * amt )                 ! Translated from CSHED call and complexity of CSRSNK and BGLCIS.  Not sure if interpretted correctly.
         amt = uptake(N_FIX) * euf(BELOW)
         Rng(icell)%nitrogen_source_sink = Rng(icell)%nitrogen_source_sink - amt
         Rng(icell)%fine_root_nitrogen(H_FACET) = Rng(icell)%fine_root_nitrogen(H_FACET) + amt
       end if
     end if
     ! Update lignin in plant parts incorporating new carbon contributions
     Rng(icell)%lignin_leaf(H_FACET) = Rng(icell)%leaf_carbon(H_FACET) * &
                                       Rng(icell)%plant_lignin_fraction(H_FACET, SURFACE_INDEX)
     Rng(icell)%lignin_fine_root(H_FACET) = Rng(icell)%fine_root_carbon(H_FACET) * &
                                            Rng(icell)%plant_lignin_fraction(H_FACET, SOIL_INDEX)
   else                               ! else no production this month
     Rng(icell)%total_pot_production(H_LYR) = 0.0
     Rng(icell)%total_pot_prod_limited_by_n(H_LYR) = 0.0
     Rng(icell)%total_pot_production(H_S_LYR) = 0.0
     Rng(icell)%total_pot_prod_limited_by_n(H_S_LYR) = 0.0
     Rng(icell)%total_pot_production(H_T_LYR) = 0.0
     Rng(icell)%total_pot_prod_limited_by_n(H_T_LYR) = 0.0
   end if                               ! If it is too cold or plants are dormant

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'HRB_GRW')

end subroutine


subroutine Woody_Growth (icell)
!**** Woody growth builds from potential production estimates to calculate actual woody growth, as limited by nutrients
!****
!**** R. Boone                 Last modified:  May 29, 2014
!                              A small change from S_LAYER to S_T_LAYER in EUC calculation.
   use Parameter_Vars
   use Structures
   implicit none

   real uptake(3), available_nitrogen, rimpct, cprod_left
   real Leaf_Allocation, line
   real gnfrac, tm, resp_temp_effect, euf(5), amt, fsol, calcup
   real rem_c_frac, tot_cup, mfprd(5)
   real cfrac(WOODY_PARTS)
   real shrub_c_sum, tree_c_sum, avg_pot_production
   real fmrspflux(5), site_potential
   integer icell, iunit, ilayer, ifacet, ipart, i

   iunit = Rng(icell)%range_type

   ! Determine nitrogen available for growth ... Woody plants can draw from all four layers
   available_nitrogen = 0.0
   do ilayer = 1, SOIL_LAYERS              
     available_nitrogen = available_nitrogen + Rng(icell)%mineral_nitrogen(ilayer)
   end do

   ! Century's site potential
   if (Globe(Rng(icell)%x, Rng(icell)%y)%precip .lt. 20.0/6.0) then
      site_potential = 1500.0
   else if (Globe(Rng(icell)%x, Rng(icell)%y)%precip .gt. 90.0/6.0) then
      site_potential = 3250.0
   else
      site_potential = line(Globe(Rng(icell)%x, Rng(icell)%y)%precip, 20.0/6.0, 1500.0, 90.0/6.0, 3250.0)
   end if
   site_potential = site_potential * Parms(iunit)%tree_site_potential   

   ! Century includes a downward correction of available minerals for savanna trees.  A module that was in Growth, but shifted here because in Growth in Century it was specific to trees.
   tm = min(available_nitrogen, 1.5)
   gnfrac = exp(-1.664 * exp(-0.00102 * tm * site_potential) * &
            Parms(iunit)%tree_basal_area_to_grass_nitrogen * Rng(icell)%tree_basal_area)
   if (gnfrac .lt. 0.0 .or. gnfrac .gt. 1.0) gnfrac = 0.0
   available_nitrogen = available_nitrogen * gnfrac
   
   ! Century includes flags set in the schedular to turn growth on or off.  I don't want to do that.  
   ! I wish to use degree-days, or topsoil available water to pet ratio, or temperature limits
   ! Phenology and proportion deciduous are now used to account for no growth by scenscent trees.  Could use daylength.   Right now, heat accumulation, presumably more sensitive to climate change
   if (Rng(icell)%ratio_water_pet > 0.0 .and. Globe(Rng(icell)%x, Rng(icell)%y)%min_temp .gt. 0.0 ) then
     do ifacet = S_FACET, T_FACET
       ! Calculate actual production values, and impact of root biomass on available nitrogen
       if ((Parms(iunit)%root_intercept_on_nutrients * Rng(icell)%fine_root_carbon(ifacet) * 2.5) .gt. 33) then
         rimpct = 1.0
       else
         rimpct = (1.0 - Parms(iunit)%root_intercept_on_nutrients * &
                 exp(-Parms(iunit)%root_effect_on_nutrients * Rng(icell)%leaf_carbon(ifacet) * 2.5))
       end if
       ! Determine actual production, restricted based on carbon to nitrogen ratios
       ! The following uses TREE_CFRAC, which is modified through dynamic carbon allocation (TREEDYNC). I would like to skip that, so using a method like in Growth.
       ! I am going to fill TREE_CFRAC with the values that come from initial distributions.  They will be static.  Those will be done in Each_Year, and may be updated in a given year if appropriate.
       ! This TREE_CFAC that is static may be inappropriate for deciduous trees.  May need two sets or make it dynamic.
       do i = 1, WOODY_PARTS
         cfrac(i) = Rng(icell)%carbon_allocation(ifacet, i)
       end do
       call Restrict_Production (icell, ifacet, 2, available_nitrogen, rimpct, cfrac, uptake)           ! 2 is correct here, just looking at leaves and fine roots.
       ! If the deciduous plants are in scenescence then decrease production by the proportion that are deciduous.  Those trees will not be growing.
       if (Rng(icell)%phenology(ifacet) .ge. 3.95) then                     ! Comparison to 4.0 exactly may be causing an error.
         if (ifacet .eq. S_FACET) then
           Rng(icell)%total_pot_production(S_LYR) = Rng(icell)%total_pot_production(S_LYR) * & 
                                                    ( 1.0 - Rng(icell)%prop_annual_decid(ifacet) )
           Rng(icell)%total_pot_production(S_T_LYR) = Rng(icell)%total_pot_production(S_T_LYR) * &
                                                    ( 1.0 - Rng(icell)%prop_annual_decid(ifacet) )
         else
           Rng(icell)%total_pot_production(T_LYR) = Rng(icell)%total_pot_production(T_LYR) * &
                                                    ( 1.0 - Rng(icell)%prop_annual_decid(ifacet) )
         end if
       end if
     end do     
   else
     Rng(icell)%total_pot_production(S_LYR) = 0.0
     Rng(icell)%total_pot_production(S_T_LYR) = 0.0
     Rng(icell)%total_pot_production(T_LYR) = 0.0
   end if

   do ifacet = S_FACET, T_FACET
     select case (ifacet)
       case (S_FACET)
          avg_pot_production = ( Rng(icell)%total_pot_production(S_LYR) + Rng(icell)%total_pot_production(S_T_LYR) ) / 2.0
       case (T_FACET)
          avg_pot_production = Rng(icell)%total_pot_production(T_LYR)
     end select
                          
     ! If growth occurs ...
     if (avg_pot_production .gt. 0.0) then
       ! Compute carbon allocation fraction for each woody part.
       ! ** RECALL that production is in biomass units ** but the bulk is proportions allocated, so units are ok.
       ! Portion left after some is taken by leaves and fine roots.  These get priority.  This doesn't shift material, just does the preliminary calculations
       Rng(icell)%carbon_allocation(ifacet, FINE_ROOT_INDEX) = Parms(iunit)%fraction_carbon_to_roots(ifacet)
       cprod_left = avg_pot_production - ( avg_pot_production * Rng(icell)%carbon_allocation(ifacet, FINE_ROOT_INDEX))
       Rng(icell)%carbon_allocation(ifacet, LEAF_INDEX) = &
                Leaf_Allocation(icell, ifacet, cprod_left, avg_pot_production)
       rem_c_frac = 1.0 - Rng(icell)%carbon_allocation(ifacet, FINE_ROOT_INDEX) - Rng(icell)%carbon_allocation(ifacet, LEAF_INDEX)
       if (rem_c_frac .lt. 1.0E-05) then
         do ipart = FINE_BRANCH_INDEX, COARSE_ROOT_INDEX                     ! No carbon left, so the remaining parts get 0 new growth.
           Rng(icell)%carbon_allocation(ifacet, ipart) = 0.0
         end do
       else
         ! A change from Century ... I don't want to include 10 more parameters controlling (juvenile and mature) carbon allocation to tree parts.
         ! I am going to use the initial carbon allocation.      
         shrub_c_sum = 0.0
         tree_c_sum = 0.0
         do ipart = 1, WOODY_PARTS
           shrub_c_sum = shrub_c_sum + Rng(icell)%shrub_carbon(ipart)
           tree_c_sum = tree_c_sum + Rng(icell)%tree_carbon(ipart)
         end do
         select case (ifacet)
           case (S_FACET)
             ! Shrubs
             if (shrub_c_sum .gt. 0.0) then                  ! If a division by zero error would occur, just don't change carbon_allocation
               Rng(icell)%carbon_allocation(S_FACET, FINE_BRANCH_INDEX) = Rng(icell)%shrub_carbon(FINE_BRANCH_INDEX) / shrub_c_sum
               Rng(icell)%carbon_allocation(S_FACET, COARSE_BRANCH_INDEX) = &
                                                                      Rng(icell)%shrub_carbon(COARSE_BRANCH_INDEX) / shrub_c_sum
               Rng(icell)%carbon_allocation(S_FACET, COARSE_ROOT_INDEX) = Rng(icell)%shrub_carbon(COARSE_ROOT_INDEX) / shrub_c_sum
             end if
           case (T_FACET)
             ! Trees
             if (tree_c_sum .gt. 0.0) then                  ! If a division by zero error would occur, just don't change carbon_allocation
               Rng(icell)%carbon_allocation(T_FACET, FINE_BRANCH_INDEX) = Rng(icell)%tree_carbon(FINE_BRANCH_INDEX) / tree_c_sum
               Rng(icell)%carbon_allocation(T_FACET, COARSE_BRANCH_INDEX) = &
                                                                      Rng(icell)%tree_carbon(COARSE_BRANCH_INDEX) / tree_c_sum         
               Rng(icell)%carbon_allocation(T_FACET, COARSE_ROOT_INDEX) = Rng(icell)%tree_carbon(COARSE_ROOT_INDEX) / tree_c_sum
             end if
         end select
         tot_cup = Rng(icell)%carbon_allocation(ifacet, FINE_BRANCH_INDEX) + &
                   Rng(icell)%carbon_allocation(ifacet, COARSE_BRANCH_INDEX) + &
                   Rng(icell)%carbon_allocation(ifacet, COARSE_ROOT_INDEX) 
         if (tot_cup .gt. 0.0) then
           do ipart = FINE_BRANCH_INDEX, COARSE_ROOT_INDEX
             Rng(icell)%carbon_allocation(ifacet, ipart) = Rng(icell)%carbon_allocation(ifacet, ipart) / tot_cup * rem_c_frac
           end do
         end if
       end if
     
       ! Calculate actual production values, and impact of root biomass on available nitrogen
       if ((Parms(iunit)%root_intercept_on_nutrients * Rng(icell)%fine_root_carbon(ifacet) * 2.5) .gt. 33) then
         rimpct = 1.0
       else
         rimpct = (1.0 - Parms(iunit)%root_intercept_on_nutrients * &
                   exp(-Parms(iunit)%root_effect_on_nutrients * Rng(icell)%fine_root_carbon(ifacet) * 2.5))
         end if
       ! Determine actual production, restricted based on carbon to nitrogen ratios
       ! The following uses TREE_CFRAC, which is modified through dynamic carbon allocation (TREEDYNC). I would like to skip that, so using a method like in Growth.
       ! I am going to fill TREE_CFRAC with the values that come from initial distributions.  They will be static.  Those will be done in Each_Year, and may be updated in a given year if appropriate.
       ! This TREE_CFAC that is static may be inappropriate for deciduous trees.  May need two sets or make it dynamic.
       do i = 1, WOODY_PARTS
         cfrac(i) = Rng(icell)%carbon_allocation(ifacet, i)
       end do
       call Restrict_Production (icell, ifacet, WOODY_PARTS, available_nitrogen, rimpct, cfrac, uptake)

       ! Calculate symbiotic N fixation accumulation 
       Rng(icell)%nitrogen_fixed(ifacet) = Rng(icell)%nitrogen_fixed(ifacet) + Rng(icell)%plant_nitrogen_fixed(ifacet)
       ! Accumulator was skipped ... too many and potentially confusing, so added later.  NFIXAC, EUPACC, EUPRT, TCNPRO
     
       ! Calculate production for each tree part
       ! This section deals with maintenance respiration calculations
       do ipart = 1, WOODY_PARTS
         mfprd(ipart) = Rng(icell)%carbon_allocation(ifacet, ipart) * avg_pot_production
         Rng(icell)%respiration_flows(ifacet) = Rng(icell)%respiration_flows(ifacet) + mfprd(ipart) + &
                                                Parms(iunit)%fraction_npp_to_respiration(ifacet)
       end do
       ! Growth of forest parts, with carbon added to the part and removed from the source-sink
       do ipart = 1, WOODY_PARTS
         select case (ifacet)
           case (S_FACET)
             Rng(icell)%shrub_carbon(ipart) = Rng(icell)%shrub_carbon(ipart) + mfprd(ipart)                               ! Translated from CSCHED ... double-check as needed
           case (T_FACET)
             Rng(icell)%tree_carbon(ipart) = Rng(icell)%tree_carbon(ipart) + mfprd(ipart)                               ! Translated from CSCHED ... double-check as needed
         end select
         Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink - mfprd(ipart)
       end do   
       Rng(icell)%fine_root_carbon(ifacet) = Rng(icell)%fine_root_carbon(ifacet) + mfprd(FINE_ROOT_INDEX)
       Rng(icell)%leaf_carbon(ifacet) = Rng(icell)%leaf_carbon(ifacet) + &
                           ( ( 1.0 -  Parms(iunit)%fraction_aground_npp_to_seeds(ifacet) ) * mfprd(LEAF_INDEX) )        ! Translated from CSHED call and complexity of CSRSNK and BGLCIS.  Not sure if interpretted correctly.
       Rng(icell)%seed_carbon(ifacet) = Rng(icell)%seed_carbon(ifacet) + &
                           ( Parms(iunit)%fraction_aground_npp_to_seeds(ifacet) * mfprd(LEAF_INDEX) )                   ! Translated from CSHED call and complexity of CSRSNK and BGLCIS.  Not sure if interpretted correctly.
       Rng(icell)%fine_branch_carbon(ifacet) = Rng(icell)%fine_branch_carbon(ifacet) + mfprd(FINE_BRANCH_INDEX)
       Rng(icell)%coarse_branch_carbon(ifacet) = Rng(icell)%coarse_branch_carbon(ifacet) + mfprd(COARSE_BRANCH_INDEX)
       Rng(icell)%coarse_root_carbon(ifacet) = Rng(icell)%coarse_root_carbon(ifacet) + mfprd(COARSE_ROOT_INDEX)

       ! Add maintenance respiration flow
       Rng(icell)%maintain_respiration(ifacet) = Rng(icell)%maintain_respiration(ifacet) + Rng(icell)%respiration_flows(ifacet)
       Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink - Rng(icell)%respiration_flows(ifacet)

       ! Maintenance respiration fluxes reduce maintenance respiration storage pool
       resp_temp_effect = 0.1 * exp(0.07 * Globe(Rng(icell)%x,Rng(icell)%y)%temperature_average)
       resp_temp_effect = min(1.0, resp_temp_effect)
       resp_temp_effect = max(0.0, resp_temp_effect)
       do ipart = 1, WOODY_PARTS                                               ! Merged two looping structures, the first didn't exist in Century, but works with this logic.
         select case (ifacet)
           case (S_FACET)
             fmrspflux(ipart) = Parms(iunit)%woody_max_fraction_npp_to_respiration(ipart) * resp_temp_effect * &
                                Rng(icell)%shrub_carbon(ipart)
           case (T_FACET)
             fmrspflux(ipart) = Parms(iunit)%woody_max_fraction_npp_to_respiration(ipart) * resp_temp_effect * &
                                Rng(icell)%tree_carbon(ipart)
         end select
         Rng(icell)%respiration_annual(ifacet) = Rng(icell)%respiration_annual(ifacet) + fmrspflux(ipart)
         Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink - fmrspflux(ipart)
       end do
    
       ! Actual uptake ... using the average of the vegetation layer parts
       select case (ifacet)
         case (S_FACET)
           if ( (Rng(icell)%total_pot_prod_limited_by_n(S_LYR) + &
                 Rng(icell)%total_pot_prod_limited_by_n(S_T_LYR) / 2.0 ) .gt. 0.0) then         
             do ipart = 1, WOODY_PARTS
               euf(ipart) = Rng(icell)%e_up(ifacet, ipart) / ( ( Rng(icell)%total_pot_prod_limited_by_n(S_LYR) + &
                                                                 Rng(icell)%total_pot_prod_limited_by_n(S_T_LYR) ) / 2.0 )
             end do
           else
             do ipart = 1, WOODY_PARTS
               euf(ipart) = 1.0                                                 ! Will cause no change later in the logic.
             end do
           end if
         case (T_FACET)
           if ( Rng(icell)%total_pot_prod_limited_by_n(T_LYR) .gt. 0.0) then         
             do ipart = 1, WOODY_PARTS
               euf(ipart) = Rng(icell)%e_up(ifacet, ipart) / Rng(icell)%total_pot_prod_limited_by_n(T_LYR) 
             end do
           else
             do ipart = 1, WOODY_PARTS
               euf(ipart) = 1.0                                                 ! Will cause no change later in the logic.
             end do
           end if           
       end select

       ! Takeup nutrients from internal storage pool.  
       if (Rng(icell)%stored_nitrogen(ifacet) .gt. 0.0) then
         do ipart = 1, WOODY_PARTS
           amt = uptake(N_STORE) * euf(ipart)
           select case (ifacet)
             case (S_FACET)
               Rng(icell)%shrub_nitrogen(ipart) = Rng(icell)%shrub_nitrogen(ipart) + amt
             case (T_FACET)
               Rng(icell)%tree_nitrogen(ipart) = Rng(icell)%tree_nitrogen(ipart) + amt
           end select               
           Rng(icell)%stored_nitrogen(ifacet) = Rng(icell)%stored_nitrogen(ifacet) - amt
         end do
       end if
       Rng(icell)%fine_root_nitrogen(ifacet) = Rng(icell)%fine_root_nitrogen(ifacet) + ( uptake(N_STORE) * euf(FINE_ROOT_INDEX) )
       Rng(icell)%leaf_nitrogen(ifacet) = Rng(icell)%leaf_nitrogen(ifacet) + &
             ( 1.0 -  Parms(iunit)%fraction_aground_npp_to_seeds(ifacet) ) * ( uptake(N_STORE) * euf(LEAF_INDEX) )        
       Rng(icell)%seed_nitrogen(ifacet) = Rng(icell)%seed_nitrogen(ifacet) + &
             ( Parms(iunit)%fraction_aground_npp_to_seeds(ifacet) ) * ( uptake(N_STORE) * euf(LEAF_INDEX) )
       Rng(icell)%fine_branch_nitrogen(ifacet) = Rng(icell)%fine_branch_nitrogen(ifacet) + ( uptake(N_STORE) * &
                                                    euf(FINE_BRANCH_INDEX) )
       Rng(icell)%coarse_branch_nitrogen(ifacet) = Rng(icell)%coarse_branch_nitrogen(ifacet) + ( uptake(N_STORE) * &
                                                    euf(COARSE_BRANCH_INDEX) )
       Rng(icell)%coarse_root_nitrogen(ifacet) = Rng(icell)%coarse_root_nitrogen(ifacet) + ( uptake(N_STORE) * &
                                                    euf(COARSE_ROOT_INDEX) )

       ! Takeup nutrients from the soil. ... Woody plants take nutrients from all four layers
       do ilayer = 1, SOIL_LAYERS
         if (Rng(icell)%mineral_nitrogen(ilayer) .gt. 0.00001) then
           fsol = 1.0
           if (available_nitrogen .gt. 0.00001) then
             calcup = uptake(N_SOIL) * Rng(icell)%mineral_nitrogen(ilayer) * fsol / available_nitrogen
           else
             calcup = 0.0 
           end if
           do ipart = 1, WOODY_PARTS
             amt = calcup * euf(ipart)
             select case (ifacet)
               case (S_FACET)
                 Rng(icell)%shrub_nitrogen(ipart) = Rng(icell)%shrub_nitrogen(ipart) + amt
               case (T_FACET)
                 Rng(icell)%tree_nitrogen(ipart) = Rng(icell)%tree_nitrogen(ipart) + amt
             end select
             Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) - amt
           end do
           Rng(icell)%fine_root_nitrogen(ifacet) = Rng(icell)%fine_root_nitrogen(ifacet) + ( calcup * euf(FINE_ROOT_INDEX) )
           Rng(icell)%leaf_nitrogen(ifacet) = Rng(icell)%leaf_nitrogen(ifacet) + ( calcup * euf(LEAF_INDEX) )
           Rng(icell)%fine_branch_nitrogen(ifacet) = Rng(icell)%fine_branch_nitrogen(ifacet) + ( calcup * &
                                                    euf(FINE_BRANCH_INDEX) )
           Rng(icell)%coarse_branch_nitrogen(ifacet) = Rng(icell)%coarse_branch_nitrogen(ifacet) + ( calcup * &
                                                    euf(COARSE_BRANCH_INDEX) )
           Rng(icell)%coarse_root_nitrogen(ifacet) = Rng(icell)%coarse_root_nitrogen(ifacet) + ( calcup * &
                                                    euf(COARSE_ROOT_INDEX) )  
         end if
       end do

       ! Takeup nutrients from nitrogen fixation.
       if (Rng(icell)%plant_nitrogen_fixed(ifacet) .gt. 0.000001) then
         do ipart = 1, WOODY_PARTS
           amt = uptake(N_FIX) * euf(ipart)
           select case (ifacet)
             case (S_FACET)
               Rng(icell)%shrub_nitrogen(ipart) = Rng(icell)%shrub_nitrogen(ipart) + amt
             case (T_FACET)
               Rng(icell)%tree_nitrogen(ipart) = Rng(icell)%tree_nitrogen(ipart) + amt
           end select
           Rng(icell)%nitrogen_source_sink = Rng(icell)%nitrogen_source_sink - amt
         end do
         Rng(icell)%fine_root_nitrogen(ifacet) = Rng(icell)%fine_root_nitrogen(ifacet) + ( uptake(N_STORE) * euf(FINE_ROOT_INDEX) )
         Rng(icell)%leaf_nitrogen(ifacet) = Rng(icell)%leaf_nitrogen(ifacet) + ( uptake(N_STORE) * euf(LEAF_INDEX) )
         Rng(icell)%fine_branch_nitrogen(ifacet) = Rng(icell)%fine_branch_nitrogen(ifacet) + ( uptake(N_STORE) * &
                                                    euf(FINE_BRANCH_INDEX) )
         Rng(icell)%coarse_branch_nitrogen(ifacet) = Rng(icell)%coarse_branch_nitrogen(ifacet) + ( uptake(N_STORE) * &
                                                    euf(COARSE_BRANCH_INDEX) )
         Rng(icell)%coarse_root_nitrogen(ifacet) = Rng(icell)%coarse_root_nitrogen(ifacet) + ( uptake(N_STORE) * &
                                                    euf(COARSE_ROOT_INDEX) )
       end if
       ! Update lignin in plant parts incorporating new carbon contributions
       Rng(icell)%lignin_leaf(ifacet) = Rng(icell)%leaf_carbon(ifacet) * &
                                        Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX)
       Rng(icell)%lignin_fine_root(ifacet) = Rng(icell)%fine_root_carbon(ifacet) * &
                                             Rng(icell)%plant_lignin_fraction(ifacet, SOIL_INDEX)
       Rng(icell)%lignin_fine_branch(ifacet) = Rng(icell)%fine_branch_carbon(ifacet) * &
                                               Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX)
       Rng(icell)%lignin_coarse_branch(ifacet) = Rng(icell)%coarse_branch_carbon(ifacet) * &
                                                 Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX)
       Rng(icell)%lignin_coarse_root(ifacet) = Rng(icell)%coarse_root_carbon(ifacet) * &
                                               Rng(icell)%plant_lignin_fraction(ifacet, SOIL_INDEX)
     else
     ! There is no production this month
       select case (ifacet)
         case (S_FACET)  
           Rng(icell)%total_pot_prod_limited_by_n(S_LYR) = 0.0
           Rng(icell)%total_pot_production(S_LYR) = 0.0
           Rng(icell)%total_pot_prod_limited_by_n(S_T_LYR) = 0.0
           Rng(icell)%total_pot_production(S_T_LYR) = 0.0
         case (T_FACET)
           Rng(icell)%total_pot_prod_limited_by_n(T_LYR) = 0.0
           Rng(icell)%total_pot_production(T_LYR) = 0.0     
       end select
       do ipart = 1, WOODY_PARTS
         Rng(icell)%e_up(ifacet, ipart) = 0.0
       end do 
     end if
  
   end do
  
   if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'WDY_GRW')
  
end subroutine



subroutine Restrict_Production (icell, ifacet, nparts, available_nitrogen, rimpct, cfrac, uptake)
!**** Restrict actual production for plants based on carbon to nitrogen ratios
!**** Note, currently (and it could be changed for clarity), cfrac stores carbon allocation in ABOVE and BELOW for herbs,
!**** and in plant parts for woody plants.  
!****
!**** R. Boone                 Last modified:  February 16, 2013
   use Parameter_Vars
   use Structures
   implicit none

   real ctob, available_nitrogen, n_available, rimpct, cfrac(WOODY_PARTS), max_n, uptake(3)    ! min_n
   real min_n_ci(WOODY_PARTS), max_n_ci(WOODY_PARTS), ustorg, temp_prod
   integer icell, ifacet, nparts, iunit, ipart

   iunit = Rng(icell)%range_type

   if ((available_nitrogen .le. 1E-4) .and. (Parms(iunit)%max_symbiotic_n_fixation_ratio .eq. 0.0)) return            ! Won't things go unset?  Set something to zero?

   uptake(N_STORE) = 0.0
   uptake(N_SOIL) = 0.0
   uptake(N_FIX) = 0.0

   ! Calculate available nitrogen based on maximum fraction and impact of root biomass
   n_available = (available_nitrogen * Parms(iunit)%fraction_nitrogen_available * rimpct) + Rng(icell)%stored_nitrogen(ifacet)

   ! Compute weighted average carbon to biomass conversion factor.  
   ! The structure is a little odd here, because this section can be called for grasses or trees.  
   ! ctob is a weighted average carbon to biomass conversion factor  
   ! I don't trust the structure, replacing with a simplier structure (but the other now appears correct as well)
   if (ifacet .eq. H_FACET) then
     ctob = 2.5                                                             ! The same conversion is used above and below ground
   else
     ctob = ( cfrac(LEAF_INDEX) * 2.5 ) + ( cfrac(FINE_ROOT_INDEX) * 2.5 ) + ( cfrac(FINE_BRANCH_INDEX) * 2.0 ) + &
            ( cfrac(COARSE_BRANCH_INDEX) * 2.0 ) + ( cfrac(COARSE_ROOT_INDEX) * 2.0 )       ! Note conversion from carbon to biomass for wood parts is 2.0, rather than 2.5
   end if
   
   ! Calculate average N/C of whole plant (grass, shrub, or tree)
   max_n = 0.0   
   do ipart = 1, nparts
     min_n_ci(ipart) = 1.0 / Parms(iunit)%maximum_c_n_ratio(ifacet, ipart)     ! CHECK THIS.  IS IT IN ERROR?  Unusual formatting in CENTURY
     max_n_ci(ipart) = 1.0 / Parms(iunit)%minimum_c_n_ratio(ifacet, ipart)     ! Note that indicators maximum and minimum are flipped
   end do
   ! The following bases results on shoots and roots only.   It works for herbs and woody as well given that the two values of interest are 1 and 2 regardless.
   max_n = max_n + (cfrac(FINE_ROOT_INDEX) * max_n_ci(FINE_ROOT_INDEX))        ! CENTURY includes a biomass to carbon convertion on CFRAC(FROOT) * MAXECI.  I'm not sure why.   I've removed them for now.
   max_n = max_n + ((1.0 - cfrac(FINE_ROOT_INDEX)) * max_n_ci(LEAF_INDEX))       ! CENTURY includes a biomass to carbon convertion on CFRAC(FROOT) * MAXECI.  I'm not sure why.   I've removed them for now.
   ! Calculate average nutrient content
!   max_n = max_n * ctob      ! Skipping this for now (counter to CENTURY RSTRP.F).  So MAX_N is being passed as a weighted nitrogen concentration.   Converting what I have to biomass doesn't make any sense.

   ! Compute the limitation on nutrients.  Min_N need not be passed.  It used in Century only for automatic fertilization.
   call Nutrient_Limitation (icell, ifacet, nparts, max_n, min_n_ci, max_n_ci, cfrac, ctob, n_available)

   ! Calculate relative yield skipped for now.  Not used in module.

   ! Calculate uptake from all sources (storage, soil, plant n fixed)
   select case (ifacet)                                                   ! Need to use the average of production.  
     case (H_FACET)
       temp_prod = ( Rng(icell)%total_pot_prod_limited_by_n(H_LYR) + Rng(icell)%total_pot_prod_limited_by_n(H_S_LYR) + &
                     Rng(icell)%total_pot_prod_limited_by_n(H_T_LYR) ) / 3.0
     case (S_FACET)
       temp_prod = ( Rng(icell)%total_pot_prod_limited_by_n(S_LYR) + Rng(icell)%total_pot_prod_limited_by_n(S_T_LYR) ) / 2.0
     case (T_FACET)
       temp_prod = Rng(icell)%total_pot_prod_limited_by_n(T_LYR)
   end select
   ustorg = min(Rng(icell)%stored_nitrogen(ifacet), temp_prod)
   ! If the storage pool contains all that is needed for uptake, then ...
   if (temp_prod .le. ustorg) then
     uptake(N_STORE) = temp_prod
     uptake(N_SOIL) = 0.0
   ! Otherwise take what is needed from the storage pool          (unneeded elseif in Century, skipped)
   else 
     uptake(N_STORE) = Rng(icell)%stored_nitrogen(ifacet)
     uptake(N_SOIL) = temp_prod - Rng(icell)%stored_nitrogen(ifacet) - Rng(icell)%plant_nitrogen_fixed(ifacet)
   end if
   uptake(N_FIX) = Rng(icell)%plant_nitrogen_fixed(ifacet)

   if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'RESTRCT')

end subroutine


subroutine Nutrient_Limitation (icell, ifacet, nparts, max_n, min_n_ci, max_n_ci, cfrac, ctob, available_nitrogen)
!**** Compute nutrient limitation on growth  (NUTRLM.F in Century)
!****
!**** R. Boone                 Last modified:  July 16, 2014
   use Parameter_Vars
   use Structures
   implicit none
   
   real max_n, min_n_ci(WOODY_PARTS), max_n_ci(WOODY_PARTS), cfrac(WOODY_PARTS), ctob, available_nitrogen
   real lfrac(V_LYRS)
   real max_n_fix, demand, totaln, ecfor(WOODY_PARTS), cpbe, totaln_used
   real total_pot_production_biomass, total_pot_production_carbon
   integer icell, ifacet, nparts, ipart, iunit
   
   iunit = Rng(icell)%range_type
   ! A lengthy section within NUTRLM deals with P and S and was skipped

   ! Get the total production for the facet.  Include the different layers (i.e., do not take the average, as that will reduce total production by default)
   max_n_fix = 0.0
   select case (ifacet) 
     case (H_FACET)
       total_pot_production_biomass = Rng(icell)%total_pot_production(H_LYR) + Rng(icell)%total_pot_production(H_S_LYR) + &
                                      Rng(icell)%total_pot_production(H_T_LYR) 
     case (S_FACET)
       total_pot_production_biomass = Rng(icell)%total_pot_production(S_LYR) + Rng(icell)%total_pot_production(S_T_LYR) 
     case (T_FACET)
       total_pot_production_biomass = Rng(icell)%total_pot_production(T_LYR) 
   end select
   ! Convert to carbon
   total_pot_production_carbon = total_pot_production_biomass / ctob
   ! Demand based on the maximum nitrogen / carbon ratio
   ! Max_n is the maximum nitrogen to carbon ratio ... convert to carbon then multiply by N:C ratio.  so X:C yields demand
   demand = total_pot_production_carbon * max_n                            ! Weighted concentration of nitrogen in plant parts
   max_n_fix = Parms(iunit)%max_symbiotic_n_fixation_ratio * total_pot_production_carbon
   totaln = available_nitrogen + max_n_fix
   
   ! Calculation of a2drat(n) skipped.  It appears associated with dynamic carbon allocation, not implemented here, and not used in this module.

   ! New N/C ratios based on nitrogen available
   if (totaln .gt. demand) then
     do ipart = 1, nparts
       ecfor(ipart) = max_n_ci(ipart)                        ! Nitrogen is readily available, and the maximum concentration per part is appropriate
     end do
     totaln_used = demand                                         ! Setting a division equal to 1, essentially, used later
   else
     if (demand .eq. 0.0) then
!       write(*,*) 'Error in Nutrient_Limitation, demand = 0.0'        ! Disabling warning ... in a global model, a cell with no demand is reasonable, in the high arctic, say.
!       stop                             DEBUG
!       The following is part of the DEBUG, to avoid errors    ! DEBUG
       demand = 0.001                                                     ! DEBUG
     end if
     do ipart = 1, nparts
       ecfor(ipart) = min_n_ci(ipart) + ( ( max_n_ci(ipart) - min_n_ci(ipart) ) * (totaln / demand))           ! Nitrogen is limited, and so a fractional portion is assigned
       totaln_used = totaln
     end do
   end if
   
   cpbe = 0.0
   ctob = 0.0
   ! Total potential production with nutrient limitation.  Here CPBE remains a N to C ratio, but adjusted for demand.
   do ipart = 1, nparts
     if (ipart .lt. 3) then           
       cpbe = cpbe + ((cfrac(ipart) * ecfor(ipart)) / 2.5)                 ! Leaves and fine roots
       ctob = ctob + ( cfrac(ipart) * 2.5 )                                ! From RESTRP.F
     else
       cpbe = cpbe + ((cfrac(ipart) * ecfor(ipart)) / 2.0)                 ! Fine and coarse branches and coarse roots
       ctob = ctob + ( cfrac(ipart) * 2.0 )                                ! From RESTRP.F
     end if
   end do
   
   ! Increase the nitrogen estimate in line with an increase in carbon to biomass conversion (i.e., between 2 and 2.5 based on plant parts involved)
   ! Send the cpbe value from carbon to biomass   
   cpbe = cpbe * ctob
   if (cpbe .eq. 0.0) then
!     write(*,*) 'Error in Nutrient Limitation, CPBE = 0.0'    ! Disabling warning ... in a global model, a cell with no demand is reasonable, in the high arctic, say.
!       stop                             DEBUG
!       The following is part of the DEBUG, to avoid errors    ! DEBUG
     cpbe = 0.001
   end if
   
   ! Calculate potential production for the nutrient limitation
   cpbe = totaln_used / cpbe
   
   ! Automatic fertilization methods skipped.
   ! See if production is limited by nutrients  (works to compute limiting nutrient in Century, but just one here)
   if (total_pot_production_biomass .gt. cpbe) then
     total_pot_production_biomass = total_pot_production_biomass * ( totaln_used / demand )
   end if
   
   ! Adjustments considering P skipped.  P not considered here.
   ! Total potential production with nitrogen limitation
   ! First store nitrogen uptake per plant part.  Note use of carbon here, rather than biomass as in Century.   Carbon is being related to nitrogen here.
   do ipart = 1,nparts
     Rng(icell)%e_up(ifacet,ipart) = total_pot_production_carbon * cfrac(ipart) * ecfor(ipart)
     if (Rng(icell)%e_up(ifacet,ipart) .lt. 0.0) then
       Rng(icell)%e_up(ifacet,ipart) = 0.0
     end if
   end do
   ! Then put in limited production estimates
   select case (ifacet)
     case (H_FACET)                       
       ! Assuming production is in-line with existing biomass, so using a weighted average.
       ! Total_pot_production is calculated from first principles, so I will use that as a guide.  Using a little brute force, for speed
       ! Total_pot_production isn't limited by n, but assuming all the layers are equally affected by n limitation is appropriate.
       if ( ( Rng(icell)%total_pot_production(H_LYR) + Rng(icell)%total_pot_production(H_S_LYR) + &
            Rng(icell)%total_pot_production(H_T_LYR) ) .ge. 0.00001 ) then
          lfrac(H_LYR) = Rng(icell)%total_pot_production(H_LYR) / ( Rng(icell)%total_pot_production(H_LYR) + &
              Rng(icell)%total_pot_production(H_S_LYR) + Rng(icell)%total_pot_production(H_T_LYR) ) 
          lfrac(H_S_LYR) = Rng(icell)%total_pot_production(H_S_LYR) / ( Rng(icell)%total_pot_production(H_LYR) + &
              Rng(icell)%total_pot_production(H_S_LYR) + Rng(icell)%total_pot_production(H_T_LYR) ) 
          lfrac(H_T_LYR) = Rng(icell)%total_pot_production(H_T_LYR) / ( Rng(icell)%total_pot_production(H_LYR) + &
              Rng(icell)%total_pot_production(H_S_LYR) + Rng(icell)%total_pot_production(H_T_LYR) ) 
       else
          lfrac(H_LYR) = 0.0
          lfrac(H_S_LYR) = 0.0
          lfrac(H_T_LYR) = 0.0
       end if
       
       Rng(icell)%total_pot_prod_limited_by_n(H_LYR) = total_pot_production_biomass * lfrac(H_LYR)
       Rng(icell)%total_pot_prod_limited_by_n(H_S_LYR) = total_pot_production_biomass * lfrac(H_S_LYR)
       Rng(icell)%total_pot_prod_limited_by_n(H_T_LYR) = total_pot_production_biomass * lfrac(H_T_LYR)       
     case (S_FACET)
       if ( ( Rng(icell)%total_pot_production(S_LYR) + Rng(icell)%total_pot_production(S_T_LYR) ) .ge. 0.00001 ) then
         lfrac(S_LYR) = Rng(icell)%total_pot_production(S_LYR) / ( Rng(icell)%total_pot_production(S_LYR) + &
             Rng(icell)%total_pot_production(S_T_LYR) ) 
         lfrac(S_T_LYR) = Rng(icell)%total_pot_production(S_T_LYR) / ( Rng(icell)%total_pot_production(S_LYR) + &
             Rng(icell)%total_pot_production(S_T_LYR) ) 
       else
         lfrac(S_LYR) = 0.0
         lfrac(S_T_LYR) = 0.0
       end if

       Rng(icell)%total_pot_prod_limited_by_n(S_LYR) = total_pot_production_biomass * lfrac(S_LYR)
       Rng(icell)%total_pot_prod_limited_by_n(S_T_LYR) = total_pot_production_biomass * lfrac(S_T_LYR)
     case (T_FACET)
       Rng(icell)%total_pot_prod_limited_by_n(T_LYR) = total_pot_production_biomass
   end select  

   ! Compute nitrogen fixation that actually occurs  (Using average for plant nitrogen fixed in shrubs
   select case (ifacet)
     case (H_FACET)
       Rng(icell)%plant_nitrogen_fixed(ifacet) = max( totaln_used - available_nitrogen, 0.0)
     case (S_FACET)
       Rng(icell)%plant_nitrogen_fixed(ifacet) = max( totaln_used - available_nitrogen, 0.0)
     case (T_FACET)
       Rng(icell)%plant_nitrogen_fixed(ifacet) = max( totaln_used - available_nitrogen, 0.0) 
   end select

   if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'NUT_LIM')
   
end subroutine


real function Leaf_Allocation (icell, ifacet, cprod_left, cprodf)
!**** Compute optimum leaf area index, based on a maximum and available production
!****
!**** R. Boone                 Last modified:  October 5, 2013 ... changed * to + in optimum leaf area index calculation
   use Parameter_Vars
   use Structures
   implicit none

   real cprod_left, cprodf
   real leafprod, rleavc_opt
   integer icell, iunit, ifacet
   
   iunit = Rng(icell)%range_type
   
   if (Rng(icell)%coarse_branch_carbon(ifacet) .gt. 0.0 ) then
     Rng(icell)%optimum_leaf_area_index(ifacet) = Parms(iunit)%maximum_leaf_area_index * &
         ( Rng(icell)%coarse_branch_carbon(ifacet) * 2.0 ) / ( Parms(iunit)%k_leaf_area_index + &
         ( Rng(icell)%coarse_branch_carbon(ifacet) * 2.0 ) )
     if (Rng(icell)%optimum_leaf_area_index(ifacet) .lt. 0.1) then
       Rng(icell)%optimum_leaf_area_index(ifacet) = 0.1
     end if
   else
     Rng(icell)%optimum_leaf_area_index(ifacet) = 0.1
   end if

   rleavc_opt = Rng(icell)%optimum_leaf_area_index(ifacet) / (2.5 * Parms(iunit)%biomass_to_leaf_area_index_factor)
   if (rleavc_opt .gt. Rng(icell)%leaf_carbon(ifacet)) then
     leafprod = min((rleavc_opt - Rng(icell)%leaf_carbon(ifacet)), cprod_left)
   else
     leafprod = 0.0
   end if

   if (cprodf .gt. 0.0) then
     Leaf_Allocation = leafprod / cprodf
   else
     !!NOTE!! CENTURY has this as 0, which is likely appropriate once trees get larger and optimum lia is more reasonable.
     Leaf_Allocation = 0.01
   end if
   if (Leaf_Allocation .lt. 0.01) Leaf_Allocation = 0.01                                  ! This trims where Century throws errors.
   if (Leaf_Allocation .gt. 1.0) Leaf_Allocation = 1.0

   if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'LEAF_AL')

end function


subroutine Grazing (icell)
!**** Remove forage that is grazed (GREM in Century).
!****
!**** R. Boone                 Last modified:  February 4, 2011
   use Parameter_Vars
   use Structures
   implicit none
   
   real carbon_removed, nitrogen_removed, total_carbon, total_nitrogen, carbon_returned, nitrogen_returned
   real fraction_nitrogen_grazed_returned, urine, feces, avg_lignin
   integer icell, ifacet, iunit
   
   iunit = Rng(icell)%range_type
   
   total_carbon = 0.0
   total_nitrogen = 0.0
   do ifacet = 1, FACETS
     ! Shoots removed.  Moving away from Century somewhat.    
     if (Rng(icell)%leaf_carbon(ifacet) .gt. 0.0) then
       carbon_removed = Rng(icell)%leaf_carbon(ifacet) * ( Rng(icell)%fraction_live_removed_grazing * &
                                                         Parms(iunit)%fraction_grazed_by_facet(ifacet) )
     else
       carbon_removed = 0.0
     end if
     if (Rng(icell)%leaf_nitrogen(ifacet) .gt. 0.0) then
       nitrogen_removed = Rng(icell)%leaf_nitrogen(ifacet) * ( Rng(icell)%fraction_live_removed_grazing * &
                                                         Parms(iunit)%fraction_grazed_by_facet(ifacet) )
     else
       nitrogen_removed = 0.0
     end if
     Rng(icell)%leaf_carbon(ifacet) = Rng(icell)%leaf_carbon(ifacet) - carbon_removed
     Rng(icell)%leaf_nitrogen(ifacet) = Rng(icell)%leaf_nitrogen(ifacet) - nitrogen_removed
     Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink + carbon_removed
     Rng(icell)%nitrogen_source_sink = Rng(icell)%nitrogen_source_sink + nitrogen_removed
     total_carbon = total_carbon + carbon_removed
     total_nitrogen = total_nitrogen + nitrogen_removed
   
     ! Standing dead removed.  
     if (Rng(icell)%dead_standing_carbon(ifacet) .gt. 0.0) then
       carbon_removed = Rng(icell)%dead_standing_carbon(ifacet) * ( Rng(icell)%fraction_dead_removed_grazing * &
                                                         Parms(iunit)%fraction_grazed_by_facet(ifacet) )
     else
       carbon_removed = 0.0
     end if
     if (Rng(icell)%dead_standing_nitrogen(ifacet) .gt. 0.0) then
       nitrogen_removed = Rng(icell)%dead_standing_nitrogen(ifacet) * ( Rng(icell)%fraction_dead_removed_grazing * &
                                                         Parms(iunit)%fraction_grazed_by_facet(ifacet) )
     else
       nitrogen_removed = 0.0
     end if
     Rng(icell)%dead_standing_carbon(ifacet) = Rng(icell)%dead_standing_carbon(ifacet) - carbon_removed
     Rng(icell)%dead_standing_nitrogen(ifacet) = Rng(icell)%dead_standing_nitrogen(ifacet) - nitrogen_removed
     Rng(icell)%carbon_source_sink = Rng(icell)%carbon_source_sink + carbon_removed
     Rng(icell)%nitrogen_source_sink = Rng(icell)%nitrogen_source_sink + nitrogen_removed
     total_carbon = total_carbon + carbon_removed
     total_nitrogen = total_nitrogen + nitrogen_removed
   end do
   
   ! Return portions of the carbon and nitrogen grazed back to the environment.  Carbon in the form of feces, urine includes nitrogen.
   carbon_returned = Parms(iunit)%fraction_carbon_grazed_returned * total_carbon
   if (carbon_returned .le. 0.0) then
     carbon_returned = 0.0
     nitrogen_returned = 0.0
     fraction_nitrogen_grazed_returned = 0.0
   else
     ! The portion of nitrogen returned is a function of clay content  (CENTURY GREM.F)
     if (Rng(icell)%clay(SURFACE_INDEX) .lt. 0.0) then
       fraction_nitrogen_grazed_returned = 0.7
     else if (Rng(icell)%clay(SURFACE_INDEX) .gt. 0.3) then
       fraction_nitrogen_grazed_returned = 0.85
     else
       fraction_nitrogen_grazed_returned = (0.85 - 0.7) / (0.3 - 0.0) * (Rng(icell)%clay(SURFACE_INDEX) - 0.3) + 0.85
     end if
   end if
   nitrogen_returned = fraction_nitrogen_grazed_returned * total_nitrogen
   urine = ( 1. - Parms(iunit)%fraction_excreted_nitrogen_in_feces) * nitrogen_returned
   feces = Parms(iunit)%fraction_excreted_nitrogen_in_feces * nitrogen_returned
   
   ! Do flows
   Rng(icell)%mineral_nitrogen(SURFACE_INDEX) = Rng(icell)%mineral_nitrogen(SURFACE_INDEX) + urine
   Rng(icell)%nitrogen_source_sink = Rng(icell)%nitrogen_source_sink - urine
   
   Rng(icell)%volatized_n = Rng(icell)%volatized_n * ( Parms(iunit)%fraction_urine_volatized * urine )
   urine = urine - Rng(icell)%volatized_n
   
   ! Move materials into litter
   avg_lignin = ( Rng(icell)%lignin_leaf(H_FACET) + Rng(icell)%lignin_leaf(S_FACET) + &
                  Rng(icell)%lignin_leaf(T_FACET) ) / 3.0
   avg_lignin = max(0.02, avg_lignin)                                              ! From Century CmpLig.f
   avg_lignin = min(0.50, avg_lignin)                  
   call Partition_Litter (icell, SURFACE_INDEX, feces, urine, avg_lignin, 'grazing')  
   
   if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'GRAZING')

end subroutine