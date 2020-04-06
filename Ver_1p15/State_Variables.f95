subroutine Write_State
!**** Saves state variable results to a file name that is provided.
!****
!**** R. Boone   Last modified: December 28, 2013  
   use Parameter_Vars
   use Structures
   implicit none

   integer  i, j, icell

   open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//Sim_Parm%state_var_file_out, &
       STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)
     
   if (ioerr == 0) then
     ! Material will be written at the top of the state variable file.   The file will be in ASCII for now.
     ! It may have to be converted to binary if too large.  The heading material will provide information, but
     ! if it doesn't agree with what is read-in, there may not be an error provided.  That will allow for scenarios
     ! to use information that differs from the original.  So some information will provide errors (e.g., the 
     ! world dimension), and others won't.
     ! 
     ! Generally, entries in Each_Year that are reset to 0 are not stored here.  They would be zeroed-out regardless.
     write(SHORT_USE_FILE,*) 'x_dim                    : ', x_dim                                                      
     write(SHORT_USE_FILE,*) 'y_dim                    : ', y_dim                                                      
     write(SHORT_USE_FILE,*) 'lower_x                  : ', lower_x                                                    
     write(SHORT_USE_FILE,*) 'lower_y                  : ', lower_y                                                    
     write(SHORT_USE_FILE,*) 'upper_x                  : ', upper_x                                                    
     write(SHORT_USE_FILE,*) 'upper_y                  : ', upper_y                                                    
     write(SHORT_USE_FILE,*) 'cellsize                 : ', cellsize                                                   
     write(SHORT_USE_FILE,*) 'range_cells              : ', range_cells                                                
     write(SHORT_USE_FILE,*) 'bin_path                 : ', bin_path                                                   
     write(SHORT_USE_FILE,*) 'app_path                 : ', app_path                                                   
     write(SHORT_USE_FILE,*) 'parm_path                : ', parm_path                                                  
     write(SHORT_USE_FILE,*) 'out_path                 : ', out_path                                                   
     write(SHORT_USE_FILE,*) 'zone_map                 : ', Sim_Parm%zone_map              
     write(SHORT_USE_FILE,*) 'land_map                 : ', Sim_Parm%land_map                                          
     write(SHORT_USE_FILE,*) 'elev_map                 : ', Sim_Parm%elev_map                                          
     write(SHORT_USE_FILE,*) 'class_map                : ', Sim_Parm%class_map                                         
     write(SHORT_USE_FILE,*) 'class_legend             : ', Sim_Parm%class_legend                                      
     write(SHORT_USE_FILE,*) 'latitude_map             : ', Sim_Parm%latitude_map                                      
     write(SHORT_USE_FILE,*) 'landscape_type_map       : ', Sim_Parm%landscape_type_map                                
     write(SHORT_USE_FILE,*) 'top_sand_map             : ', Sim_Parm%top_sand_map                                      
     write(SHORT_USE_FILE,*) 'top_silt_map             : ', Sim_Parm%top_silt_map                                      
     write(SHORT_USE_FILE,*) 'top_clay_map             : ', Sim_Parm%top_clay_map                                      
     write(SHORT_USE_FILE,*) 'top_gravel_map           : ', Sim_Parm%top_gravel_map                                    
     write(SHORT_USE_FILE,*) 'top_bulk_density_map     : ', Sim_Parm%top_bulk_density_map                              
     write(SHORT_USE_FILE,*) 'top_organic_carbon_map   : ', Sim_Parm%top_organic_carbon_map                            
     write(SHORT_USE_FILE,*) 'sub_sand_map             : ', Sim_Parm%sub_sand_map                                      
     write(SHORT_USE_FILE,*) 'sub_silt_map             : ', Sim_Parm%sub_silt_map                                      
     write(SHORT_USE_FILE,*) 'sub_clay_map             : ', Sim_Parm%sub_clay_map                                      
     write(SHORT_USE_FILE,*) 'sub_gravel_map           : ', Sim_Parm%sub_gravel_map                                    
     write(SHORT_USE_FILE,*) 'sub_bulk_density_map     : ', Sim_Parm%sub_bulk_density_map                              
     write(SHORT_USE_FILE,*) 'sub_organic_carbon_map   : ', Sim_Parm%sub_organic_carbon_map                            
     write(SHORT_USE_FILE,*) 'deciduous_tree_cover_map : ', Sim_Parm%deciduous_tree_cover_map                          
     write(SHORT_USE_FILE,*) 'evergreen_tree_cover_map : ', Sim_Parm%evergreen_tree_cover_map                          
     write(SHORT_USE_FILE,*) 'shrub_cover_map          : ', Sim_Parm%shrub_cover_map                                   
     write(SHORT_USE_FILE,*) 'herb_cover_map           : ', Sim_Parm%herb_cover_map                                    
     write(SHORT_USE_FILE,*) 'precip_average_map       : ', Sim_Parm%precip_average_map                                
     write(SHORT_USE_FILE,*) 'temperature_average_map  : ', Sim_Parm%temperature_average_map                           
     write(SHORT_USE_FILE,*) 'precip_path_prefix       : ', Sim_Parm%precip_path_prefix                                
     write(SHORT_USE_FILE,*) 'max_temp_path_prefix     : ', Sim_Parm%max_temp_path_prefix                              
     write(SHORT_USE_FILE,*) 'min_temp_path_prefix     : ', Sim_Parm%min_temp_path_prefix                              
     write(SHORT_USE_FILE,*) 'precip_temp_suffix       : ', Sim_Parm%precip_temp_suffix                                
     write(SHORT_USE_FILE,*) 'parms_file_name          : ', Sim_Parm%parms_file_name                                   
     write(SHORT_USE_FILE,*) 'start_yr                 : ', Sim_Parm%start_yr                                          
     write(SHORT_USE_FILE,*) 'end_yr                   : ', Sim_Parm%end_yr                                            
     write(SHORT_USE_FILE,*) 'state_var_flag           : ', Sim_Parm%state_var_flag                                    
     write(SHORT_USE_FILE,*) 'state_var_file_out       : ', Sim_Parm%state_var_file_out                                
     write(SHORT_USE_FILE,*) 'state_var_file_in        : ', Sim_Parm%state_var_file_in                                 
                                                                                                                       
     do icell = 1, range_cells                                                                                         
       write(SHORT_USE_FILE,*) 'zone       : ', Globe(Rng(icell)%x, Rng(icell)%y)%zone                                 
       write(SHORT_USE_FILE,*) 'x          : ', Rng(icell)%x                                                           
       write(SHORT_USE_FILE,*) 'y          : ', Rng(icell)%y                                                           
       write(SHORT_USE_FILE,*) 'range_type : ', Rng(icell)%range_type                                                  
       write(SHORT_USE_FILE,*) Rng(icell)%last_month_day_length     ! The day length of the previous month, to know when spring and fall come.
       write(SHORT_USE_FILE,*) Rng(icell)%day_length_increasing     ! Increasing or decreasing day length, comparing the current to previous day lengths.
       write(SHORT_USE_FILE,*) Rng(icell)%day_length                ! Day length, calculated based on latitude and month
       write(SHORT_USE_FILE,*) Rng(icell)%heat_accumulation         ! Heat accumulation above a base temperature (e.g., 4.4 C in Boone (1999))
       write(SHORT_USE_FILE,*) (Rng(icell)%facet_cover(i),i=1,FACETS)       ! The proportion occupied by each facet
       write(SHORT_USE_FILE,*) (Rng(icell)%total_population(i),i=1,V_LYRS)  ! The total population of each vegetation layer
       write(SHORT_USE_FILE,*) Rng(icell)%bare_cover                ! Bare cover stored, rather than continually summing the three facets.
       write(SHORT_USE_FILE,*) (Rng(icell)%prop_annual_decid(i),i=1,FACETS) ! Proportion of facet that is annual plants (H_FACET) or deciduous (S_FACET and T_FACET)
       write(SHORT_USE_FILE,*) Rng(icell)%pot_evap                ! Potential evapotranspiration for the cell (cm/month)
       write(SHORT_USE_FILE,*) Rng(icell)%evaporation             ! Water evaporated from the soil and vegetation (cm/month)
       write(SHORT_USE_FILE,*) Rng(icell)%snow                    ! Snowpack, in cm
       write(SHORT_USE_FILE,*) Rng(icell)%snow_liquid             ! Snowpack liquid water.
       write(SHORT_USE_FILE,*) Rng(icell)%melt                    ! Snow that melts from snowpack (cm water)
       write(SHORT_USE_FILE,*) Rng(icell)%pet_remaining           ! Potential evaporation decremented as steps are calculated.  Appears to be a bookkeeping tool.
       write(SHORT_USE_FILE,*) Rng(icell)%ppt_soil                ! Precipitation adjusted for snow accumulation and melt, and available to infiltrate the soil (cm)
       write(SHORT_USE_FILE,*) Rng(icell)%runoff                  ! Runoff from the rangeland cell
       write(SHORT_USE_FILE,*) Rng(icell)%ratio_water_pet         ! Ratio of available water to potential evapotranspiration
       write(SHORT_USE_FILE,*) Rng(icell)%pet_top_soil            ! Potential evaporation from top soil (cm/day)
       write(SHORT_USE_FILE,*) (Rng(icell)%n_leached(i),i=1,SOIL_LAYERS)  ! Nitrogen leached from soil (AMTLEA in Century)
       write(SHORT_USE_FILE,*) (Rng(icell)%asmos(i),i=1,SOIL_LAYERS)      ! Used in summing water
       write(SHORT_USE_FILE,*) (Rng(icell)%amov(i),i=1,SOIL_LAYERS)       ! Used in summing water movement
       write(SHORT_USE_FILE,*) Rng(icell)%storm_flow              ! Storm flow
       write(SHORT_USE_FILE,*) Rng(icell)%holding_tank            ! Stores water temporarily.  Was asmos(layers+1) in H2OLos
       write(SHORT_USE_FILE,*) Rng(icell)%transpiration           ! Transpiration water loss
       write(SHORT_USE_FILE,*) (Rng(icell)%relative_water_content(i),i=1,SOIL_LAYERS) ! Used to initialize and during simulation in CENTURY. Here, only during simulation
       write(SHORT_USE_FILE,*) (Rng(icell)%water_available(i),i=1,3)                  ! Water available to plants, available for growth (1), survival (2), and in the two top layers (3)
       write(SHORT_USE_FILE,*) Rng(icell)%annual_evapotranspiration           ! Annual actual evapotranspiration
       write(SHORT_USE_FILE,*) Rng(icell)%total_aground_live_biomass      ! Total aboveground green biomass (g/m^2)
       write(SHORT_USE_FILE,*) Rng(icell)%total_bground_live_biomass      ! Total belowground green biomass (g/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%total_litter_carbon(i),i=1,2)              ! Average monthly litter carbon (g/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%total_litter_nitrogen(i),i=1,2)            ! Average monthly litter carbon (g/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%root_shoot_ratio(i),i=1,FACETS)         ! Root shoot ratio
       write(SHORT_USE_FILE,*) Rng(icell)%tree_basal_area                  ! Basal area for trees
       write(SHORT_USE_FILE,*) Rng(icell)%soil_surface_temperature         ! Average soil surface temperature (C)
       write(SHORT_USE_FILE,*) (Rng(icell)%sand(i),i=1,4)                 ! The percent sand in the soil
       write(SHORT_USE_FILE,*) (Rng(icell)%silt(i),i=1,4)                 ! The percent silt in the soil
       write(SHORT_USE_FILE,*) (Rng(icell)%clay(i),i=1,4)                 ! The percent clay in the soil
       write(SHORT_USE_FILE,*) (Rng(icell)%mineral_nitrogen(i),i=1,4)     ! Mineral nitrogen content for layer  (g/m2)
       write(SHORT_USE_FILE,*) (Rng(icell)%field_capacity(i),i=1,4)       ! Field capacity for four soils layers shown above.
       write(SHORT_USE_FILE,*) (Rng(icell)%wilting_point(i),i=1,4)        ! Wilting point for four soil layers shown above.
       write(SHORT_USE_FILE,*) Rng(icell)%soil_total_carbon       ! grams per square meter
       write(SHORT_USE_FILE,*) (Rng(icell)%tree_carbon(i),i=1,WOODY_PARTS)      ! Tree carbon in its components.   These must all be merged or otherwise crosswalked at some point.
       write(SHORT_USE_FILE,*) (Rng(icell)%tree_nitrogen(i),i=1,WOODY_PARTS)    ! Tree nitrogen in its components.   These must all be merged or otherwise crosswalked at some point.
       write(SHORT_USE_FILE,*) (Rng(icell)%shrub_carbon(i),i=1,WOODY_PARTS)      ! Shrub carbon in its components.   These must all be merged or otherwise crosswalked at some point.
       write(SHORT_USE_FILE,*) (Rng(icell)%shrub_nitrogen(i),i=1,WOODY_PARTS)    ! Shrub nitrogen in its components.   These must all be merged or otherwise crosswalked at some point.
       write(SHORT_USE_FILE,*) (Rng(icell)%carbon_nitrogen_ratio(i),i=1,2)        ! Carbon to nitrogen ratio, SURFACE, SOIL
       write(SHORT_USE_FILE,*) (Rng(icell)%fast_soil_carbon(i),i=1,2)             ! Soil organic matter carbon, surface and soil  g/m2  (SOM1C in Century)
       write(SHORT_USE_FILE,*) Rng(icell)%intermediate_soil_carbon        ! Intermediate soil carbon   g/m2  (SOMC2 in Century)
       write(SHORT_USE_FILE,*) Rng(icell)%passive_soil_carbon             ! Passive soil carbon    g/m2   (SOMC3 in Century)
       write(SHORT_USE_FILE,*) (Rng(icell)%fast_soil_nitrogen(i),i=1,2)           ! Soil organic matter nitrogen, surface and soil  g/m2  (SOM1E in Century and SSOM1E in Savanna)
       write(SHORT_USE_FILE,*) Rng(icell)%intermediate_soil_nitrogen      ! Intermediate soil nitrogen   g/m2  (SOM2E in Century)
       write(SHORT_USE_FILE,*) Rng(icell)%passive_soil_nitrogen           ! Passive soil nitrogen   g/m2   (SOM3E in Century)
       write(SHORT_USE_FILE,*) Rng(icell)%potential_production                    ! Calculated potential production for the cell, an index.  Based on soil temperature, so not specific to facets.
       write(SHORT_USE_FILE,*) (Rng(icell)%belowground_pot_production(i),i=1,V_LYRS)      ! BIOMASS, Belowground potential production in g/m2
       write(SHORT_USE_FILE,*) (Rng(icell)%aboveground_pot_production(i),i=1,V_LYRS)      ! BIOMASS, Abovegroudn potential production in g/m2
       write(SHORT_USE_FILE,*) (Rng(icell)%total_pot_production(i),i=1,V_LYRS)            ! BIOMASS, Calculate total potential production, in g/m2 with all the corrections in place.
       write(SHORT_USE_FILE,*) (Rng(icell)%co2_effect_on_production(i),i=1,FACETS)        ! Calculated effect of CO2 increasing from 350 to 700 ppm on grassland production, per facet
       write(SHORT_USE_FILE,*) (Rng(icell)%total_pot_prod_limited_by_n(i),i=1,V_LYRS)     ! Total potential production reflecting limits due to nitrogen in place in g/m2 (EPRODL)
       write(SHORT_USE_FILE,*) Rng(icell)%monthly_net_primary_production                 ! Monthly net primary production in g/m2
       write(SHORT_USE_FILE,*) Rng(icell)%fraction_live_removed_grazing   ! Fraction of live forage removed by grazing  (FLGREM in CENTURY)
       write(SHORT_USE_FILE,*) Rng(icell)%fraction_dead_removed_grazing   ! Fraction of dead forage removed by grazing  (FDGREM in CENTURY)
       write(SHORT_USE_FILE,*) Rng(icell)%temp_effect_on_decomp           ! Temperature effect on decomposition (TFUNC in CENTURY Cycle.f)  (index)
       write(SHORT_USE_FILE,*) Rng(icell)%water_effect_on_decomp          ! Water effect on decomposition (index)  (Aboveground and belowground entries in CENTURY set to equal, so distinction not made here)
       write(SHORT_USE_FILE,*) Rng(icell)%anerobic_effect_on_decomp       ! Anerobic effects on decomposition  (index)  (EFFANT in Savanna)
       write(SHORT_USE_FILE,*) Rng(icell)%all_effects_on_decomp           ! Combined effects on decomposition, which in Savanna includes anerobic (CYCLE.F)  (index)
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_fine_root_carbon(i),i=1,FACETS)       !,4)    ! Dead fine root carbon of the four types cited above.
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_fine_root_nitrogen(i),i=1,FACETS)     !,4)    ! Dead fine root nitrogen of the four types cited above.
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_standing_carbon(i),i=1,FACETS)        !,4)    ! Standing dead carbon of leaf and stem, of the four types cited above.
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_standing_nitrogen(i),i=1,FACETS)      !,4)    ! Standing dead nitrogen of leaf and stem, of the four types cited above.
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_seed_carbon(i),i=1,FACETS)            !,4)    ! Dead seed carbon of the four types cited above.   (gC/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_seed_nitrogen(i),i=1,FACETS)          !,4)    ! Dead seed nitrogen of the four types cited above.           (units?)
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_leaf_carbon(i),i=1,FACETS)            !,4)    ! Dead leaf carbon of the four types cited above.   (gC/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_leaf_nitrogen(i),i=1,FACETS)          !,4)    ! Dead leaf nitrogen of the four types cited above.
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_fine_branch_carbon(i),i=1,FACETS)     !,4)    ! Dead fine branch carbon of the four types cited above.   (gC/m^2)
       write(SHORT_USE_FILE,*) Rng(icell)%dead_total_fine_branch_carbon               ! Dead fine branch carbon, summed across facets
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_fine_branch_nitrogen(i),i=1,FACETS)   !,4)    ! Dead fine branch nitrogen of the four types cited above.
       write(SHORT_USE_FILE,*) Rng(icell)%dead_total_fine_branch_nitrogen             ! Dead fine branch nitrogen, summed across facets
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_coarse_root_carbon(i),i=1,FACETS)     !,4)    ! Dead coarse root carbon of the four types cited above.   (gC/m^2)
       write(SHORT_USE_FILE,*) Rng(icell)%dead_total_coarse_root_carbon               ! Dead total coarse root carbon, summed across facets
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_coarse_root_nitrogen(i),i=1,FACETS)   !,4)    ! Dead coarse root nitrogen of the four types cited above.
       write(SHORT_USE_FILE,*) Rng(icell)%dead_total_coarse_root_nitrogen             ! Dead total coarse root nitrogen, summed across facets
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_coarse_branch_carbon(i),i=1,FACETS)     !,4)    ! Dead coarse wood carbon of the four types cited above.   (gC/m^2)
       write(SHORT_USE_FILE,*) Rng(icell)%dead_total_coarse_branch_carbon               ! Dead total coarse wood carbon, summed across facets
       write(SHORT_USE_FILE,*) (Rng(icell)%dead_coarse_branch_nitrogen(i),i=1,FACETS)   !,4)    ! Dead coarse wood nitrogen of the four types cited above.
       write(SHORT_USE_FILE,*) Rng(icell)%dead_total_coarse_branch_nitrogen              ! Dead total coarse wood nitrogen, summed across facets
       write(SHORT_USE_FILE,*) (Rng(icell)%lignin_fine_root(i),i=1,FACETS)               ! Fine root lignin concentration
       write(SHORT_USE_FILE,*) (Rng(icell)%lignin_coarse_root(i),i=1,FACETS)             ! Coarse root lignin concentration
       write(SHORT_USE_FILE,*) (Rng(icell)%lignin_fine_branch(i),i=1,FACETS)             ! Fine branch lignin concentration
       write(SHORT_USE_FILE,*) (Rng(icell)%lignin_coarse_branch(i),i=1,FACETS)           ! Coarse branch lignin concentration
       write(SHORT_USE_FILE,*) (Rng(icell)%lignin_leaf(i),i=1,FACETS)                    ! Leaf lignin concentration
       write(SHORT_USE_FILE,*) ((Rng(icell)%plant_lignin_fraction(i,j),i=1,FACETS),j=1,2)        ! Lignin in structural residue, at the surface (1) and in the soil (2)  (STRLIG)
       write(SHORT_USE_FILE,*) (Rng(icell)%litter_structural_carbon(i),i=1,2)            ! Litter structural carbon at the surface (1) and in the soil (2)  (STRCIS, or in Savanna, SSTRCIS, with unlabeled and labeled merged)
       write(SHORT_USE_FILE,*) (Rng(icell)%litter_metabolic_carbon(i),i=1,2)             ! Litter metabolic carbon at the surface (1) and in the soil (2)  (METCIS, or in Savanna, SMETCIS)
       write(SHORT_USE_FILE,*) (Rng(icell)%litter_structural_nitrogen(i),i=1,2)          ! Litter structural nitrogen at the surface (1) and in the soil (2)  (STRUCE, or in Savanna, SSTRUCE, with STRUCE named for "elements"  I am only including nitrogen, as in Savanna, so dropping the name)
       write(SHORT_USE_FILE,*) (Rng(icell)%litter_metabolic_nitrogen(i),i=1,2)           ! Litter structural nitrogen at the surface (1) and in the soil (2)  (METABE, or in Savanna, SSTRUCE, with STRUCE named for "elements"  I am only including nitrogen, as in Savanna, so dropping the name)
       write(SHORT_USE_FILE,*) (Rng(icell)%maintain_respiration(i),i=1,FACETS)           ! Maintainence respiration
       write(SHORT_USE_FILE,*) (Rng(icell)%phenology(i),i=1,FACETS)                     ! Phenological stage, a continuous variable from 0 to 4.
       write(SHORT_USE_FILE,*) (Rng(icell)%fine_root_carbon(i),i=1,FACETS)              ! Fine root carbon    (gC/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%fine_root_nitrogen(i),i=1,FACETS)            ! Fine root nitrogen
       write(SHORT_USE_FILE,*) (Rng(icell)%seed_carbon(i),i=1,FACETS)                   ! Seed carbon         (gC/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%seed_nitrogen(i),i=1,FACETS)                 ! Seed nitrogen       (units?)
       write(SHORT_USE_FILE,*) (Rng(icell)%leaf_carbon(i),i=1,FACETS)                   ! Leaf carbon         (gC/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%leaf_nitrogen(i),i=1,FACETS)                 ! Leaf nitrogen
       write(SHORT_USE_FILE,*) (Rng(icell)%fine_branch_carbon(i),i=1,FACETS)            ! Fine branch carbon  (gC/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%fine_branch_nitrogen(i),i=1,FACETS)          ! Fine branch nitrogen
       write(SHORT_USE_FILE,*) (Rng(icell)%coarse_root_carbon(i),i=1,FACETS)            ! Coarse root carbon  (gC/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%coarse_root_nitrogen(i),i=1,FACETS)          ! Coarse root nitrogen
       write(SHORT_USE_FILE,*) (Rng(icell)%coarse_branch_carbon(i),i=1,FACETS)          ! Coarse branch carbon  (gC/m^2)
       write(SHORT_USE_FILE,*) (Rng(icell)%coarse_branch_nitrogen(i),i=1,FACETS)        ! Coarse branch nitrogen
       write(SHORT_USE_FILE,*) (Rng(icell)%stored_nitrogen(i),i=1,FACETS)               ! Stored nitrogen (CRPSTG in Century GROWTH, STORAGE in RESTRP).  I can't find where this is initialized, except for a gridded system input.  Assumed 0 for now, but here as a placeholder.
       write(SHORT_USE_FILE,*) (Rng(icell)%plant_nitrogen_fixed(i),i=1,FACETS)          ! Plant nitrogen fixed
       write(SHORT_USE_FILE,*) (Rng(icell)%nitrogen_fixed(i),i=1,FACETS)                ! Nitrogen fixed.  Not sure what components distinguish it, just yet.  (NFIX)
       write(SHORT_USE_FILE,*) (Rng(icell)%respiration_flows(i),i=1,FACETS)             ! Maintenance respiration flows to storage pool  (MRSPSTG)
       write(SHORT_USE_FILE,*) (Rng(icell)%respiration_annual(i),i=1,FACETS)            ! Maintenance respiration flows for year         (MRSPANN)
       write(SHORT_USE_FILE,*) Rng(icell)%carbon_source_sink                    ! Carbon pool.  (g/m2)   (CSRSNK)    I don't know the utility of this, but incorporating it.
       write(SHORT_USE_FILE,*) Rng(icell)%nitrogen_source_sink                  ! Nitrogen pool.  (g/m2)   (ESRSNK)
       write(SHORT_USE_FILE,*) ((Rng(icell)%carbon_allocation(i,j),i=1,FACETS), j=1,WOODY_PARTS) ! Shrub carbon allocation, by proportion  (TREE_CFAC in Century, except statis here)  Brought into this array, even though it requires more memory, in case I do incorporate dynamic allocation at some point.
       write(SHORT_USE_FILE,*) (Rng(icell)%optimum_leaf_area_index(i),i=1,FACETS)       ! Optimum leaf area index
       write(SHORT_USE_FILE,*) (Rng(icell)%leaf_area_index(i),i=1,FACETS)               ! Leaf area index
       write(SHORT_USE_FILE,*) Rng(icell)%water_function                        ! Water function influencing mortality (AGWFUNC and BGWFUNC in Century, merged here since CYCLE assigns them equal and they start with the same value)
       write(SHORT_USE_FILE,*) Rng(icell)%fire_severity                         ! A score from 0 to 1 reflecting fire intensity
       write(SHORT_USE_FILE,*) Rng(icell)%burned_carbon                         ! The sum of carbon burned, only on the 1 m plots, not whole plant death
       write(SHORT_USE_FILE,*) Rng(icell)%burned_nitrogen                       ! The sum of nitrogen burned, only on the 1 m plots, not whole plant death
       write(SHORT_USE_FILE,*) Rng(icell)%fertilized_nitrogen_added             ! Total fertilized nitrogen added  (g/m2)
       write(SHORT_USE_FILE,*) Rng(icell)%fertilized_carbon_added               ! Total fertilized carbon added  (g/m2)
     end do                                                                                                            
  else
    write(ECHO_FILE,*) 'Unable to open the state variable output file for writing: ', Sim_Parm%state_var_file_out
  end if
  close(SHORT_USE_FILE)

end subroutine



subroutine Read_State
!**** Reads state variable results from a file name that is provided.
!****
!**** R. Boone   Last modified: December 28, 2013.
   use Parameter_Vars
   use Structures
   implicit none

   character(255)  str
   integer  i, j, icell, xtemp, ytemp, celltemp, lxtemp, lytemp, uxtemp, uytemp, rangetemp
   integer  zonetemp
   real     cellsizetemp

   open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//Sim_Parm%state_var_file_in, &
       STATUS='OLD', ACTION='READ', IOSTAT=ioerr)
     
   if (ioerr == 0) then
     ! See Write_State for additional comments.
     str = ''
     read(SHORT_USE_FILE,'(A255)') str                     ! x_dim    
     str = str(28:200)
     read(str,'(I10)') xtemp
     if (xtemp .ne. x_dim) then
       write(*,*) 'The x dimension in the simulation does not match that in'
       write(*,*) 'the state variable file. The value in the simulation is: ',x_dim
       write(*,*) 'whereas the value in the state variable file is: ', xtemp
       write(*,*) 'The simulation cannot continue.'
       stop
     end if
     str = ''
     read(SHORT_USE_FILE,'(A255)') str                     ! y_dim                                                      
     str = str(28:200)
     read(str,'(i10)') ytemp
     if (ytemp .ne. y_dim) then
       write(*,*) 'The y dimension in the simulation does not match that in'
       write(*,*) 'the state variable file. The value in the simulation is: ',y_dim
       write(*,*) 'whereas the value in the state variable file is: ', ytemp
       write(*,*) 'The simulation cannot continue.'
       stop
     end if
     str = ''
     read(SHORT_USE_FILE,'(A255)') str                     ! lower_x
     str = str(28:200)
     read(str,'(i10)') lxtemp
     if (lxtemp .ne. lower_x) then
       write(*,*) 'The lower x dimension in the simulation does not match that in'
       write(*,*) 'the state variable file. The value in the simulation is: ',lower_x
       write(*,*) 'whereas the value in the state variable file is: ', lxtemp
       write(*,*) 'The simulation cannot continue.'
       stop
     end if
     str = ''
     read(SHORT_USE_FILE,'(A255)') str                     ! lower_y
     str = str(28:200)
     read(str,'(i10)') lytemp
     if (lytemp .ne. lower_y) then
       write(*,*) 'The lower y dimension in the simulation does not match that in'
       write(*,*) 'the state variable file. The value in the simulation is: ',lower_y
       write(*,*) 'whereas the value in the state variable file is: ', lytemp
       write(*,*) 'The simulation cannot continue.'
       stop
     end if
     str = ''
     read(SHORT_USE_FILE,'(A255)') str                     ! upper_x
     str = str(28:200)
     read(str,'(i10)') uxtemp
     if (uxtemp .ne. upper_x) then
       write(*,*) 'The upper x dimension in the simulation does not match that in'
       write(*,*) 'the state variable file. The value in the simulation is: ',upper_x
       write(*,*) 'whereas the value in the state variable file is: ', uxtemp
       write(*,*) 'The simulation cannot continue.'
       stop
     end if
     str = ''
     read(SHORT_USE_FILE,'(A255)') str                     ! upper_y
     str = str(28:200)
     read(str,'(i10)') uytemp
     if (uytemp .ne. upper_y) then
       write(*,*) 'The lower y dimension in the simulation does not match that in'
       write(*,*) 'the state variable file. The value in the simulation is: ',lower_y
       write(*,*) 'whereas the value in the state variable file is: ', lytemp
       write(*,*) 'The simulation cannot continue.'
       stop
     end if
     str = ''
     read(SHORT_USE_FILE,'(A255)') str                     ! cellsize
     str = str(28:200)
     read(str,'(F10.8)') cellsizetemp
     if ( int(cellsizetemp * 10000) .ne. int(cellsize * 10000) ) then             ! Avoiding equality test of floating point 
       write(*,*) 'The cellsize in the simulation does not match that in'
       write(*,*) 'the state variable file. The value in the simulation is: ',cellsize
       write(*,*) 'whereas the value in the state variable file is: ', celltemp
       write(*,*) 'The simulation cannot continue.'
       stop
     end if
     str = ''
     read(SHORT_USE_FILE,'(A255)') str                     ! range_cells
     str = str(28:200)
     read(str,'(i10)') rangetemp
     if (rangetemp .ne. range_cells) then
       write(*,*) 'The number of rangeland cells in the simulation does not match that in'
       write(*,*) 'the state variable file. The value in the simulation is: ',range_cells
       write(*,*) 'whereas the value in the state variable file is: ', rangetemp
       write(*,*) 'The simulation cannot continue.'
       stop
     end if
     read(SHORT_USE_FILE,'(A255)') str                     ! bin_path
     read(SHORT_USE_FILE,'(A255)') str                     ! app_path
     read(SHORT_USE_FILE,'(A255)') str                     ! parm_path
     read(SHORT_USE_FILE,'(A255)') str                     ! out_path
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%zone_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%land_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%elev_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%class_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%class_legend
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%latitude_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%landscape_type_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%top_sand_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%top_silt_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%top_clay_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%top_gravel_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%top_bulk_density_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%top_organic_carbon_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%sub_sand_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%sub_silt_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%sub_clay_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%sub_gravel_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%sub_bulk_density_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%sub_organic_carbon_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%deciduous_tree_cover_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%evergreen_tree_cover_map     
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%shrub_cover_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%herb_cover_map
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%precip_average_map                                
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%temperature_average_map                           
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%precip_path_prefix                                
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%max_temp_path_prefix                              
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%min_temp_path_prefix                              
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%precip_temp_suffix                                
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%parms_file_name                                   
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%start_yr          Just for completeness, so a user knows the years used in a spin-up file                                    
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%end_yr                                            
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%state_var_flag                                    
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%state_var_file_out                                
     read(SHORT_USE_FILE,'(A255)') str                     ! Sim_Parm%state_var_file_in                                 
                                                                                                                       
     do icell = 1, range_cells                                                                                         
       ! Agreement between the STATE VARIABLE file and the simulation will be kept in-sync by checking the first four 
       ! values.  That will ensure that the file and simulation are in sync.  Simulations comparing 10 year simulations
       ! to two five year simulations can identify if the state variable files are yielding equivalent results.
       str = ''
       read(SHORT_USE_FILE,'(A255)') str                     ! Zone ID
       str = str(15:200)
       read(str,'(i10)') zonetemp
       if (zonetemp .ne. Globe(Rng(icell)%x, Rng(icell)%y)%zone) then
         write(*,*) 'The zonal id does not match what was read in from the map.'
         write(*,*) 'The value in the simulation is: ',Globe(Rng(icell)%x, Rng(icell)%y)%zone
         write(*,*) 'whereas the value in the state variable file is: ', zonetemp
         write(*,*) 'The simulation cannot continue.'
         stop
       end if
       str = ''
       read(SHORT_USE_FILE,'(A255)') str                     ! X
       str = str(15:200)
       read(str,'(i10)') xtemp
       if (xtemp .ne. Rng(icell)%x) then
         write(*,*) 'The x dimension does not match what was read in from the map.'
         write(*,*) 'The value in the simulation is: ',Rng(icell)%x
         write(*,*) 'whereas the value in the state variable file is: ', xtemp
         write(*,*) 'The simulation cannot continue.'
         stop
       end if
       str = ''
       read(SHORT_USE_FILE,'(A255)') str                     ! Y
       str = str(15:200)
       read(str,'(i10)') ytemp
       if (ytemp .ne. Rng(icell)%y) then
         write(*,*) 'The y dimension does not match what was read in from the map.'
         write(*,*) 'The value in the simulation is: ',Rng(icell)%y
         write(*,*) 'whereas the value in the state variable file is: ', zonetemp
         write(*,*) 'The simulation cannot continue.'
         stop
       end if
       str = ''
       read(SHORT_USE_FILE,'(A255)') str                     ! range_type
       str = str(15:200)
       read(str,'(i10)') rangetemp
       if (rangetemp .ne. Rng(icell)%range_type) then
         write(*,*) 'The range type does not match what was read in from the map.'
         write(*,*) 'The value in the simulation is: ',Rng(icell)%range_type
         write(*,*) 'whereas the value in the state variable file is: ', rangetemp
         write(*,*) 'The simulation cannot continue.'
         stop
       end if
       read(SHORT_USE_FILE,*) Rng(icell)%last_month_day_length     ! The day length of the previous month, to know when spring and fall come.
       read(SHORT_USE_FILE,*) Rng(icell)%day_length_increasing     ! Increasing or decreasing day length, comparing the current to previous day lengths.
       read(SHORT_USE_FILE,*) Rng(icell)%day_length                ! Day length, calculated based on latitude and month
       read(SHORT_USE_FILE,*) Rng(icell)%heat_accumulation         ! Heat accumulation above a base temperature (e.g., 4.4 C in Boone (1999))
       read(SHORT_USE_FILE,*) (Rng(icell)%facet_cover(i),i=1,FACETS)       ! The proportion occupied by each facet
       read(SHORT_USE_FILE,*) (Rng(icell)%total_population(i),i=1,V_LYRS)  ! The total population of each vegetation layer
       read(SHORT_USE_FILE,*) Rng(icell)%bare_cover                ! Bare cover stored, rather than continually summing the three facets.
       read(SHORT_USE_FILE,*) (Rng(icell)%prop_annual_decid(i),i=1,FACETS) ! Proportion of facet that is annual plants (H_FACET) or deciduous (S_FACET and T_FACET)
       read(SHORT_USE_FILE,*) Rng(icell)%pot_evap                ! Potential evapotranspiration for the cell (cm/month)
       read(SHORT_USE_FILE,*) Rng(icell)%evaporation             ! Water evaporated from the soil and vegetation (cm/month)
       read(SHORT_USE_FILE,*) Rng(icell)%snow                    ! Snowpack, in cm
       read(SHORT_USE_FILE,*) Rng(icell)%snow_liquid             ! Snowpack liquid water.
       read(SHORT_USE_FILE,*) Rng(icell)%melt                    ! Snow that melts from snowpack (cm water)
       read(SHORT_USE_FILE,*) Rng(icell)%pet_remaining           ! Potential evaporation decremented as steps are calculated.  Appears to be a bookkeeping tool.
       read(SHORT_USE_FILE,*) Rng(icell)%ppt_soil                ! Precipitation adjusted for snow accumulation and melt, and available to infiltrate the soil (cm)
       read(SHORT_USE_FILE,*) Rng(icell)%runoff                  ! Runoff from the rangeland cell
       read(SHORT_USE_FILE,*) Rng(icell)%ratio_water_pet         ! Ratio of available water to potential evapotranspiration
       read(SHORT_USE_FILE,*) Rng(icell)%pet_top_soil            ! Potential evaporation from top soil (cm/day)
       read(SHORT_USE_FILE,*) (Rng(icell)%n_leached(i),i=1,SOIL_LAYERS)  ! Nitrogen leached from soil (AMTLEA in Century)
       read(SHORT_USE_FILE,*) (Rng(icell)%asmos(i),i=1,SOIL_LAYERS)      ! Used in summing water
       read(SHORT_USE_FILE,*) (Rng(icell)%amov(i),i=1,SOIL_LAYERS)       ! Used in summing water movement
       read(SHORT_USE_FILE,*) Rng(icell)%storm_flow              ! Storm flow
       read(SHORT_USE_FILE,*) Rng(icell)%holding_tank            ! Stores water temporarily.  Was asmos(layers+1) in H2OLos
       read(SHORT_USE_FILE,*) Rng(icell)%transpiration           ! Transpiration water loss
       read(SHORT_USE_FILE,*) (Rng(icell)%relative_water_content(i),i=1,SOIL_LAYERS) ! Used to initialize and during simulation in CENTURY. Here, only during simulation
       read(SHORT_USE_FILE,*) (Rng(icell)%water_available(i),i=1,3)                  ! Water available to plants, available for growth (1), survival (2), and in the two top layers (3)
       read(SHORT_USE_FILE,*) Rng(icell)%annual_evapotranspiration           ! Annual actual evapotranspiration
       read(SHORT_USE_FILE,*) Rng(icell)%total_aground_live_biomass          ! Total aboveground live biomass (g/m^2)
       read(SHORT_USE_FILE,*) Rng(icell)%total_bground_live_biomass          ! Total belowground live biomass (g/m^2)       
       read(SHORT_USE_FILE,*) (Rng(icell)%total_litter_carbon(i),i=1,2)      ! Average monthly litter carbon (g/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%total_litter_nitrogen(i),i=1,2)    ! Average monthly litter carbon (g/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%root_shoot_ratio(i),i=1,FACETS)    ! Root shoot ratio
       read(SHORT_USE_FILE,*) Rng(icell)%tree_basal_area                  ! Basal area for trees
       read(SHORT_USE_FILE,*) Rng(icell)%soil_surface_temperature         ! Average soil surface temperature (C)
       read(SHORT_USE_FILE,*) (Rng(icell)%sand(i),i=1,4)                 ! The percent sand in the soil
       read(SHORT_USE_FILE,*) (Rng(icell)%silt(i),i=1,4)                 ! The percent silt in the soil
       read(SHORT_USE_FILE,*) (Rng(icell)%clay(i),i=1,4)                 ! The percent clay in the soil
       read(SHORT_USE_FILE,*) (Rng(icell)%mineral_nitrogen(i),i=1,4)     ! Mineral nitrogen content for layer  (g/m2)
       read(SHORT_USE_FILE,*) (Rng(icell)%field_capacity(i),i=1,4)       ! Field capacity for four soils layers shown above.
       read(SHORT_USE_FILE,*) (Rng(icell)%wilting_point(i),i=1,4)        ! Wilting point for four soil layers shown above.
       read(SHORT_USE_FILE,*) Rng(icell)%soil_total_carbon       ! grams per square meter
       read(SHORT_USE_FILE,*) (Rng(icell)%tree_carbon(i),i=1,WOODY_PARTS)      ! Tree carbon in its components.   These must all be merged or otherwise crosswalked at some point.
       read(SHORT_USE_FILE,*) (Rng(icell)%tree_nitrogen(i),i=1,WOODY_PARTS)    ! Tree nitrogen in its components.   These must all be merged or otherwise crosswalked at some point.
       read(SHORT_USE_FILE,*) (Rng(icell)%shrub_carbon(i),i=1,WOODY_PARTS)      ! Shrub carbon in its components.   These must all be merged or otherwise crosswalked at some point.
       read(SHORT_USE_FILE,*) (Rng(icell)%shrub_nitrogen(i),i=1,WOODY_PARTS)    ! Shrub nitrogen in its components.   These must all be merged or otherwise crosswalked at some point.
       read(SHORT_USE_FILE,*) (Rng(icell)%carbon_nitrogen_ratio(i),i=1,2)        ! Carbon to nitrogen ratio, SURFACE, SOIL
       read(SHORT_USE_FILE,*) (Rng(icell)%fast_soil_carbon(i),i=1,2)             ! Soil organic matter carbon, surface and soil  g/m2  (SOM1C in Century)
       read(SHORT_USE_FILE,*) Rng(icell)%intermediate_soil_carbon        ! Intermediate soil carbon   g/m2  (SOMC2 in Century)
       read(SHORT_USE_FILE,*) Rng(icell)%passive_soil_carbon             ! Passive soil carbon    g/m2   (SOMC3 in Century)
       read(SHORT_USE_FILE,*) (Rng(icell)%fast_soil_nitrogen(i),i=1,2)           ! Soil organic matter nitrogen, surface and soil  g/m2  (SOM1E in Century and SSOM1E in Savanna)
       read(SHORT_USE_FILE,*) Rng(icell)%intermediate_soil_nitrogen      ! Intermediate soil nitrogen   g/m2  (SOM2E in Century)
       read(SHORT_USE_FILE,*) Rng(icell)%passive_soil_nitrogen           ! Passive soil nitrogen   g/m2   (SOM3E in Century)
       read(SHORT_USE_FILE,*) Rng(icell)%potential_production                    ! Calculated potential production for the cell, an index.  Based on soil temperature, so not specific to facets.
       read(SHORT_USE_FILE,*) (Rng(icell)%belowground_pot_production(i),i=1,V_LYRS)      ! BIOMASS, Belowground potential production in g/m2
       read(SHORT_USE_FILE,*) (Rng(icell)%aboveground_pot_production(i),i=1,V_LYRS)      ! BIOMASS, Abovegroudn potential production in g/m2
       read(SHORT_USE_FILE,*) (Rng(icell)%total_pot_production(i),i=1,V_LYRS)            ! BIOMASS, Calculate total potential production, in g/m2 with all the corrections in place.
       read(SHORT_USE_FILE,*) (Rng(icell)%co2_effect_on_production(i),i=1,FACETS)        ! Calculated effect of CO2 increasing from 350 to 700 ppm on grassland production, per facet
       read(SHORT_USE_FILE,*) (Rng(icell)%total_pot_prod_limited_by_n(i),i=1,V_LYRS)     ! Total potential production reflecting limits due to nitrogen in place in g/m2 (EPRODL)
       read(SHORT_USE_FILE,*) Rng(icell)%monthly_net_primary_production                 ! Monthly net primary production in g/m2
       read(SHORT_USE_FILE,*) Rng(icell)%fraction_live_removed_grazing   ! Fraction of live forage removed by grazing  (FLGREM in CENTURY)
       read(SHORT_USE_FILE,*) Rng(icell)%fraction_dead_removed_grazing   ! Fraction of dead forage removed by grazing  (FDGREM in CENTURY)
       read(SHORT_USE_FILE,*) Rng(icell)%temp_effect_on_decomp           ! Temperature effect on decomposition (TFUNC in CENTURY Cycle.f)  (index)
       read(SHORT_USE_FILE,*) Rng(icell)%water_effect_on_decomp          ! Water effect on decomposition (index)  (Aboveground and belowground entries in CENTURY set to equal, so distinction not made here)
       read(SHORT_USE_FILE,*) Rng(icell)%anerobic_effect_on_decomp       ! Anerobic effects on decomposition  (index)  (EFFANT in Savanna)
       read(SHORT_USE_FILE,*) Rng(icell)%all_effects_on_decomp           ! Combined effects on decomposition, which in Savanna includes anerobic (CYCLE.F)  (index)
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_fine_root_carbon(i),i=1,FACETS)       !,4)    ! Dead fine root carbon of the four types cited above.
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_fine_root_nitrogen(i),i=1,FACETS)     !,4)    ! Dead fine root nitrogen of the four types cited above.
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_standing_carbon(i),i=1,FACETS)        !,4)    ! Standing dead carbon of leaf and stem, of the four types cited above.
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_standing_nitrogen(i),i=1,FACETS)      !,4)    ! Standing dead nitrogen of leaf and stem, of the four types cited above.
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_seed_carbon(i),i=1,FACETS)            !,4)    ! Dead seed carbon of the four types cited above.   (gC/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_seed_nitrogen(i),i=1,FACETS)          !,4)    ! Dead seed nitrogen of the four types cited above.           (units?)
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_leaf_carbon(i),i=1,FACETS)            !,4)    ! Dead leaf carbon of the four types cited above.   (gC/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_leaf_nitrogen(i),i=1,FACETS)          !,4)    ! Dead leaf nitrogen of the four types cited above.
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_fine_branch_carbon(i),i=1,FACETS)     !,4)    ! Dead fine branch carbon of the four types cited above.   (gC/m^2)
       read(SHORT_USE_FILE,*) Rng(icell)%dead_total_fine_branch_carbon               ! Dead fine branch carbon, summed across facets
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_fine_branch_nitrogen(i),i=1,FACETS)   !,4)    ! Dead fine branch nitrogen of the four types cited above.
       read(SHORT_USE_FILE,*) Rng(icell)%dead_total_fine_branch_nitrogen             ! Dead fine branch nitrogen, summed across facets
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_coarse_root_carbon(i),i=1,FACETS)     !,4)    ! Dead coarse root carbon of the four types cited above.   (gC/m^2)
       read(SHORT_USE_FILE,*) Rng(icell)%dead_total_coarse_root_carbon               ! Dead total coarse root carbon, summed across facets
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_coarse_root_nitrogen(i),i=1,FACETS)   !,4)    ! Dead coarse root nitrogen of the four types cited above.
       read(SHORT_USE_FILE,*) Rng(icell)%dead_total_coarse_root_nitrogen             ! Dead total coarse root nitrogen, summed across facets
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_coarse_branch_carbon(i),i=1,FACETS)     !,4)    ! Dead coarse wood carbon of the four types cited above.   (gC/m^2)
       read(SHORT_USE_FILE,*) Rng(icell)%dead_total_coarse_branch_carbon               ! Dead total coarse wood carbon, summed across facets
       read(SHORT_USE_FILE,*) (Rng(icell)%dead_coarse_branch_nitrogen(i),i=1,FACETS)   !,4)    ! Dead coarse wood nitrogen of the four types cited above.
       read(SHORT_USE_FILE,*) Rng(icell)%dead_total_coarse_branch_nitrogen              ! Dead total coarse wood nitrogen, summed across facets
       read(SHORT_USE_FILE,*) (Rng(icell)%lignin_fine_root(i),i=1,FACETS)               ! Fine root lignin concentration
       read(SHORT_USE_FILE,*) (Rng(icell)%lignin_coarse_root(i),i=1,FACETS)             ! Coarse root lignin concentration
       read(SHORT_USE_FILE,*) (Rng(icell)%lignin_fine_branch(i),i=1,FACETS)             ! Fine branch lignin concentration
       read(SHORT_USE_FILE,*) (Rng(icell)%lignin_coarse_branch(i),i=1,FACETS)           ! Coarse branch lignin concentration
       read(SHORT_USE_FILE,*) (Rng(icell)%lignin_leaf(i),i=1,FACETS)                    ! Leaf lignin concentration
       read(SHORT_USE_FILE,*) ((Rng(icell)%plant_lignin_fraction(i,j),i=1,FACETS),j=1,2)        ! Lignin in structural residue, at the surface (1) and in the soil (2)  (STRLIG)
       read(SHORT_USE_FILE,*) (Rng(icell)%litter_structural_carbon(i),i=1,2)            ! Litter structural carbon at the surface (1) and in the soil (2)  (STRCIS, or in Savanna, SSTRCIS, with unlabeled and labeled merged)
       read(SHORT_USE_FILE,*) (Rng(icell)%litter_metabolic_carbon(i),i=1,2)             ! Litter metabolic carbon at the surface (1) and in the soil (2)  (METCIS, or in Savanna, SMETCIS)
       read(SHORT_USE_FILE,*) (Rng(icell)%litter_structural_nitrogen(i),i=1,2)          ! Litter structural nitrogen at the surface (1) and in the soil (2)  (STRUCE, or in Savanna, SSTRUCE, with STRUCE named for "elements"  I am only including nitrogen, as in Savanna, so dropping the name)
       read(SHORT_USE_FILE,*) (Rng(icell)%litter_metabolic_nitrogen(i),i=1,2)           ! Litter structural nitrogen at the surface (1) and in the soil (2)  (METABE, or in Savanna, SSTRUCE, with STRUCE named for "elements"  I am only including nitrogen, as in Savanna, so dropping the name)
       read(SHORT_USE_FILE,*) (Rng(icell)%maintain_respiration(i),i=1,FACETS)           ! Maintainence respiration
       read(SHORT_USE_FILE,*) (Rng(icell)%phenology(i),i=1,FACETS)                     ! Phenological stage, a continuous variable from 0 to 4.
       read(SHORT_USE_FILE,*) (Rng(icell)%fine_root_carbon(i),i=1,FACETS)              ! Fine root carbon    (gC/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%fine_root_nitrogen(i),i=1,FACETS)            ! Fine root nitrogen
       read(SHORT_USE_FILE,*) (Rng(icell)%seed_carbon(i),i=1,FACETS)                   ! Seed carbon         (gC/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%seed_nitrogen(i),i=1,FACETS)                 ! Seed nitrogen       (units?)
       read(SHORT_USE_FILE,*) (Rng(icell)%leaf_carbon(i),i=1,FACETS)                   ! Leaf carbon         (gC/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%leaf_nitrogen(i),i=1,FACETS)                 ! Leaf nitrogen
       read(SHORT_USE_FILE,*) (Rng(icell)%fine_branch_carbon(i),i=1,FACETS)            ! Fine branch carbon  (gC/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%fine_branch_nitrogen(i),i=1,FACETS)          ! Fine branch nitrogen
       read(SHORT_USE_FILE,*) (Rng(icell)%coarse_root_carbon(i),i=1,FACETS)            ! Coarse root carbon  (gC/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%coarse_root_nitrogen(i),i=1,FACETS)          ! Coarse root nitrogen
       read(SHORT_USE_FILE,*) (Rng(icell)%coarse_branch_carbon(i),i=1,FACETS)          ! Coarse branch carbon  (gC/m^2)
       read(SHORT_USE_FILE,*) (Rng(icell)%coarse_branch_nitrogen(i),i=1,FACETS)        ! Coarse branch nitrogen
       read(SHORT_USE_FILE,*) (Rng(icell)%stored_nitrogen(i),i=1,FACETS)               ! Stored nitrogen (CRPSTG in Century GROWTH, STORAGE in RESTRP).  I can't find where this is initialized, except for a gridded system input.  Assumed 0 for now, but here as a placeholder.
       read(SHORT_USE_FILE,*) (Rng(icell)%plant_nitrogen_fixed(i),i=1,FACETS)          ! Plant nitrogen fixed
       read(SHORT_USE_FILE,*) (Rng(icell)%nitrogen_fixed(i),i=1,FACETS)                ! Nitrogen fixed.  Not sure what components distinguish it, just yet.  (NFIX)
       read(SHORT_USE_FILE,*) (Rng(icell)%respiration_flows(i),i=1,FACETS)             ! Maintenance respiration flows to storage pool  (MRSPSTG)
       read(SHORT_USE_FILE,*) (Rng(icell)%respiration_annual(i),i=1,FACETS)            ! Maintenance respiration flows for year         (MRSPANN)
       read(SHORT_USE_FILE,*) Rng(icell)%carbon_source_sink                    ! Carbon pool.  (g/m2)   (CSRSNK)    I don't know the utility of this, but incorporating it.
       read(SHORT_USE_FILE,*) Rng(icell)%nitrogen_source_sink                  ! Nitrogen pool.  (g/m2)   (ESRSNK)
       read(SHORT_USE_FILE,*) ((Rng(icell)%carbon_allocation(i,j),i=1,FACETS), j=1,WOODY_PARTS) ! Shrub carbon allocation, by proportion  (TREE_CFAC in Century, except statis here)  Brought into this array, even though it requires more memory, in case I do incorporate dynamic allocation at some point.
       read(SHORT_USE_FILE,*) (Rng(icell)%optimum_leaf_area_index(i),i=1,FACETS)       ! Optimum leaf area index
       read(SHORT_USE_FILE,*) (Rng(icell)%leaf_area_index(i),i=1,FACETS)               ! Leaf area index
       read(SHORT_USE_FILE,*) Rng(icell)%water_function                        ! Water function influencing mortality (AGWFUNC and BGWFUNC in Century, merged here since CYCLE assigns them equal and they start with the same value)
       read(SHORT_USE_FILE,*) Rng(icell)%fire_severity                         ! A score from 0 to 1 reflecting fire intensity
       read(SHORT_USE_FILE,*) Rng(icell)%burned_carbon                         ! The sum of carbon burned, only on the 1 m plots, not whole plant death
       read(SHORT_USE_FILE,*) Rng(icell)%burned_nitrogen                       ! The sum of nitrogen burned, only on the 1 m plots, not whole plant death
       read(SHORT_USE_FILE,*) Rng(icell)%fertilized_nitrogen_added             ! Total fertilized nitrogen added  (g/m2)
       read(SHORT_USE_FILE,*) Rng(icell)%fertilized_carbon_added               ! Total fertilized carbon added  (g/m2)
     end do                                                                                                            
  else
    write(ECHO_FILE,*) 'Unable to open the state variable output file for reading: ', Sim_Parm%state_var_file_out
  end if
  close(SHORT_USE_FILE)

end subroutine