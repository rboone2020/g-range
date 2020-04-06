module Structures
  ! The structures used in G_RANGE are defined here.  The primary structure defines the cells comprising the landscape.
  ! This will store the information for each cell, including oceans, which is not particularly efficient, but with
  ! machine memory so plentiful now, clarity in programming is most important.
  ! R. Boone    Last modified:  September 13, 2014
  use Parameter_Vars

  !********************************************************
  ! Stores attributes of cells comprising the globe
  ! (I went through all the structures, and commented-out those entries that are rarely used, and removed a few
  !  that don't appear at all.  For those that are rarely used, swap them for more common entries.  For example, if biomass
  !  is rarely use, replace it with carbon * 2.5 for leaves and shoots, and carbon * 2.0 for wood.)
  type GlobalCell
    integer    :: zone                    ! A unique ID, from 1 to N, for each cell.  ZONE is a name that comes from ARC GRID.
    ! integer    :: elev                    ! Elevation, in meters ... commented out, not used.
    logical    :: land                    ! Defining land versus ocean
    integer    :: class                   ! Land cover / land use class
    integer    :: landscape_type          ! Landscape type identifying the units for which parameters are given
    logical    :: rangeland               ! A flag showing whether land is rangeland or not
    real       :: latitude                ! The latitude of the center of the cell
    real       :: precip                  ! The total precipitation in the month, in mm
    real       :: max_temp                ! The average maximum temperature in the month, in C
    real       :: min_temp                ! The average minimum temperature in the month, in C
    real       :: precip_average          ! The average annual precipitation, in mm / yr
    real       :: temperature_average     ! The average annual temperature, in C
    real       :: top_sand                ! The percent sand in the top-soil
    real       :: top_silt                ! The percent silt in the top-soil
    real       :: top_clay                ! The percent clay in the top-soil
    real       :: top_gravel              ! The percent rock in the top-soil
    real       :: top_bulk_density        ! The bulk density of the top-soil
    real       :: top_organic_carbon      ! The percent organic matter carbon content in the top-soil
    real       :: sub_sand                ! The percent sand in the sub-soil
    real       :: sub_silt                ! The percent silt in the sub-soil
    real       :: sub_clay                ! The percent clay in the sub-soil
    real       :: sub_gravel              ! The percent rock in the sub-soil
    real       :: sub_bulk_density        ! The bulk density of the sub-soil
    real       :: sub_organic_carbon      ! The percent organic matter carbon content in the sub-soil
    real       :: decid_tree_cover        ! Deciduous tree cover, in percent, perhaps from DeFries, for example, which is from data from 1993-94, but only trees.
    real       :: egreen_tree_cover       ! Evergreen tree cover, in percent, perhaps from DeFries, for example, which is from data from 1993-94, but only trees.
    real       :: shrub_cover             ! DERIVED from other layers using MODIS 44B products and other layers.  VERY poorly known.
    real       :: herb_cover              ! Herbaceous cover, from MODIS 44B product
    real       :: prop_burned             ! The proportion of the cell burned, as shown in maps  (NOTE SCALE DEPENDENCE)
    real       :: prop_fertilized         ! The proportion of the cell fertilized, as shown in maps  (NOTE SCALE DEPENDENCE)
    real       :: temporary               ! A temporary map location, used to read-in values
  end type
  type (GlobalCell) :: Globe(MAX_X_DIM, MAX_Y_DIM)

  !********************************************************
  ! Stores the identifiers of cells that have been shown as RANGELAND
  type RangeCell
    integer    :: x                         ! X dimension of rangeland cell
    integer    :: y                         ! Y dimension of rangeland cell
    integer    :: range_type                ! Identifier storing the type of rangeland cell, used as a key to the Parms strcuture

    real       :: last_month_day_length     ! The day length of the previous month, to know when spring and fall come.
    logical    :: day_length_increasing     ! Increasing or decreasing day length, comparing the current to previous day lengths.
    real       :: day_length                ! Day length, calculated based on latitude and month
    real       :: heat_accumulation         ! Heat accumulation above a base temperature (e.g., 4.4 C in Boone (1999))

    real       :: facet_cover(FACETS)       ! The proportion occupied by each facet
    real       :: total_population(V_LYRS)  ! The total population of each vegetation layer
    real       :: bare_cover                ! Bare cover stored, rather than continually summing the three facets.
    real       :: prop_annual_decid(FACETS) ! Proportion of facet that is annual plants (H_FACET) or deciduous (S_FACET and T_FACET)

    real       :: pot_evap                ! Potential evapotranspiration for the cell (cm/month)
    real       :: evaporation             ! Water evaporated from the soil and vegetation (cm/month)
    real       :: snow                    ! Snowpack, in cm
    real       :: snow_liquid             ! Snowpack liquid water.
    ! Editing the model 01/02/2014 to prevent snow and snow liquid from skyrocketing.   Adding an field for OLD SNOW and ICE, prior to clearing out snow each year.
    real       :: old_snow                ! Snow from past years, including glacial build up.   This will be an accumulator, but essentially outside of active process modeling
    real       :: old_snow_liquid         ! Ditto
                                          ! This won't be output until there is some need for it.
    ! End of addition
    real       :: melt                    ! Snow that melts from snowpack (cm water)
    real       :: pet_remaining           ! Potential evaporation decremented as steps are calculated.  Appears to be a bookkeeping tool.
    real       :: ppt_soil                ! Precipitation adjusted for snow accumulation and melt, and available to infiltrate the soil (cm)
    real       :: runoff                  ! Runoff from the rangeland cell
    real       :: ratio_water_pet         ! Ratio of available water to potential evapotranspiration
!    real       :: co2_value               ! CO2 effect on evapotranspiration, per facet
    real       :: pet_top_soil            ! Potential evaporation from top soil (cm/day)
    real       :: n_leached(SOIL_LAYERS)  ! Nitrogen leached from soil (AMTLEA in Century)
!    real       :: stream(8)               ! Used in stream flow accumulations.  Not sure yet what the 8 cells are for.
    real       :: asmos(SOIL_LAYERS)      ! Used in summing water
    real       :: amov(SOIL_LAYERS)       ! Used in summing water movement
    real       :: storm_flow              ! Storm flow
    real       :: holding_tank            ! Stores water temporarily.  Was asmos(layers+1) in H2OLos
    real       :: transpiration           ! Transpiration water loss
    real       :: relative_water_content(SOIL_LAYERS) ! Used to initialize and during simulation in CENTURY. Here, only during simulation
    real       :: water_available(3)                  ! Water available to plants, available for growth (1), survival (2), and in the two top layers (3)
    real       :: annual_evapotranspiration           ! Annual actual evapotranspiration

    real       :: total_aground_live_biomass          ! Total aboveground live biomass (g/m^2)
    real       :: total_bground_live_biomass          ! Total belowground live biomass (g/m^2)

    real       :: total_litter_carbon(2)              ! Average monthly litter carbon (g/m^2)
    real       :: total_litter_nitrogen(2)            ! Average monthly litter carbon (g/m^2)

!    real       :: non_symbiotic_soil_n_fixing     ! Non-symbiotic soil N fixation (g/m^2)
!    real       :: atmosphere_n_fixing             ! N fixation in atmosphere (g/m^2)
!    real       :: base_n_deposition               ! Base N deposition, updated each year  (g/m^2)

    real       :: root_shoot_ratio(FACETS)         ! Root shoot ratio
    real       :: tree_basal_area                  ! Basal area for trees

    real       :: soil_surface_temperature         ! Average soil surface temperature (C)

                                          ! Soils as in Century 4.5 NLayer=4, 0-15, 15-30, 30-45, 45-60 cm.
                                          ! These will be initialized using approximations and weighted averages from HWSD soils database, which is 0-30 for TOP, 30-100 for SUB.
    real       :: sand(4)                 ! The percent sand in the soil
    real       :: silt(4)                 ! The percent silt in the soil
    real       :: clay(4)                 ! The percent clay in the soil
    ! Some of the following are only used to set field capacities and wilting point.  Perhaps replace with Globe values instead.
!    real       :: gravel(4)               ! The percent rock in the soil
!    real       :: bulk_density(4)         ! The bulk density of the soil
!    real       :: organic_carbon(4)       ! The percent organic matter carbon content in the soil (USED?)
    real       :: mineral_nitrogen(4)     ! Mineral nitrogen content for layer  (g/m2)
!     real, dimension(4) :: soil_depth = (/15.,30.,45.,60./) ! The depth of soils, in cm.   Appears hardwired in some parts of CENTURY, flexible, and up to 9 layers, in other parts of CENTURY.  Likely I noted some values from an early version, but this is a simplification, so...
    real, dimension(4) :: soil_depth = (/15.,15.,15.,15./) ! The depth of soils, in cm.   Appears hardwired in some parts of CENTURY, flexible, and up to 9 layers, in other parts of CENTURY.  Likely I noted some values from an early version, but this is a simplification, so...

    real       :: field_capacity(4)       ! Field capacity for four soils layers shown above.
    real       :: wilting_point(4)        ! Wilting point for four soil layers shown above.
    real       :: soil_total_carbon       ! grams per square meter

    real       :: tree_carbon(WOODY_PARTS)      ! Tree carbon in its components.   These must all be merged or otherwise crosswalked at some point.
    real       :: tree_nitrogen(WOODY_PARTS)    ! Tree nitrogen in its components.   These must all be merged or otherwise crosswalked at some point.
    real       :: shrub_carbon(WOODY_PARTS)      ! Shrub carbon in its components.   These must all be merged or otherwise crosswalked at some point.
    real       :: shrub_nitrogen(WOODY_PARTS)    ! Shrub nitrogen in its components.   These must all be merged or otherwise crosswalked at some point.

    real       :: carbon_nitrogen_ratio(2)        ! Carbon to nitrogen ratio, SURFACE, SOIL
    real       :: fast_soil_carbon(2)             ! Soil organic matter carbon, surface and soil  g/m2  (SOM1C in Century)
    real       :: intermediate_soil_carbon        ! Intermediate soil carbon   g/m2  (SOMC2 in Century)
    real       :: passive_soil_carbon             ! Passive soil carbon    g/m2   (SOMC3 in Century)
    real       :: fast_soil_nitrogen(2)           ! Soil organic matter nitrogen, surface and soil  g/m2  (SOM1E in Century and SSOM1E in Savanna)
    real       :: intermediate_soil_nitrogen      ! Intermediate soil nitrogen   g/m2  (SOM2E in Century)
    real       :: passive_soil_nitrogen           ! Passive soil nitrogen   g/m2   (SOM3E in Century)

    real       :: potential_production                    ! Calculated potential production for the cell, an index.  Based on soil temperature, so not specific to facets.
    real       :: belowground_pot_production(V_LYRS)      ! BIOMASS, Belowground potential production in g/m2
    real       :: aboveground_pot_production(V_LYRS)      ! BIOMASS, Abovegroudn potential production in g/m2
    real       :: total_pot_production(V_LYRS)            ! BIOMASS, Calculate total potential production, in g/m2 with all the corrections in place.
    real       :: co2_effect_on_production(FACETS)        ! Calculated effect of CO2 increasing from 350 to 700 ppm on grassland production, per facet

    real       :: total_pot_prod_limited_by_n(V_LYRS)     ! Coefficient on total potential production reflecting limits due to nitrogen in place  (EPRODL)
    real       :: monthly_net_primary_production          ! Monthly net primary production in g/m2, summed from total_pot_prod_limited_by_n

    real       :: fraction_live_removed_grazing   ! Fraction of live forage removed by grazing  (FLGREM in CENTURY)
    real       :: fraction_dead_removed_grazing   ! Fraction of dead forage removed by grazing  (FDGREM in CENTURY)

    ! Facets are used here.  Facets are: 1 - Herb, 2 - Shrub, 3 - Tree
    ! NOT USED RIGHT NOW:  The array index here is:  1 - Phenological death, 2 - Incremental death, 3 - Herbivory, 4 - Fire
    real       :: temp_effect_on_decomp           ! Temperature effect on decomposition (TFUNC in CENTURY Cycle.f)  (index)
    real       :: water_effect_on_decomp          ! Water effect on decomposition (index)  (Aboveground and belowground entries in CENTURY set to equal, so distinction not made here)
    real       :: anerobic_effect_on_decomp       ! Anerobic effects on decomposition  (index)  (EFFANT in Savanna)
    real       :: all_effects_on_decomp           ! Combined effects on decomposition, which in Savanna includes anerobic (CYCLE.F)  (index)

    real       :: dead_fine_root_carbon(FACETS)       !,4)    ! Dead fine root carbon of the four types cited above.
    real       :: dead_fine_root_nitrogen(FACETS)     !,4)    ! Dead fine root nitrogen of the four types cited above.
    real       :: dead_standing_carbon(FACETS)        !,4)    ! Standing dead carbon of leaf and stem, of the four types cited above.
    real       :: dead_standing_nitrogen(FACETS)      !,4)    ! Standing dead nitrogen of leaf and stem, of the four types cited above.
    real       :: dead_seed_carbon(FACETS)            !,4)    ! Dead seed carbon of the four types cited above.   (gC/m^2)
    real       :: dead_seed_nitrogen(FACETS)          !,4)    ! Dead seed nitrogen of the four types cited above.           (units?)
    real       :: dead_leaf_carbon(FACETS)            !,4)    ! Dead leaf carbon of the four types cited above.   (gC/m^2)
    real       :: dead_leaf_nitrogen(FACETS)          !,4)    ! Dead leaf nitrogen of the four types cited above.
    real       :: dead_fine_branch_carbon(FACETS)     !,4)    ! Dead fine branch carbon of the four types cited above.   (gC/m^2)
    real       :: dead_total_fine_branch_carbon               ! Dead fine branch carbon, summed across facets
    real       :: dead_fine_branch_nitrogen(FACETS)   !,4)    ! Dead fine branch nitrogen of the four types cited above.
    real       :: dead_total_fine_branch_nitrogen             ! Dead fine branch nitrogen, summed across facets
    real       :: dead_coarse_root_carbon(FACETS)     !,4)    ! Dead coarse root carbon of the four types cited above.   (gC/m^2)
    real       :: dead_total_coarse_root_carbon               ! Dead total coarse root carbon, summed across facets
    real       :: dead_coarse_root_nitrogen(FACETS)   !,4)    ! Dead coarse root nitrogen of the four types cited above.
    real       :: dead_total_coarse_root_nitrogen             ! Dead total coarse root nitrogen, summed across facets
    real       :: dead_coarse_branch_carbon(FACETS)     !,4)    ! Dead coarse wood carbon of the four types cited above.   (gC/m^2)
    real       :: dead_total_coarse_branch_carbon               ! Dead total coarse wood carbon, summed across facets
    real       :: dead_coarse_branch_nitrogen(FACETS)   !,4)    ! Dead coarse wood nitrogen of the four types cited above.
    real       :: dead_total_coarse_branch_nitrogen              ! Dead total coarse wood nitrogen, summed across facets

    real       :: lignin_fine_root(FACETS)               ! Fine root lignin concentration
    real       :: lignin_coarse_root(FACETS)             ! Coarse root lignin concentration
    real       :: lignin_fine_branch(FACETS)             ! Fine branch lignin concentration
    real       :: lignin_coarse_branch(FACETS)           ! Coarse branch lignin concentration
    real       :: lignin_leaf(FACETS)                    ! Leaf lignin concentration

    real       :: plant_lignin_fraction(FACETS,2)        ! Lignin in structural residue, at the surface (1) and in the soil (2)  (STRLIG)
    real       :: litter_structural_carbon(2)            ! Litter structural carbon at the surface (1) and in the soil (2)  (STRCIS, or in Savanna, SSTRCIS, with unlabeled and labeled merged)
    real       :: litter_metabolic_carbon(2)             ! Litter metabolic carbon at the surface (1) and in the soil (2)  (METCIS, or in Savanna, SMETCIS)
    real       :: litter_structural_nitrogen(2)          ! Litter structural nitrogen at the surface (1) and in the soil (2)  (STRUCE, or in Savanna, SSTRUCE, with STRUCE named for "elements"  I am only including nitrogen, as in Savanna, so dropping the name)
    real       :: litter_metabolic_nitrogen(2)           ! Litter structural nitrogen at the surface (1) and in the soil (2)  (METABE, or in Savanna, SSTRUCE, with STRUCE named for "elements"  I am only including nitrogen, as in Savanna, so dropping the name)

    ! Temporary storage places, but used across swaths of DECOMP.  If memory is limited, devise an alternative.
    real       :: tnetmin(2)                             ! Temporary storage
    real       :: tminup(2)                              ! Temporary storage
    real       :: grossmin(2)                            ! Temporary storage
    real       :: volitn(2)                              ! Temporary storage, volitized nitrogen
    real       :: fixnit                                 ! Temporary storage, total fixed nitrogen
    real       :: runoffn                                ! Temporary storage, runoff nitrogen
    real       :: e_up(FACETS, WOODY_PARTS)              ! Temporary storage, eup().  Woody parts dimensions for that, but includes ABOVE and BELOW in 1 and 2 for herbaceous material
    real       :: volatized_n                     ! Accumulator for monthy volatilization of N

    ! Growth parameters and others
    real       :: maintain_respiration(FACETS)           ! Maintainence respiration
    real       :: phenology(FACETS)                     ! Phenological stage, a continuous variable from 0 to 4.
    real       :: fine_root_carbon(FACETS)              ! Fine root carbon    (gC/m^2)
    real       :: fine_root_nitrogen(FACETS)            ! Fine root nitrogen  (gN/m^2)
    real       :: seed_carbon(FACETS)                   ! Seed carbon         (gC/m^2)
    real       :: seed_nitrogen(FACETS)                 ! Seed nitrogen       (gN/m^2)
    real       :: leaf_carbon(FACETS)                   ! Leaf carbon         (gC/m^2)
    real       :: leaf_nitrogen(FACETS)                 ! Leaf nitrogen       (gN/m^2)
    real       :: fine_branch_carbon(FACETS)            ! Fine branch carbon  (gC/m^2)
    real       :: fine_branch_nitrogen(FACETS)          ! Fine branch nitrogen(gN/m^2)
    real       :: coarse_root_carbon(FACETS)            ! Coarse root carbon  (gC/m^2)
    real       :: coarse_root_nitrogen(FACETS)          ! Coarse root nitrogen(gN/m^2)
    real       :: coarse_branch_carbon(FACETS)          ! Coarse branch carbon  (gC/m^2)
    real       :: coarse_branch_nitrogen(FACETS)        ! Coarse branch nitrogen(gN/m^2)

    real       :: stored_nitrogen(FACETS)               ! Stored nitrogen (CRPSTG in Century GROWTH, STORAGE in RESTRP).  I can't find where this is initialized, except for a gridded system input.  Assumed 0 for now, but here as a placeholder.
    real       :: plant_nitrogen_fixed(FACETS)          ! Plant nitrogen fixed
    real       :: nitrogen_fixed(FACETS)                ! Nitrogen fixed.  Not sure what components distinguish it, just yet.  (NFIX)

    real       :: respiration_flows(FACETS)             ! Maintenance respiration flows to storage pool  (MRSPSTG)
    real       :: respiration_annual(FACETS)            ! Maintenance respiration flows for year         (MRSPANN)
    real       :: carbon_source_sink                    ! Carbon pool.  (g/m2)   (CSRSNK)    I don't know the utility of this, but incorporating it.
    real       :: nitrogen_source_sink                  ! Nitrogen pool.  (g/m2)   (ESRSNK)

    real       :: carbon_allocation(FACETS, WOODY_PARTS) ! Shrub carbon allocation, by proportion  (TREE_CFAC in Century, except statis here)  Brought into this array, even though it requires more memory, in case I do incorporate dynamic allocation at some point.

    real       :: optimum_leaf_area_index(FACETS)       ! Optimum leaf area index
    real       :: leaf_area_index(FACETS)               ! Leaf area index

    real       :: water_function                        ! Water function influencing mortality (AGWFUNC and BGWFUNC in Century, merged here since CYCLE assigns them equal and they start with the same value)

    real       :: fire_severity                         ! A score from 0 to 1 reflecting fire intensity
    real       :: burned_carbon                         ! The sum of carbon burned, only on the 1 m plots, not whole plant death
    real       :: burned_nitrogen                       ! The sum of nitrogen burned, only on the 1 m plots, not whole plant death

    real       :: fertilized_nitrogen_added             ! Total fertilized nitrogen added  (g/m2)
    real       :: fertilized_carbon_added               ! Total fertilized carbon added  (g/m2)

    integer    :: large_error_count                     ! The count of cells being reset because their values were very very large
    integer    :: neg_error_count                       ! The count of cell being reset because values were below zero

  end type
  type (RangeCell) :: Rng(MAX_RANGE_CELLS)        ! Rng = Rangeland, but "Range" is a reserved word.  In the grids at 0.083 degree resolution, 1,083,507 of those were rangeland.

  !********************************************************
  ! Stores the contents of the SIM_PARM.GRG parameter file
  type SimParm
    logical             :: echo_maps                                  ! Whether to echo spatial maps (1) or not (0)
    character(100)      :: zone_map                                   ! Zonal layer file name
    character(100)      :: land_map                                   ! Land vs. ocean layer file name
    character(100)      :: elev_map                                   ! Elevation layer file name
    character(100)      :: class_map                                  ! Land cover / land use map
    character(30)       :: class_legend                               ! Land cover legend file name, with a rangeland flag stored
    character(100)      :: latitude_map                               ! Map of latitude of the center of each cell
    character(100)      :: landscape_type_map                         ! Map of landscape type
    character(100)      :: top_sand_map                               ! Map of topsoil sand percentage
    character(100)      :: top_silt_map                               ! Map of topsoil silt percentage
    character(100)      :: top_clay_map                               ! Map of topsoil clay percentage
    character(100)      :: top_gravel_map                             ! Map of topsoil gravel volume
    character(100)      :: top_bulk_density_map                       ! Map of topsoil bulk density
    character(100)      :: top_organic_carbon_map                     ! Map of topsoil organic carbon
    character(100)      :: sub_sand_map                               ! Map of sub-soil sand percentage
    character(100)      :: sub_silt_map                               ! Map of sub-soil silt percentage
    character(100)      :: sub_clay_map                               ! Map of sub-soil clay percentage
    character(100)      :: sub_gravel_map                             ! Map of sub-soil gravel volume
    character(100)      :: sub_bulk_density_map                       ! Map of sub-soil bulk density
    character(100)      :: sub_organic_carbon_map                     ! Map of sub-soil organic carbon
    character(100)      :: deciduous_tree_cover_map                   ! Map of deciduous tree cover
    character(100)      :: evergreen_tree_cover_map                   ! Map of evergreen tree cover
    character(100)      :: shrub_cover_map                            ! Map of shrub cover, (e.g., MOD44B as modified in some way to isolate shrubs)
    character(100)      :: herb_cover_map                             ! Map of herbaceous cover (MOD44B)  (BARE GROUND is available if helpful, in MOD44B)
    character(100)      :: precip_average_map                         ! Annual average precipitation map
    character(100)      :: temperature_average_map                    ! Annual average temperature map
    character(140)      :: precip_path_prefix                         ! The path and prefix to the precipitation data
    character(140)      :: max_temp_path_prefix                       ! The path and prefix to the average maximum temperature data
    character(140)      :: min_temp_path_prefix                       ! The path and prefix to the average minimum temperature data
    character(20)       :: precip_temp_suffix                         ! The suffix to the precipitation and temperature data
    character(100)      :: parms_file_name                            ! The name of the file storing landscape unit parameters
    integer             :: fire_maps_used                             ! Whether fire maps will be used (1) or fire will be based on frequencies in the land unit parameter file (0)
    character(140)      :: fire_path                                  ! The path and prefix to any fire maps, if used
    character(100)      :: fire_file_name                             ! The name of the file storing fire history, if used
    integer             :: fertilize_maps_used                        ! Whether fertilize maps will be used (1) or fertilizing will be based on frequencies in the land unit parameter file (0)
    character(140)      :: fertilization_path                         ! The path and prefix to any fertilization maps, if used
    character(100)      :: fertilization_file_name                    ! The name of the file storing fertilization history, if used
    integer             :: start_yr                                   ! Starting year of simulation
    integer             :: end_yr                                     ! Ending year of simulation
    character(100)      :: co2effect_file_name                        ! The name of the file storing the effect of CO2 on production, by year and by facet
    integer             :: echo_level                                 ! A flag representing how much information to echo to the screen ( 0 - After initialization, echo only dates, 1 - Echo weather map names)
    integer             :: state_var_flag                             ! How to treat state variables, 0 = Nothing, 1 = Write-out, 2 = Read, 3 = Read and Write-out
    character(100)      :: state_var_file_out                         ! The output name used for the state variable, if used
    character(100)      :: state_var_file_in                          ! The input name used for the state variable, if used
  end type
  type (SimParm)    :: Sim_Parm

  !********************************************************
  ! Stores parameters unique to each landscape unit

  type UnitParm
    real                :: melting_temp                               ! Melting temperature for snow, in C degrees
    real                :: melting_slope                              ! Slope of melting equation (cm snow / degree C)
    real                :: prcp_threshold                             ! Precipitation required, in cm, before runoff occurs
    real                :: prcp_threshold_fraction                    ! The fraction of precipitation that is runoff
    real                :: base_flow_fraction                         ! The fraction of water classed as base flow
    real                :: soil_transpiration_fraction(SOIL_LAYERS)   ! The water transpired for each depth
    real                :: init_soil_c_n_ratio                        ! Initial soil carbon to nitrogen ratio, from Potter and Klooster (1997) or similar source
    real                :: init_lignin_n_ratio                        ! Initial lignin to nitrogen ratio in litter
    real                :: tree_carbon(WOODY_PARTS,2)                 ! Initial tree carbon g/m2 (RLEAVC, RLVCIS, FBRCIS ... others?  in Century)  ALIVE and DEAD, with LEAVES and FINE ROOTS DEAD NOT USED
    real                :: shrub_carbon(WOODY_PARTS,2)                ! Initial tree carbon g/m2 (RLEAVC, RLVCIS, FBRCIS ... others?  in Century)  ALIVE and DEAD, with LEAVES and FINE ROOTS DEAD NOT USED

    real                :: plant_dimension(FACETS)                    ! Dimension, in meters, of the length or width (i.e., square) of the root volume of plant

    real                :: litter_effect_on_soil_temp                 ! Effect of litter on soil temperature relative to live and standing dead biomass (ELITST in Century)
    real                :: biomass_effect_on_min_soil_temp            ! Effect of biomass on minimum soil surface temperature (PMNTMP in Century)
    real                :: maximum_biomass_soil_temp                  ! Maximum biomass for soil temperature calculations (PMXBIO in Century)
    real                :: biomass_effect_on_max_soil_temp            ! Effect of biomass on maximum soil surface temperature (PMXTMP in Century)

    real                :: ppt_regression_points(3)                   ! Controls on the shape of the regression line connecting available water to PET ratio and plant production.
    real                :: temperature_production(4)                  ! Parameters describing the effect of temperature on potential production.  See the Century Parameterization Workbook for examples.
    real                :: standing_dead_production_halved            ! Level of aboveground standing dead + 10% surface structural C that reduces production by half due to phyiscal obstruction (BIOK5 in CENTURY)
    real                :: radiation_production_coefficient           ! Coefficient for calculating potential aboveground monthly production as a function of solar radiation.  PRDX in CENTURY, with its meaning defined in the CENTURY code, and the online material outdated.
    real                :: fraction_carbon_to_roots(FACETS)           ! Fraction of carbon production allocated to roots (likely from CENTURY 4.  CENTURY 4.5 uses a complex dynamic carbon allocation between roots and shoots)

    integer             :: grazing_effect                             ! A flag from 1 to 6, describing grazing effect responses
    real                :: grazing_effect_multiplier                  ! A multiplier applied to grazing effects 4, 5, and 6

    real                 :: temperature_effect_decomposition(4)        ! Four values, 1) x location of inflection point, 2) y location of inflection point, 3) setp size (distance from the maximum point to the minimum point), 4) slope of line at inflection point.  Default values are:  15.40, 11.75, 29.70, 0.031
    real                 :: anerobic_effect_decomposition(3)           ! Three values, 1) ratio below which there is no effect, 2) ratio below which there is a maximum effect, 3) minimum value of impact of precipitation on anerobic decomposition.
    real                 :: effect_of_co2_on_production(FACETS)        ! Effect of CO2 concentration on production, with 1 meaning no effect, per facet

    real                 :: decomp_rate_structural_litter(2)           ! Decomposition rate of structural litter (per year)				
    real                 :: decomp_rate_metabolic_litter(2)            ! Decomposition rate of metabolic litter (per year)				
    real                 :: decomp_rate_fast_som(2)                    ! Decomposition rate of the fast SOM (soil organic matter) pool (per year)				
    real                 :: decomp_rate_slow_som                       ! Decomposition rate of the slow SOM (soil organic matter) pool (per year)				
    real                 :: decomp_rate_inter_som                      ! Decomposition rate of the intermediate SOM (soil organic matter) pool (per year)				
    real                 :: decomp_rate_fine_branch                    ! Decomposition rate of fine branches (per year)				
    real                 :: decomp_rate_coarse_branch                  ! Decomposition rate of coarse branches (per year)				
    real                 :: decomp_rate_coarse_root                    ! Decomposition rate of coarse roots (per year)				
    real                 :: decomp_rate_structural_litter_inverts(2)   ! Decomposition rate of structural litter in layers 1 and 2 due to invertebrates				
    real                 :: drainage_affecting_anaerobic_decomp        ! Drainage affecting the rate of anaerobic decomposition, spanning from 0 to 1.
    real                 :: feces_lignin                               ! Feces lignin content				
    real                 :: fraction_urine_volatized                   ! Fraction of urine volatilized				
    real                 :: fraction_gross_n_mineral_volatized         ! Fraction of gross N mineral volitalized				
    real                 :: rate_volatization_mineral_n	               ! Rate of volitalization of mineral N				
    real                 :: precip_n_deposition(2)                     ! Parameters relating precipitation to deposition N rate				
!    real                 :: precip_n_symbiotic(2)                      ! Parameters relating precipitation to symbiotic N fixation rate				
    real                 :: decomp_litter_mix_facets                   ! Degree of mixing of litter fall among facets (0 = none, 1 = complete)				
    real                 :: lignin_content_fraction_and_precip(2,2)    ! 1,1 = Intercept, aboveground, 1,2 = Slope, aboveground, 2,1 = Intercept, belowground, 2,2 = Slope, belowground

    real                 :: degree_days_phen(FACETS,10)                ! Degree days and the relationship to plant phenology, by FACET, by 10 values shaping the curve.
    real                 :: degree_days_reset(FACETS)                  ! Total degree days to reset phenology
    real                 :: root_effect_on_nutrients                   ! Effect of root biomass on available nutrients used to determine growth (RICTRL)
    real                 :: root_intercept_on_nutrients                ! Intercept of relationship of effect of root biomass on available nutrients (RIINT)
    real                 :: tree_site_potential                        ! Site potential for trees, which adjusts nitrogen availability in savannas
    real                 :: tree_basal_area_to_grass_nitrogen          ! Correction relating tree basal area to grass nitrogen fraction
    real                 :: tree_basal_area_to_wood_biomass            ! Correction relating tree basal area to grass nitrogen fraction
    real                 :: max_symbiotic_n_fixation_ratio             ! Symbiotic nitrogen fixation maximum for grassland (g N fixed / g C new growth)
    real                 :: fraction_nitrogen_available                ! Fraction of nitrogen available to plants
    real                 :: minimum_c_n_ratio(FACETS,WOODY_PARTS)      ! Minimum carbon/nitrogen ratio (Parts > 2 for grasses will be empty)  (CERCRP set in FLTCE, BUT SET HERE, not doing dynamic carbon as in Century)
    real                 :: maximum_c_n_ratio(FACETS,WOODY_PARTS)      ! Maximum carbon/nitrogen ratio (Parts > 2 for grasses will be empty)  (CERCRP set in FLTCE, BUT SET HERE, not doing dynamic carbon as in Century)
    real                 :: fraction_npp_to_respiration(FACETS)        ! Fraction of net primary production that goes to maintenance respiration
    real                 :: fraction_seeds_not_germinated(FACETS)      ! Fraction of seeds that do not germinate.  This is not used in population dynamics, but rather in decomposition
    real                 :: herb_max_fraction_npp_to_respiration(2)    ! Grass maximum fraction of net primary production that goes to maintenance respiration
    real                 :: woody_max_fraction_npp_to_respiration(5)   ! Woody maximum fraction of net primary production that goes to maintenance respiration
    real                 :: maximum_leaf_area_index                    ! Maximum leaf area index for trees
    real                 :: k_leaf_area_index                          ! I don't know what this is ... ASK.  Not documented on the web or in code.
    real                 :: biomass_to_leaf_area_index_factor          ! Biomass to leaf area index factor
    real                 :: annual_fraction_volatilized_n              ! Annual fraction of nitrogen volatilized
    real                 :: max_herb_root_death_rate                   ! Maximum herbaceous root death rate per month 
    real                 :: fraction_n_absorbed_by_residue(2)          ! Fraction of nitrogen absorbed by residue, for surface (1) and soil (2)
    real                 :: shoot_death_rate(4)                        ! Shoot death rate due to 1) water stress, 2) phenology, 3) shading, according to carbon centration in 4.  (FSDETH)
    real                 :: prop_annuals                               ! The proportion of annual plants in the herbaceous facet.
    real                 :: month_to_remove_annuals                    ! Month to remove annual plants, following their standing dead for some time, contributing to litter, etc.
    real                 :: relative_seed_production(FACETS)           ! Annual seed production, in relative number per year.  Increase one group to favor it.  Increase all groups to increase general establishment.
    real                 :: fraction_aground_npp_to_seeds(FACETS)      ! Fraction of aboveground net primary productivity that goes to seeds, by facets.  For woody plants, it is the proportion of carbon for leaf growth diverted to seeds.
    real                 :: water_effect_on_establish(FACETS,4)        ! Available water:PET ratio effect on establishment, per facet, and with 2 pairs of values used in a regression.
    real                 :: herb_root_effect_on_establish(FACETS,4)    ! Herbaceous root biomass effect on establishment, per facet, and with 2 pairs of values used in a regression.
    real                 :: litter_effect_on_establish(FACETS,4)       ! Litter cover effect on establishment, per facet, and with 2 pairs of values used in a regression.
    real                 :: woody_cover_effect_on_establish(FACETS,4)  ! Woody cover effect on establishment, per facet, and with 2 pairs of values used in a regression.  
    real                 :: nominal_plant_death_rate(FACETS)           ! Nominal plant death rate.  This may be increased by various factors.
    real                 :: water_effect_on_death_rate(FACETS,4)       ! Available water:PET ratio effect on plant death rate, per facet, and with 2 pairs per value used in regression.
    real                 :: grazing_effect_on_death_rate(FACETS,4)     ! Grazing rate effect on plant death rate, per facet, and with 2 pairs per value used in regression.
    real                 :: shading_effect_on_death_rate(FACETS,4)     ! Effect of shading, associated with LAI, on death rate.  Suitable for younger age classes of trees as well (although not explicitly modeled here).
    real                 :: fall_rate_of_standing_dead(FACETS)         ! Rate per month of standing dead to fall to litter
    real                 :: leaf_death_rate(FACETS)                    ! Leaf death rate per month, per facet.
    real                 :: death_rate_of_deciduous_leaves             ! The rate of death of leaves after fall has arrived.
    real                 :: temperature_leaf_out_and_fall(2)           ! Temperature in C at leaf on (1) in spring and leaf fall (2) in fall
    real                 :: drought_deciduous(FACETS)                  ! Whether the deciduous fraction of the plants are typical deciduous (FALSE) or drought deciduous (TRUE)
    real                 :: fraction_woody_leaf_n_translocated         ! The fraction of nitrogen in dead leaves that is translocated into storage.
    real                 :: fine_root_death_rate(FACETS)               ! Death rate of fine root in woody plants, per facet (herbs are a placeholder)
    real                 :: fine_branch_death_rate(FACETS)             ! Death rate of fine branches in woody plants, per facet (herbs are a placeholder)
    real                 :: coarse_branch_death_rate(FACETS)             ! Death rate of coarse wood in woody plants, per facet (herbs are a placeholder)
    real                 :: coarse_root_death_rate(FACETS)             ! Death rate of coarse root in woody plants, per facet (herbs are a placeholder)
    real                 :: fraction_carbon_grazed_returned            ! Fraction of carbon in grazed material that is returned to the system (the rest is in carcasses or milk or the like)
    real                 :: fraction_excreted_nitrogen_in_feces        ! Fraction of nitrogen excreted that is in feces.  The remainder is in urine.
    real                 :: fraction_grazed_by_facet(FACETS)           ! The fraction of grazing that comes from each facet.  Sum to 100%.
    real                 :: fraction_grazed                            ! The annual fraction of plant material that is removed.  This includes both live and dead material.

    real                 :: frequency_of_fire                          ! The probability of fire per year for any given cell within the landscape unit (NOTE SCALE DEPENDENCE, USE DEPENDS ON fire_maps_used), set to 0 for no fire  (unitless)
    real                 :: fraction_burned                            ! The proportion of a landscape cell that burns, in the case of a fire event (NOTE SCALE DEPENDENCE. USE DEPENDS ON fire_maps_used. ALSO ONE FIRE PER YEAR MAX)  (unitless)
    integer              :: burn_month                                 ! The month in which patches will be burned, in the case of a fire event (ONE FIRE PER YEAR MAX, USE DEPENDS ON fire_maps_used)  (month)
    real                 :: fuel_vs_intensity(FIRE_SEVERITIES)                     ! The fuel load as related to low and high intensity fires  (g biomass / m2)
    real                 :: green_vs_intensity(2, FIRE_SEVERITIES)                    ! The proportion of aboveground vegetation that is green versus fire intensity (unitless)
    real                 :: fraction_shoots_burned(FACETS, FIRE_SEVERITIES)        ! The proportion of live leaves and shoots removed by a fire event, by facet, for low and high intensity fire  (unitless)
    real                 :: fraction_standing_dead_burned(FACETS, FIRE_SEVERITIES) ! The proportion of standing dead removed by a fire event, by facet, for low and high intensity fire  (unitless)
    real                 :: fraction_plants_burned_dead(FACETS, FIRE_SEVERITIES)   ! The proportion of plants that are burned that die, by facet, for low and high intensity fire  (unitless)
    real                 :: fraction_litter_burned(FACETS, FIRE_SEVERITIES)        ! The proportion of litter removed by a fire event, by facet, for low and high intensity fire  (unitless)
    real                 :: fraction_burned_carbon_as_ash                          ! The proportion of carbon in burned aboveground material that is ash, going to structural litter  (unitless)
    real                 :: fraction_burned_nitrogen_as_ash                        ! The proportion of nitrogen in burned aboveground material that is ash, going to soil mineral nitrogen  (unitless)
    real                 :: frequency_of_fertilization                 ! The probability of fertlization per year in the landscape unit (USE DEPENDS ON fertilize_maps_used) (unitless)
    real                 :: fraction_fertilized                        ! The proportion of a landscape cell that is fertilized, in the case of a fertilization event  (NOTE SCALE DEPENDENCE.  USE DEPENDS ON fertilize_maps_used)  (unitless)
    integer              :: fertilize_month                            ! The month in which fertilization occurs (one event per year per landscape unit)  (month)
    real                 :: fertilize_nitrogen_added                   ! Amount of inorganic nitrogen added during a fertilization event  (g / m2)
    real                 :: fertilize_carbon_added                     ! Amount of carbon added as part of organic matter fertilizer (g / m2)
    real                 :: fertilize_carbon_nitrogen_ratio            ! Ratio of carbon to nitrogen in organic matter fertilizer added  (unitless)
    ! The following variables are at the landscape unit level, but are calculated, rather than read in.
    integer              :: pot_population(FACETS)                     ! The number of plants that can be supported on 1 x 1 km of land
    real                 :: indiv_plant_area(FACETS)                   ! For brevity, going to store the area of plants
  end type
  type (UnitParm)   :: Parms(2000)        ! Parms will store landscape-specific parameters

  !********************************************************
  ! Stores output variables and a flag as to whether or not they should be written to file
  type OutLayers
    character(60)        :: whole_map_output                           ! Strings identifying
    logical              :: write_out                                  ! Either write-out the entry or do not
  end type
  type (OutLayers)  :: Outs(500)


  !********************************************************
  ! Stores entries in RNG and will tally the number of cell/months where values are reset after exceeding a specified range,
  ! which is usually 0 to V_LARGE.  Some entries in RNG will be skipped, given they cannot exceed the range and still have
  ! the simulation function (e.g., x, y, range_type).  Some others should never prove to be a problem, but they are included
  ! as a check of methods (e.g., day_length).  See Rng() for a description of the entries.
  ! ** 1 = EXCEEDED V_LARGE   2 = WAS NEGATIVE
  type ExceedCnts
    integer       :: last_month_day_length(2)
    integer       :: day_length(2)
    integer       :: heat_accumulation(2)
    integer       :: facet_cover(2)
    integer       :: total_population(2)
    integer       :: bare_cover(2)
    integer       :: prop_annual_decid(2)
    integer       :: pot_evap(2)
    integer       :: evaporation(2)
    integer       :: snow(2)
    integer       :: snow_liquid(2)
    integer       :: melt(2)
    integer       :: pet_remaining(2)
    integer       :: ppt_soil(2)
    integer       :: runoff(2)
    integer       :: ratio_water_pet(2)
    integer       :: co2_value(2)
    integer       :: pet_top_soil(2)
    integer       :: n_leached(2)
    integer       :: asmos(2)
    integer       :: amov(2)
    integer       :: storm_flow(2)
    integer       :: holding_tank(2)
    integer       :: transpiration(2)
    integer       :: relative_water_content(2)
    integer       :: water_available(2)
    integer       :: annual_evapotranspiration(2)
    integer       :: total_aground_live_biomass(2)
    integer       :: total_bground_live_biomass(2)
    integer       :: total_litter_carbon(2)
    integer       :: total_litter_nitrogen(2)
    integer       :: root_shoot_ratio(2)
    integer       :: tree_basal_area(2)
    integer       :: soil_surface_temperature(2)
    integer       :: mineral_nitrogen(2)
    integer       :: field_capacity(2)
    integer       :: wilting_point(2)
    integer       :: soil_total_carbon(2)
    integer       :: tree_carbon(2)
    integer       :: tree_nitrogen(2)
    integer       :: shrub_carbon(2)
    integer       :: shrub_nitrogen(2)
    integer       :: carbon_nitrogen_ratio(2)
    integer       :: structural_carbon(2)
    integer       :: metabolic_carbon(2)
    integer       :: fast_soil_carbon(2)
    integer       :: intermediate_soil_carbon(2)
    integer       :: passive_soil_carbon(2)
    integer       :: fast_soil_nitrogen(2)
    integer       :: intermediate_soil_nitrogen(2)
    integer       :: passive_soil_nitrogen(2)
    integer       :: potential_production(2)
    integer       :: belowground_pot_production(2)
    integer       :: aboveground_pot_production(2)
    integer       :: total_pot_production(2)
    integer       :: co2_effect_on_production(2)
    integer       :: total_pot_prod_limited_by_n(2)
    integer       :: monthly_net_primary_production(2)
    integer       :: fraction_live_removed_grazing(2)
    integer       :: fraction_dead_removed_grazing(2)
    integer       :: temp_effect_on_decomp(2)
    integer       :: water_effect_on_decomp(2)
    integer       :: anerobic_effect_on_decomp(2)
    integer       :: all_effects_on_decomp(2)
    integer       :: dead_fine_root_carbon(2)
    integer       :: dead_fine_root_nitrogen(2)
    integer       :: dead_standing_carbon(2)
    integer       :: dead_standing_nitrogen(2)
    integer       :: dead_seed_carbon(2)
    integer       :: dead_seed_nitrogen(2)
    integer       :: dead_leaf_carbon(2)
    integer       :: dead_leaf_nitrogen(2)
    integer       :: dead_fine_branch_carbon(2)
    integer       :: dead_total_fine_branch_carbon(2)
    integer       :: dead_fine_branch_nitrogen(2)
    integer       :: dead_total_fine_branch_nitrogen(2)
    integer       :: dead_coarse_root_carbon(2)
    integer       :: dead_total_coarse_root_carbon(2)
    integer       :: dead_coarse_root_nitrogen(2)
    integer       :: dead_total_coarse_root_nitrogen(2)
    integer       :: dead_coarse_branch_carbon(2)
    integer       :: dead_total_coarse_branch_carbon(2)
    integer       :: dead_coarse_branch_nitrogen(2)
    integer       :: dead_total_coarse_branch_nitrogen(2)
    integer       :: lignin_fine_root(2)
    integer       :: lignin_coarse_root(2)
    integer       :: lignin_fine_branch(2)
    integer       :: lignin_coarse_branch(2)
    integer       :: lignin_leaf(2)
    integer       :: plant_lignin_fraction(2)
    integer       :: litter_structural_carbon(2)
    integer       :: litter_metabolic_carbon(2)
    integer       :: litter_structural_nitrogen(2)
    integer       :: litter_metabolic_nitrogen(2)
    integer       :: tnetmin(2)
    integer       :: tminup(2)
    integer       :: grossmin(2)
    integer       :: volitn(2)
    integer       :: fixnit(2)
    integer       :: runoffn(2)
    integer       :: e_up(2)
    integer       :: volatized_n(2)
    integer       :: maintain_respiration(2)
    integer       :: phenology(2)
    integer       :: fine_root_carbon(2)
    integer       :: fine_root_nitrogen(2)
    integer       :: seed_carbon(2)
    integer       :: seed_nitrogen(2)
    integer       :: leaf_carbon(2)
    integer       :: leaf_nitrogen(2)
    integer       :: fine_branch_carbon(2)
    integer       :: fine_branch_nitrogen(2)
    integer       :: coarse_root_carbon(2)
    integer       :: coarse_root_nitrogen(2)
    integer       :: coarse_branch_carbon(2)
    integer       :: coarse_branch_nitrogen(2)
    integer       :: stored_nitrogen(2)
    integer       :: plant_nitrogen_fixed(2)
    integer       :: nitrogen_fixed(2)
    integer       :: respiration_flows(2)
    integer       :: respiration_annual(2)
    integer       :: carbon_source_sink(2)
    integer       :: nitrogen_source_sink(2)
    integer       :: carbon_allocation(2)
    integer       :: optimum_leaf_area_index(2)
    integer       :: leaf_area_index(2)
    integer       :: water_function(2)
    integer       :: fire_severity(2)
    integer       :: burned_carbon(2)
    integer       :: burned_nitrogen(2)
  end type
  type (ExceedCnts) :: Exceed

end module
