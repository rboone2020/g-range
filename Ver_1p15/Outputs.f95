subroutine Initialize_Outputs
  !**** Output a file that contains the zonal IDs, X, and Y values for each of the rangeland cells, plus their types.
  !**** This file may be used by a browsing program to plot rangeland maps.  Initially, I had GRange producing global
  !**** maps, but those will be 10x larger than necessary, so now only rangeland cells are produced.
  !****
  !**** R.B. Boone       Last modified:  March 7, 2011
  use Structures
  use Parameter_Vars
  implicit none
  
  integer icell, ix, iy
  character*1  c
  character*140  clear_file
  
  c = ', '
  ! Open an output data file
  open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//OUT_DATA_NAME, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)
  if (ioerr == 0) then
    write(SHORT_USE_FILE,*) range_cells, x_dim, y_dim, lower_x, lower_y, upper_x, upper_y, cellsize, &
      Sim_Parm%start_yr, Sim_Parm%end_yr, &
      '    // Range cells, XDim, YDim, LX, LY, HX, HY, Cell size, First yr, Last Yr. Below: Cell ID, Zone, X, Y, Landscape'
!      '    // Range cells, XDim, YDim, LX, LY, HX, HY, Cell size, First yr, Last Yr. Below: Cell ID, Zone, X, Y, Landscape, Cover'
    do icell = 1, range_cells
      ix = Rng(icell)%x
      iy = Rng(icell)%y
      write(SHORT_USE_FILE,1111) icell,c, Globe(ix,iy)%zone,c, Rng(icell)%x,c, Rng(icell)%y,c, Rng(icell)%range_type 
    end do
1111  format(I10,A,I8,A,I6,A,I6,A,I4,A,I4)
  else
    write(*,*) 'The file describing output cells cannot be opened: ', app_path(1:len_trim(app_path))//OUT_DATA_NAME
    stop
  end if
  close(SHORT_USE_FILE)
  
  ! Outputing a global map, for use in mapping.  
  open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//OUT_MAP_NAME, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)
  if (ioerr == 0) then
    do iy=1,y_dim
      write(SHORT_USE_FILE,1001) (Globe(ix,iy)%land,ix=1,x_dim)
    end do       
1001  format(5000L1)                                    ! This value is hardwired
  else
    write(*,*) 'The file describing output cells cannot be opened: ', app_path(1:len_trim(app_path))//OUT_DATA_NAME
    stop
  end if
  close(SHORT_USE_FILE)
  
  ! Opening and then closing files that were being appended to, to clear-out their contents
  ! May explore STATUS='REPLACE' and STATUS='SCRATCH'
  do ix = 1, WHOLE_MAP_COUNT 
    clear_file = out_path(1:len_trim(out_path))//Outs(ix)%whole_map_output(1:&
                 len_trim(Outs(ix)%whole_map_output))//'.gof'
    if (Outs(ix)%write_out .eq. .TRUE.) then
      open(SHORT_USE_FILE, FILE=clear_file, STATUS='REPLACE')
      close(SHORT_USE_FILE)
    end if
  end do
  
end subroutine


subroutine Output_Surfaces
  !**** Output the surfaces from GRange every month.  Selected surfaces are produced, so that
  !**** users may turn on or off any surface they wish.
  !**** This is a very bulky format, although not inefficient.  I do not know how to refer to elements of a structure
  !**** by their name.  I had a posting on Adobe FORTRAN discussion board for weeks but no response, so apparently it is 
  !**** not trivial.  For now, I am using this very bulky technique, with a set of lines per variable. 
  !
  !**** R.B. Boone       Last modified:  March 11, 2011
  use Structures
  use Parameter_Vars
  implicit none

  real           :: rng_out(MAX_RANGE_CELLS)                        
  integer        :: i, icell, ilayer, ifacet, ipart
  logical        :: write_it

  ! I don't know how to identify an entry in a FORTRAN structure.  If I could, this would be much briefer.  

  ! X
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'x_loc' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = float(Rng(icell)%x)
    end do
    call Output_One_Surface ('x_loc', 5, rng_out)
  end if

  ! Y
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'y_loc' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = float(Rng(icell)%y)
    end do
    call Output_One_Surface ('y_loc', 5, rng_out)
  end if

  ! Rangeland type
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'range_type' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = float(Rng(icell)%range_type)
    end do
    call Output_One_Surface ('range_type', 10, rng_out)
  end if



  ! *********************
  ! *********************
  ! The following three entries are different than all the others.  They output precipitation, maximum, and minimum temperature.
  ! These are the only variables stored in Globe that are output.  
  ! Precipitation
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'precip' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Globe(Rng(icell)%x, Rng(icell)%y)%precip
    end do
    call Output_One_Surface ('precip', 6, rng_out)
  end if

  ! Maximum temperature
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'max_temp' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Globe(Rng(icell)%x, Rng(icell)%y)%max_temp
    end do
    call Output_One_Surface ('max_temp', 8, rng_out)
  end if

  ! Minimum temperature
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'min_temp' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Globe(Rng(icell)%x, Rng(icell)%y)%min_temp
    end do
    call Output_One_Surface ('min_temp', 8, rng_out)
  end if

  ! End of special section
  ! *********************

  ! Day Length
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'day_length' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%day_length
    end do
    call Output_One_Surface ('day_length', 10, rng_out)
  end if

  ! Heat accumulation
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'heat_accumulation' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%heat_accumulation
    end do
    call Output_One_Surface ('heat_accumulation', 17, rng_out)
  end if

  ! Potential evapotranspiration
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'pot_evap' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%pot_evap
    end do
    call Output_One_Surface ('pot_evap', 8, rng_out)
  end if

  ! Evaporation
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'evaporation' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%evaporation
    end do
    call Output_One_Surface ('evaporation', 11, rng_out)
  end if

  ! Snow
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'snow' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%snow
    end do
    call Output_One_Surface ('snow', 4, rng_out)
  end if

  ! Snow liquid
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'snow_liquid' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%snow_liquid
    end do
    call Output_One_Surface ('snow_liquid', 11, rng_out)
  end if

  ! Melt
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'melt' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%melt
    end do
    call Output_One_Surface ('melt', 4, rng_out)
  end if

  ! Potential evapotranspiration remaining
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'pet_remaining' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%pet_remaining
    end do
    call Output_One_Surface ('pet_remaining', 13, rng_out)
  end if

  ! Precipitation in soil
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'ppt_soil' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%ppt_soil
    end do
    call Output_One_Surface ('ppt_soil', 8, rng_out)
  end if

  ! Runoff 
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'runoff' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%runoff
    end do
    call Output_One_Surface ('runoff', 6, rng_out)
  end if
 
  ! Ratio of water to potential evapotranspiration
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'ratio_water_pet' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%ratio_water_pet
    end do
    call Output_One_Surface ('ratio_water_pet', 15, rng_out)
  end if

  ! Potential evapotranspiration from top soil
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'pet_top_soil' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%pet_top_soil
    end do
    call Output_One_Surface ('pet_top_soil', 12, rng_out)
  end if

  ! Nitrogen leached from soil (SOIL_LAYERS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'n_leached' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, SOIL_LAYERS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%n_leached(ilayer)
      end do
      call Output_One_Surface ('n_leached', 9, rng_out)
    end do
  end if

  ! Holding_tank                           
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'holding_tank' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%holding_tank
    end do
    call Output_One_Surface ('holding_tank', 12, rng_out)
  end if

  ! Transpiration
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'transpiration' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%transpiration
    end do
    call Output_One_Surface ('transpiration', 13, rng_out)
  end if

  ! Relative_water_content(SOIL_LAYERS)    
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'relative_water_content' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, SOIL_LAYERS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%relative_water_content(ilayer)
      end do
      call Output_One_Surface ('relative_water_content', 22, rng_out)
    end do
  end if

  ! Water_available(3)                     
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'water_available' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 3
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%water_available(ilayer)
      end do
      call Output_One_Surface ('water_available', 15, rng_out)
    end do
  end if

  ! Annual_evapotranspiration              
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'annual_evapotranspiration' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%annual_evapotranspiration
    end do
    call Output_One_Surface ('annual_evapotranspiration', 25, rng_out)
  end if

  ! Facet_cover(FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'facet_cover' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%facet_cover(ifacet)
      end do
      call Output_One_Surface ('facet_cover', 11, rng_out)
    end do
  end if

  ! Bare_cover
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'bare_cover' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%bare_cover
    end do
    call Output_One_Surface ('bare_cover', 10, rng_out)
  end if

  ! Total population (V_LAYERS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'total_population' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, V_LYRS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%total_population(ilayer)
      end do
      call Output_One_Surface ('total_population', 16, rng_out)
    end do
  end if

  ! Prop_annual_decid(FACETS)              
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'prop_annual_decid' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%prop_annual_decid(ifacet)
      end do
      call Output_One_Surface ('prop_annual_decid', 17, rng_out)
    end do
  end if

  ! total_aground_live_biomass
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'total_aground_live_biomass' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%total_aground_live_biomass
    end do
    call Output_One_Surface ('total_aground_live_biomass', 26, rng_out)
  end if

  ! total_bground_live_biomass
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'total_bground_live_biomass' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%total_bground_live_biomass
    end do
    call Output_One_Surface ('total_bground_live_biomass', 26, rng_out)
  end if

  ! total_litter_carbon(SURFACE and SOIL)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'total_litter_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 2
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%total_litter_carbon(ilayer)
      end do
      call Output_One_Surface ('total_litter_carbon', 19, rng_out)
    end do
  end if

  ! total_litter_nitrogen(SURFACE and SOIL)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'total_litter_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 2
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%total_litter_nitrogen(ilayer)
      end do
      call Output_One_Surface ('total_litter_nitrogen', 21, rng_out)
    end do
  end if

  ! Root_shoot_ratio(FACETS)               
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'root_shoot_ratio' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%root_shoot_ratio(ifacet)
      end do
      call Output_One_Surface ('root_shoot_ratio', 16, rng_out)
    end do
  end if

  ! tree_basal_area
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'tree_basal_area' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%tree_basal_area
    end do
    call Output_One_Surface ('tree_basal_area', 15, rng_out)
  end if

  ! soil_surface_temperature
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'soil_surface_temperature' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%soil_surface_temperature
    end do
    call Output_One_Surface ('soil_surface_temperature', 24, rng_out)
  end if

  ! sand (SOIL_LAYERS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'sand' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, SOIL_LAYERS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%sand(ilayer)
      end do
      call Output_One_Surface ('sand', 4, rng_out)
    end do
  end if

  ! Silt (SOIL_LAYERS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'silt' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, SOIL_LAYERS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%silt(ilayer)
      end do
      call Output_One_Surface ('silt', 4, rng_out)
    end do
  end if

  ! Clay (SOIL_LAYERS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'clay' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, SOIL_LAYERS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%clay(ilayer)
      end do
      call Output_One_Surface ('clay', 4, rng_out)
    end do
  end if

  ! Mineral_nitrogen(SOIL_LAYERS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'mineral_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, SOIL_LAYERS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%mineral_nitrogen(ilayer)
      end do
      call Output_One_Surface ('mineral_nitrogen', 16, rng_out)
    end do
  end if

  ! Field_capacity (SOIL_LAYERS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'field_capacity' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, SOIL_LAYERS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%field_capacity(ilayer)
      end do
      call Output_One_Surface ('field_capacity', 14, rng_out)
    end do
  end if

  ! Wilting_point (SOIL_LAYERS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'wilting_point' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, SOIL_LAYERS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%wilting_point(ilayer)
      end do
      call Output_One_Surface ('wilting_point', 13, rng_out)
    end do
  end if

  ! Soil_total_carbon
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'soil_total_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%soil_total_carbon
    end do
    call Output_One_Surface ('soil_total_carbon', 17, rng_out)
  end if

  ! Tree_carbon (WOODY_PARTS)               
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'tree_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ipart = 1, WOODY_PARTS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%tree_carbon(ipart)
      end do
      call Output_One_Surface ('tree_carbon', 11, rng_out)
    end do
  end if

  ! Tree_nitrogen(WOODY_PARTS)             
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'tree_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ipart = 1, WOODY_PARTS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%tree_nitrogen(ipart)
      end do
      call Output_One_Surface ('tree_nitrogen', 13, rng_out)
    end do
  end if

  ! Shrub_carbon(WOODY_PARTS)              
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'shrub_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ipart = 1, WOODY_PARTS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%shrub_carbon(ilayer)
      end do
      call Output_One_Surface ('shrub_carbon', 12, rng_out)
    end do
  end if

  ! Shrub_nitrogen (WOODY_PARTS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'shrub_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ipart = 1, WOODY_PARTS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%shrub_nitrogen(ipart)
      end do
      call Output_One_Surface ('shrub_nitrogen', 14, rng_out)
    end do
  end if

  ! Carbon_nitrogen_ratio (SURFACE AND SOIL)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'carbon_nitrogen_ratio' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 2
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%carbon_nitrogen_ratio(ilayer)
      end do
      call Output_One_Surface ('carbon_nitrogen_ratio', 21, rng_out)
    end do
  end if

  ! Fast_soil_carbon (SURFACE AND SOIL)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fast_soil_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 2
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%fast_soil_carbon(ilayer)
      end do
      call Output_One_Surface ('fast_soil_carbon', 16, rng_out)
    end do
  end if

  ! Intermediate_soil_carbon 
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'intermediate_soil_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%intermediate_soil_carbon
    end do
    call Output_One_Surface ('intermediate_soil_carbon', 24, rng_out)
  end if

  ! Passive_soil_carbon
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'passive_soil_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%passive_soil_carbon
    end do
    call Output_One_Surface ('passive_soil_carbon', 19, rng_out)
  end if

  ! Fast_soil_nitrogen (SURFACE and SOIL)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fast_soil_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 2
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%fast_soil_nitrogen(ilayer)
      end do
      call Output_One_Surface ('fast_soil_nitrogen', 18, rng_out)
    end do
  end if

  ! Intermediate_soil_nitrogen
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'intermediate_soil_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%intermediate_soil_nitrogen
    end do
    call Output_One_Surface ('intermediate_soil_nitrogen', 26, rng_out)
  end if

  ! Passive_soil_nitrogen
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'passive_soil_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%passive_soil_nitrogen
    end do
    call Output_One_Surface ('passive_soil_nitrogen', 21, rng_out)
  end if

  ! Potential_production (V_LYRS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'potential_production' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%potential_production
    end do
    call Output_One_Surface ('potential_production', 20, rng_out)
  end if

  ! Belowground_pot_production(V_LYRS)     
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'belowground_pot_production' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, V_LYRS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%belowground_pot_production(ilayer)
      end do
      call Output_One_Surface ('belowground_pot_production', 26, rng_out)
    end do
  end if

  ! Aboveground_pot_production(V_LYRS)     
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'aboveground_pot_production' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, V_LYRS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%aboveground_pot_production(ilayer)
      end do
      call Output_One_Surface ('aboveground_pot_production', 26, rng_out)
    end do
  end if

  ! Total_pot_production (V_LYRS)           
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'total_pot_production' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, V_LYRS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%total_pot_production(ilayer)
      end do
      call Output_One_Surface ('total_pot_production', 20, rng_out)
    end do
  end if

  ! CO2_effect_on_production               
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'co2_effect_on_production' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%co2_effect_on_production(ifacet)
      end do
      call Output_One_Surface ('co2_effect_on_production', 24, rng_out)
    end do
  end if

  ! Total_pot_prod_limited_by_n(V_LYRS)    
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'total_pot_prod_limited_by_n' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, V_LYRS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%total_pot_prod_limited_by_n(ilayer)
      end do
      call Output_One_Surface ('total_pot_prod_limited_by_n', 27, rng_out)
    end do                     
  end if

  ! Monthly net primary production
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'monthly_net_primary_prod' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%monthly_net_primary_production
    end do
    call Output_One_Surface ('monthly_net_primary_prod', 24, rng_out)
  end if

  ! Fraction_live_removed_grazing
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fraction_live_removed_grazing' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%fraction_live_removed_grazing
    end do
    call Output_One_Surface ('fraction_live_removed_grazing', 29, rng_out)
  end if

  ! Fraction_dead_removed_grazing
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fraction_dead_removed_grazing' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%fraction_dead_removed_grazing
    end do
    call Output_One_Surface ('fraction_dead_removed_grazing', 29, rng_out)
  end if

  ! Temp_effect_on_decomp
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'temp_effect_on_decomp' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%temp_effect_on_decomp
    end do
    call Output_One_Surface ('temp_effect_on_decomp', 21, rng_out)
  end if

  ! Water_effect_on_decomp
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'water_effect_on_decomp' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%water_effect_on_decomp
    end do
    call Output_One_Surface ('water_effect_on_decomp', 22, rng_out)
  end if

  ! Anerobic_effect_on_decomp
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'anerobic_effect_on_decomp' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%anerobic_effect_on_decomp
    end do
    call Output_One_Surface ('anerobic_effect_on_decomp', 25, rng_out)
  end if

  ! All_effects_on_decomp
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'all_effects_on_decomp' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%all_effects_on_decomp
    end do
    call Output_One_Surface ('all_effects_on_decomp', 21, rng_out)
  end if

  ! Dead_fine_root_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_fine_root_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_fine_root_carbon(ifacet)
      end do
      call Output_One_Surface ('dead_fine_root_carbon', 21, rng_out)
    end do
  end if

  ! Dead_fine_root_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_fine_root_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_fine_root_nitrogen(ifacet)
      end do
      call Output_One_Surface ('dead_fine_root_nitrogen', 23, rng_out)
    end do
  end if

  ! Dead_standing_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_standing_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_standing_carbon(ifacet)
      end do
      call Output_One_Surface ('dead_standing_carbon', 20, rng_out)
    end do
  end if

  ! Dead_standing_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_standing_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_standing_nitrogen(ifacet)
      end do
      call Output_One_Surface ('dead_standing_nitrogen', 22, rng_out)
    end do
  end if

  ! Dead_seed_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_seed_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_seed_carbon(ifacet)
      end do
      call Output_One_Surface ('dead_seed_carbon', 16, rng_out)
    end do
  end if

  ! Dead_seed_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_seed_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_seed_nitrogen(ifacet)
      end do
      call Output_One_Surface ('dead_seed_nitrogen', 18, rng_out)
    end do
  end if

  ! Dead_leaf_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_leaf_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_leaf_carbon(ifacet)
      end do
      call Output_One_Surface ('dead_leaf_carbon', 16, rng_out)
    end do
  end if

  ! Dead_leaf_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_leaf_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_leaf_nitrogen(ifacet)
      end do
      call Output_One_Surface ('dead_leaf_nitrogen', 18, rng_out)
    end do
  end if

  ! Dead_fine_branch_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_fine_branch_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_fine_branch_carbon(ifacet)
      end do
      call Output_One_Surface ('dead_fine_branch_carbon', 23, rng_out)
    end do
  end if

  ! Dead_total_fine_branch_carbon
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_total_fine_branch_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%dead_total_fine_branch_carbon
    end do
    call Output_One_Surface ('dead_total_fine_branch_carbon', 29, rng_out)
  end if

  ! Dead_fine_branch_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_fine_branch_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_fine_branch_nitrogen(ifacet)
      end do
      call Output_One_Surface ('dead_fine_branch_nitrogen', 25, rng_out)
    end do
  end if

  ! Dead_total_fine_branch_nitrogen
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_total_fine_branch_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%dead_total_fine_branch_nitrogen
    end do
    call Output_One_Surface ('dead_total_fine_branch_nitrogen', 31, rng_out)
  end if                      

  ! Dead_coarse_root_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_coarse_root_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_coarse_root_carbon(ifacet)
      end do
      call Output_One_Surface ('dead_coarse_root_carbon', 23, rng_out)
    end do
  end if

  ! Dead_total_coarse_root_carbon
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_total_coarse_root_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%dead_total_coarse_root_carbon
    end do
    call Output_One_Surface ('dead_total_coarse_root_carbon', 29, rng_out)
  end if

  ! Dead_coarse_root_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_coarse_root_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_coarse_root_nitrogen(ifacet)
      end do
      call Output_One_Surface ('dead_coarse_root_nitrogen', 25, rng_out)
    end do
  end if

  ! Dead_total_coarse_root_nitrogen 
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_total_coarse_root_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%dead_total_coarse_root_nitrogen
    end do
    call Output_One_Surface ('dead_total_coarse_root_nitrogen', 31, rng_out)
  end if
  
  ! Dead_coarse_branch_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_coarse_branch_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_coarse_branch_carbon(ifacet)
      end do
      call Output_One_Surface ('dead_coarse_branch_carbon', 25, rng_out)
    end do
  end if

  ! Dead_total_coarse_branch_carbon
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_total_coarse_branch_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%dead_total_coarse_branch_carbon
    end do
    call Output_One_Surface ('dead_total_coarse_branch_carbon', 31, rng_out)
  end if

  ! Dead_coarse_branch_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_coarse_branch_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%dead_coarse_branch_nitrogen(ifacet)
      end do
      call Output_One_Surface ('dead_coarse_branch_nitrogen', 27, rng_out)
    end do
  end if

  ! Dead_total_coarse_branch_nitrogen
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'dead_total_coarse_branch_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%dead_total_coarse_branch_nitrogen
    end do
    call Output_One_Surface ('dead_total_coarse_branch_nitrogen', 33, rng_out)
  end if
  
  ! Lignin_fine_root (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'lignin_fine_root' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%lignin_fine_root(ifacet)
      end do
      call Output_One_Surface ('lignin_fine_root', 16, rng_out)
    end do
  end if

  ! Lignin_coarse_root (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'lignin_coarse_root' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%lignin_coarse_root(ifacet)
      end do
      call Output_One_Surface ('lignin_coarse_root', 18, rng_out)
    end do
  end if

  ! Lignin_fine_branch (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'lignin_fine_branch' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%lignin_fine_branch(ifacet)
      end do
      call Output_One_Surface ('lignin_fine_branch', 18, rng_out)
    end do
  end if

  ! Lignin_coarse_branch (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'lignin_coarse_branch' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%lignin_coarse_branch(ifacet)
      end do
      call Output_One_Surface ('lignin_coarse_branch', 20, rng_out)
    end do
  end if

  ! Lignin_leaf (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'lignin_leaf' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%lignin_leaf(ifacet)
      end do
      call Output_One_Surface ('lignin_leaf', 11, rng_out)
    end do
  end if

!  ! Plant lignin fraction (SURFACE and SOIL)
!  write_it = .FALSE.
!  do i = 1, WHOLE_MAP_COUNT
!    if ( 'plant_lignin_fraction' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
!  end do
!  if ( write_it .eq. .TRUE. ) then
!    do ifacet = 1, FACETS
!      do ilayer = 1, 2
!        do icell = 1, range_cells
!          rng_out(icell) = Rng(icell)%plant_lignin_fraction(ifacet,ilayer)
!        end do
!        call Output_One_Surface ('plant_lignin_fraction', 21, rng_out)
!      end do
!    end do
!  end if

  ! Litter_structural_carbon (SURFACE and SOIL)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'litter_structural_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 2
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%litter_structural_carbon(ilayer)
      end do
      call Output_One_Surface ('litter_structural_carbon', 24, rng_out)
    end do
  end if

  ! Litter_metabolic_carbon (SURFACE and SOIL)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'litter_metabolic_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 2
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%litter_metabolic_carbon(ilayer)
      end do
      call Output_One_Surface ('litter_metabolic_carbon', 23, rng_out)
    end do
  end if

  ! Litter_structural_nitrogen (SURFACE and SOIL)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'litter_structural_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 2
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%litter_structural_nitrogen(ilayer)
      end do
      call Output_One_Surface ('litter_structural_nitrogen', 26, rng_out)
    end do
  end if

  ! Litter_metabolic_nitrogen (SURFACE and SOIL)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'litter_metabolic_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ilayer = 1, 2
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%litter_metabolic_nitrogen(ilayer)
      end do
      call Output_One_Surface ('litter_metabolic_nitrogen', 25, rng_out)
    end do
  end if

  ! Maintain_respiration (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'maintain_respiration' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%maintain_respiration(ifacet)
      end do
      call Output_One_Surface ('maintain_respiration', 20, rng_out)
    end do
  end if

  ! Phenology (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'phenology' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%phenology(ifacet)
      end do
      call Output_One_Surface ('phenology', 9, rng_out)
    end do
  end if

  ! Fine_root_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fine_root_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%fine_root_carbon(ifacet)
      end do
      call Output_One_Surface ('fine_root_carbon', 16, rng_out)
    end do
  end if

  ! Fine_root_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fine_root_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%fine_root_nitrogen(ifacet)
      end do
      call Output_One_Surface ('fine_root_nitrogen', 18, rng_out)
    end do
  end if

  ! Seed_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'seed_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%seed_carbon(ifacet)
      end do
      call Output_One_Surface ('seed_carbon', 11, rng_out)
    end do
  end if

  ! Seed_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'seed_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%seed_nitrogen(ifacet)
      end do
      call Output_One_Surface ('seed_nitrogen', 13, rng_out)
    end do
  end if

  ! Leaf_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'leaf_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%leaf_carbon(ifacet)
      end do
      call Output_One_Surface ('leaf_carbon', 11, rng_out)
    end do
  end if

  ! Leaf_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'leaf_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%leaf_nitrogen(ifacet)
      end do
      call Output_One_Surface ('leaf_nitrogen', 13, rng_out)
    end do
  end if

  ! Fine_branch_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fine_branch_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%fine_branch_carbon(ifacet)
      end do
      call Output_One_Surface ('fine_branch_carbon', 18, rng_out)
    end do
  end if

  ! Fine_branch_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fine_branch_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%fine_branch_nitrogen(ifacet)
      end do
      call Output_One_Surface ('fine_branch_nitrogen', 20, rng_out)
    end do
  end if

  ! Coarse_root_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'coarse_root_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%coarse_root_carbon(ifacet)
      end do
      call Output_One_Surface ('coarse_root_carbon', 18, rng_out)
    end do
  end if

  ! Coarse_root_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'coarse_root_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%coarse_root_nitrogen(ifacet)
      end do
      call Output_One_Surface ('coarse_root_nitrogen', 20, rng_out)
    end do
  end if         

  ! Coarse_branch_carbon (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'coarse_branch_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%coarse_branch_carbon(ifacet)
      end do
      call Output_One_Surface ('coarse_branch_carbon', 20, rng_out)
    end do
  end if

  ! Coarse_branch_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'coarse_branch_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%coarse_branch_nitrogen(ifacet)
      end do
      call Output_One_Surface ('coarse_branch_nitrogen', 22, rng_out)
    end do
  end if

  ! Stored_nitrogen (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'stored_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%stored_nitrogen(ifacet)
      end do
      call Output_One_Surface ('stored_nitrogen', 15, rng_out)
    end do
  end if

  ! Plant_nitrogen_fixed (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'plant_nitrogen_fixed' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%plant_nitrogen_fixed(ifacet)
      end do
      call Output_One_Surface ('plant_nitrogen_fixed', 20, rng_out)
    end do
  end if

  ! Nitrogen_fixed (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'nitrogen_fixed' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%nitrogen_fixed(ifacet)
      end do
      call Output_One_Surface ('nitrogen_fixed', 14, rng_out)
    end do
  end if

  ! Respiration_flows (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'respiration_flows' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%respiration_flows(ifacet)
      end do
      call Output_One_Surface ('respiration_flows', 17, rng_out)
    end do
  end if

  ! Respiration_annual (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'respiration_annual' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%respiration_annual(ifacet)
      end do
      call Output_One_Surface ('respiration_annual', 18, rng_out)
    end do
  end if

  ! Optimum_leaf_area_index (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'optimum_leaf_area_index' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%optimum_leaf_area_index(ifacet)
      end do
      call Output_One_Surface ('optimum_leaf_area_index', 23, rng_out)
    end do
  end if

  ! Leaf_area_index (FACETS)
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'leaf_area_index' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do ifacet = 1, FACETS
      do icell = 1, range_cells
        rng_out(icell) = Rng(icell)%leaf_area_index(ifacet)
      end do
      call Output_One_Surface ('leaf_area_index', 15, rng_out)
    end do
  end if

  ! Water_function
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'water_function' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%water_function
    end do
    call Output_One_Surface ('water_function', 14, rng_out)
  end if

  ! Carbon_source_sink
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'carbon_source_sink' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%carbon_source_sink
    end do
    call Output_One_Surface ('carbon_source_sink', 18, rng_out)
  end if

  ! Nitrogen_source_sink
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'nitrogen_source_sink' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%nitrogen_source_sink
    end do
    call Output_One_Surface ('nitrogen_source_sink', 20, rng_out)
  end if

  ! Fire severity
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fire_severity' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%fire_severity
    end do
    call Output_One_Surface ('fire_severity', 13, rng_out)
  end if


  ! Burned carbon
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'burned_carbon' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%burned_carbon
    end do
    call Output_One_Surface ('burned_carbon', 13, rng_out)
  end if

  ! Burned nitrogen
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'burned_nitrogen' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%burned_nitrogen
    end do
    call Output_One_Surface ('burned_nitrogen', 15, rng_out)
  end if

  ! Fertilized nitrogen added
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fertilized_nitrogen_added' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%fertilized_nitrogen_added
    end do
    call Output_One_Surface ('fertilized_nitrogen_added', 25, rng_out)
  end if
  
  ! Fertilized carbon added
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'fertilized_carbon_added' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%fertilized_carbon_added
    end do
    call Output_One_Surface ('fertilized_carbon_added', 23, rng_out)
  end if

  ! Large error count
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'large_error_count' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%large_error_count
    end do
    call Output_One_Surface ('large_error_count', 17, rng_out)
  end if

  ! Negative error count
  write_it = .FALSE.
  do i = 1, WHOLE_MAP_COUNT
    if ( 'neg_error_count' .eq. Outs(i)%whole_map_output .and. Outs(i)%write_out .eq. .TRUE. ) write_it = .TRUE.
  end do
  if ( write_it .eq. .TRUE. ) then
    do icell = 1, range_cells
      rng_out(icell) = Rng(icell)%neg_error_count
    end do
    call Output_One_Surface ('neg_error_count', 15, rng_out)
  end if

  ! Skipping carbon_allocation, but that could be added.
  
end subroutine




subroutine Output_One_Surface (file_root, n, rng_out)
  !**** Append the monthly results to an output file.  The format of the file is straightforward,
  !**** representing the rangeland cells as long strings.  
  !****
  !**** rng_out stores the contents to be written out.  
  !****
  !**** I am having trouble passing a string with 90 or more trailing blanks.  They are coming in as undefined.  I will incorporate
  !**** n, which will report the length of the file root.
  !****
  !**** R.B. Boone       Last modified:  March 7, 2011
  use Structures
  use Parameter_Vars
  implicit none

  character(100) ::   file_root
  character(140) ::   out_file
  real           ::   rng_out(MAX_RANGE_CELLS)                      
  real           ::   rc, yr, mn
  integer        ::   icell, n

  out_file = out_path(1:len_trim(out_path))//file_root(1:n)//'.gof'
  open(SHORT_USE_FILE, FILE=out_file, ACTION='WRITE', POSITION='APPEND', FORM='BINARY', IOSTAT=ioerr)
  if (ioerr == 0) then
    ! Write-out the month's information -
    ! Writing the number of range cells, year, and month, at the 'top' of every month.  This will act as a check sum.  If the values do not
    ! match when read back in in a display program, something has gone wrong, the range cells do not agree, or more likely,
    ! the programs are out of sync, and what is supposed to be X and Y at a given location will not match.
    rc = range_cells                       ! Forcing retyping for output
    yr = year
    mn = month
    write(SHORT_USE_FILE,ERR=1000) rc, yr, mn
    do icell = 1, range_cells
      write(SHORT_USE_FILE,ERR=1000) rng_out(icell)
    end do
  else
    write(*,*) 'Unable to write an output file, which is: ', out_file
    stop
  end if
  close(SHORT_USE_FILE)

  return

1000  write(*,*) 'Unable to write a record within the output file, which is: ', &
          out_path(1:len_trim(out_path))//file_root(1:n)//'.gof'
  stop

end subroutine