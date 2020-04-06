subroutine Initialize_Parms
!**** Initialize_Parms is the first effort to intialize the model.
!**** This routine uses entries from the UNIX library, which should be widely available, but is non-standard.
!**** R. Boone    Last modified: April 6, 2020   Edited to incorporate bug fixes prior to Github inclusion
   use Parameter_Vars
   use Structures
   implicit none
   integer*4 err_code, getcwd
   integer   cnt_rate, cnt_max
   integer i, temp_out, itemp1
   character(100) a_string
   character(130) b_string
   character(10) c_string
   character(20) d_string          ! Units

   write(*,*)
   write(*,*) '     G-RANGE Ver 1.15  April 6, 2020    '
   write(*,*)
   ! ***************************************************
   ! ***************************************************
   ! Change the following to TRUE to check for -NaN at the end of each procedure.  False will disable the checking.
   check_nan_flag = .TRUE.
   ! Change the following to TRUE to have the model stop when -NaN is found.  If the check-nan-flag is false, this will have no effect.
   stop_on_nan_flag = .TRUE.
   ! ***************************************************
   ! ***************************************************
   call cpu_time(start_time)
   call system_clock(clock_start_time, cnt_rate, cnt_max)
   ! ------------------------------------------------
   ! ----  SETUP THE INITIAL PATHWAYS AND FILES  ----
   ! Get the path to the binary file.  Nothing at all is known about the environment at this point.  This will fill one pathway.
   err_code = getcwd(bin_path)
   if (err_code == 0) then
     write(*,*) 'The path to the binary file is: ', trim(bin_path)
   else
     write(*,*) 'The path to the binary file cannot be determined.'
     stop
   end if
   ! Now that we know the binary path, read a specific existing file to get the application path.
   open(SHORT_USE_FILE, FILE=bin_path(1:len_trim(bin_path))//INI_FILE, ACTION='READ', IOSTAT=ioerr)
   if (ioerr == 0) then
     read(SHORT_USE_FILE,*) app_path
     app_path=trim(app_path)
     write(*,*) 'App path: ',trim(app_path)
     parm_path=app_path(1:len_trim(app_path))//'\Parms\'
     write(*,*) 'Parm path: ',trim(parm_path)
     out_path=app_path(1:len_trim(app_path))//'\Output\'
     write(*,*) 'Output path: ',trim(out_path)
   else
     write(*,*) 'The application path cannot be read from the GRANGE.INI file in the program directory: ', &
     bin_path(1:len_trim(bin_path))//INI_FILE
     stop
   end if
   close(SHORT_USE_FILE)

   ! Open an echo file, which will remain open throughout the simulation, and accept any type of echoed information.
   open(ECHO_FILE, FILE=app_path(1:len_trim(app_path))//ECHO_NAME, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)
   if (ioerr == 0) then
     ! Catching up on echoing variables
     write(ECHO_FILE,*) 'G-RANGE ECHO FILE - Capturing echoed inputs and other information'
     write(ECHO_FILE,*) ' '
     write(ECHO_FILE,*) 'The path to the binary file is: ', trim(bin_path)
     write(ECHO_FILE,*) 'The application pathway is:     ', trim(app_path)
     write(ECHO_FILE,*) 'The pathway to this echo file is:    ', app_path(1:len_trim(app_path))//ECHO_NAME
   else
     write(*,*) 'The echo file cannot be opened: ', app_path(1:len_trim(app_path))//ECHO_NAME
     stop
   end if

   ! Set all the defaults of the shared parameters
   call Set_Defaults

   ! -------------------------------------------
   ! ----  PROCESS THE MAIN PARAMETER FILE  ----
   open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//SIMPARM_NAME, ACTION='READ', IOSTAT=ioerr)
   if (ioerr == 0) then
     read(SHORT_USE_FILE,*) Sim_Parm%echo_maps
     write(ECHO_FILE,*) 'Whether or not to echo maps: ', Sim_Parm%echo_maps
     read(SHORT_USE_FILE,*) Sim_Parm%echo_level
     write(ECHO_FILE,*) 'The level of echoing to screen used: ', Sim_Parm%echo_level
     read(SHORT_USE_FILE,*) a_string                     ! A string that says 'MAP_FILES'
     read(SHORT_USE_FILE,*) Sim_Parm%zone_map
     write(ECHO_FILE,*) 'The zonal map, with unique IDs for each cell: ', trim(Sim_Parm%zone_map)
     read(SHORT_USE_FILE,*) Sim_Parm%land_map
     write(ECHO_FILE,*) 'The land versus ocean map: ', trim(Sim_Parm%land_map)
     read(SHORT_USE_FILE,*) Sim_Parm%elev_map
     write(ECHO_FILE,*) 'The elevation map: ', trim(Sim_Parm%elev_map)
     read(SHORT_USE_FILE,*) Sim_Parm%latitude_map
     write(ECHO_FILE,*) 'The latitude map: ', trim(Sim_Parm%latitude_map)
     read(SHORT_USE_FILE,*) Sim_Parm%top_sand_map
     write(ECHO_FILE,*) 'The topsoil sand map: ', trim(Sim_Parm%top_sand_map)
     read(SHORT_USE_FILE,*) Sim_Parm%top_silt_map
     write(ECHO_FILE,*) 'The topsoil silt map: ', trim(Sim_Parm%top_silt_map)
     read(SHORT_USE_FILE,*) Sim_Parm%top_clay_map
     write(ECHO_FILE,*) 'The topsoil clay map: ', trim(Sim_Parm%top_clay_map)
     read(SHORT_USE_FILE,*) Sim_Parm%top_gravel_map
     write(ECHO_FILE,*) 'The topsoil gravel map: ', trim(Sim_Parm%top_gravel_map)
     read(SHORT_USE_FILE,*) Sim_Parm%top_bulk_density_map
     write(ECHO_FILE,*) 'The topsoil bulk density map: ', trim(Sim_Parm%top_bulk_density_map)
     read(SHORT_USE_FILE,*) Sim_Parm%top_organic_carbon_map
     write(ECHO_FILE,*) 'The topsoil organic carbon map: ', trim(Sim_Parm%top_organic_carbon_map)
     read(SHORT_USE_FILE,*) Sim_Parm%sub_sand_map
     write(ECHO_FILE,*) 'The sub-soil sand map: ', trim(Sim_Parm%sub_sand_map)
     read(SHORT_USE_FILE,*) Sim_Parm%sub_silt_map
     write(ECHO_FILE,*) 'The sub-soil silt map: ', trim(Sim_Parm%sub_silt_map)
     read(SHORT_USE_FILE,*) Sim_Parm%sub_clay_map
     write(ECHO_FILE,*) 'The sub-soil clay map: ', trim(Sim_Parm%sub_clay_map)
     read(SHORT_USE_FILE,*) Sim_Parm%sub_gravel_map
     write(ECHO_FILE,*) 'The sub-soil gravel map: ', trim(Sim_Parm%sub_gravel_map)
     read(SHORT_USE_FILE,*) Sim_Parm%sub_bulk_density_map
     write(ECHO_FILE,*) 'The sub-soil bulk density map: ', trim(Sim_Parm%sub_bulk_density_map)
     read(SHORT_USE_FILE,*) Sim_Parm%sub_organic_carbon_map
     write(ECHO_FILE,*) 'The sub-soil organic carbon map: ', trim(Sim_Parm%sub_organic_carbon_map)
     read(SHORT_USE_FILE,*) Sim_Parm%deciduous_tree_cover_map
     write(ECHO_FILE,*) 'The deciduous tree cover map: ', trim(Sim_Parm%deciduous_tree_cover_map)
     read(SHORT_USE_FILE,*) Sim_Parm%evergreen_tree_cover_map
     write(ECHO_FILE,*) 'The evergreen tree cover map: ', trim(Sim_Parm%evergreen_tree_cover_map)
     read(SHORT_USE_FILE,*) Sim_Parm%shrub_cover_map
     write(ECHO_FILE,*) 'The shrub cover map: ', trim(Sim_Parm%shrub_cover_map)
     read(SHORT_USE_FILE,*) Sim_Parm%herb_cover_map
     write(ECHO_FILE,*) 'The herbaceous cover map: ', trim(Sim_Parm%herb_cover_map)
     read(SHORT_USE_FILE,*) Sim_Parm%landscape_type_map
     write(ECHO_FILE,*) 'The landscape type map: ', trim(Sim_Parm%landscape_type_map)
     read(SHORT_USE_FILE,*) Sim_Parm%class_map
     write(ECHO_FILE,*) 'The land cover and land use map: ', trim(Sim_Parm%class_map)
     read(SHORT_USE_FILE,*) Sim_Parm%class_legend
     write(ECHO_FILE,*) 'The land cover and land use classification file: ', trim(Sim_Parm%class_legend)
     read(SHORT_USE_FILE,*) Sim_Parm%precip_average_map
     write(ECHO_FILE,*) 'The average annual precipitation: ', trim(Sim_Parm%precip_average_map)
     read(SHORT_USE_FILE,*) Sim_Parm%temperature_average_map
     write(ECHO_FILE,*) 'The average annual temperature: ', trim(Sim_Parm%temperature_average_map)
     read(SHORT_USE_FILE,*) Sim_Parm%precip_path_prefix
     write(ECHO_FILE,*) 'The precipitation path and prefix: ', trim(Sim_Parm%precip_path_prefix)
     read(SHORT_USE_FILE,*) Sim_Parm%max_temp_path_prefix
     write(ECHO_FILE,*) 'The maximum temperature path and prefix: ', trim(Sim_Parm%max_temp_path_prefix)
     read(SHORT_USE_FILE,*) Sim_Parm%min_temp_path_prefix
     write(ECHO_FILE,*) 'The minimum temperature path and prefix: ', trim(Sim_Parm%min_temp_path_prefix)
     read(SHORT_USE_FILE,*) Sim_Parm%precip_temp_suffix
     write(ECHO_FILE,*) 'The suffix used for both precipitation and temperature: ', trim(Sim_Parm%precip_temp_suffix)
     read(SHORT_USE_FILE,*) Sim_Parm%parms_file_name
     write(ECHO_FILE,*) 'The landscape unit parameters file name: ', trim(Sim_Parm%parms_file_name)
     read(SHORT_USE_FILE,*) Sim_Parm%fire_maps_used
     write(ECHO_FILE,*) 'Fire maps are to be used (1) or not (0): ', Sim_Parm%fire_maps_used
     read(SHORT_USE_FILE,*) Sim_Parm%fire_path
     write(ECHO_FILE,*) 'The relative pathway to fire maps, if used: ', trim(Sim_Parm%fire_path)
     read(SHORT_USE_FILE,*) Sim_Parm%fire_file_name
     write(ECHO_FILE,*) 'The file name storing mapped fire histories, if used: ', trim(Sim_Parm%fire_file_name)
     read(SHORT_USE_FILE,*) Sim_Parm%fertilize_maps_used
     write(ECHO_FILE,*) 'Fertilize maps are to be used (1) or not (0): ', Sim_Parm%fertilize_maps_used
     read(SHORT_USE_FILE,*) Sim_Parm%fertilization_path
     write(ECHO_FILE,*) 'The relative pathway to fertilization maps, if used: ', trim(Sim_Parm%fertilization_path)
     read(SHORT_USE_FILE,*) Sim_Parm%fertilization_file_name
     write(ECHO_FILE,*) 'The file name storing mapped fertlization histories, if used: ', trim(Sim_Parm%fertilization_file_name)
     read(SHORT_USE_FILE,*) Sim_Parm%co2effect_file_name
     write(ECHO_FILE,*) 'The file name storing the effect of CO2 on production every year, by facet: ', &
                        trim(Sim_Parm%co2effect_file_name)
     read(SHORT_USE_FILE,*) Sim_Parm%state_var_flag
     write(ECHO_FILE,*) 'The state variable flag: ', Sim_Parm%state_var_flag
     read(SHORT_USE_FILE,*) Sim_Parm%state_var_file_out
     write(ECHO_FILE,*) 'The state variable output name: ', trim(Sim_Parm%state_var_file_out)
     read(SHORT_USE_FILE,*) Sim_Parm%state_var_file_in
     write(ECHO_FILE,*) 'The state variable input name: ', trim(Sim_Parm%state_var_file_in)
     read(SHORT_USE_FILE,*) Sim_Parm%start_yr
     write(ECHO_FILE,*) 'Start year for simulation: ', Sim_Parm%start_yr
     read(SHORT_USE_FILE,*) Sim_Parm%end_yr
     write(ECHO_FILE,*) 'End year for simulation: ', Sim_Parm%end_yr
   else
     write(ECHO_FILE,*) 'The main parameter file cannot be opened: ',app_path(1:len_trim(app_path))//SIMPARM_NAME
   end if
   close(SHORT_USE_FILE)

   ! -------------------------------------------
   ! ----  PROCESS THE FILE THAT SHOWS WHICH OUTPUTS TO PRODUCE  ----
   open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//OUT_LIST_NAME, ACTION='READ', IOSTAT=ioerr)
   if (ioerr == 0) then
     do i = 1, WHOLE_MAP_COUNT
       read(SHORT_USE_FILE,*) temp_out, Outs(i)%whole_map_output, itemp1, d_string, c_string, b_string            ! The temporary variable stores array indices, not needed here but used in the browser tool
       if (temp_out .gt. 0) then
         Outs(i)%write_out = .TRUE.
       else
         Outs(i)%write_out = .FALSE.
       end if
       write(ECHO_FILE,*) 'Output for: ', trim(Outs(i)%whole_map_output), '  is (1) yes or (0) no: ', Outs(i)%write_out
     end do
   else
     write(ECHO_FILE,*) 'The output selection file cannot be opened: ',app_path(1:len_trim(app_path))//OUT_LIST_NAME
   end if
   close(SHORT_USE_FILE)

end subroutine


subroutine Initialize_Globe
!**** Initialize_Globe reads the initial spatial data supporting G-Range.
!**** The program will uses GRID ASCII files.
!**** NB: In this module, echoing routines are sometimes customized to specific inputs, but easily changed.
!**** R. Boone    Last modified: July 8, 2010 ... modified to use reading subroutine(s)
   use Parameter_Vars
   use Structures
   implicit none
   integer a_val, ix, iy, i, rowarr(MAX_X_DIM), land_classes, irange, iclass
   character(12) a_string
   character(40) big_string

   ! Read in the ZONAL map.  Set the dimensions of the area as well, given cell size.  Because of that
   ! the Read_Map subroutine is not used here.
   write(*,*) 'Reading in the zonal map.'
   open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//Sim_Parm%zone_map, ACTION='READ', IOSTAT=ioerr)
   if (ioerr == 0) then
     ! Read the header
     read(SHORT_USE_FILE,*) a_string, x_dim
     read(SHORT_USE_FILE,*) a_string, y_dim
     write(ECHO_FILE,*) 'The dimensions of the world are: ', x_dim, ' cells wide by ', y_dim, ' cells high.'
     read(SHORT_USE_FILE,*) a_string, lower_x
     read(SHORT_USE_FILE,*) a_string, lower_y
     read(SHORT_USE_FILE,*) a_string, cellsize
     upper_x = ( cellsize * x_dim ) + lower_x
     upper_y = ( cellsize * y_dim ) + lower_y        ! Lower x and y will be negative, and addition is appropriate
     if (upper_x .ne. 180 .or. upper_y .ne. 90) then
       write(ECHO_FILE,*) 'WARNING: The dimensions of the world do not appear to be correct: ', lower_x, upper_x, lower_y, upper_y
     end if
     read(SHORT_USE_FILE,*) a_string, a_val
     ! Read the data
     do iy=1,y_dim
       read(SHORT_USE_FILE,*) (Globe(ix,iy)%zone, ix=1,x_dim)
     end do
     !! NOTE:  The zonal grids are not ideal for echoing out, as the values get very large.
     !!        If the remaining grids come into the program ok, we'll assume that the zonal
     !!        grid is coming in ok, until evidence suggests otherwise.
   else
     write(ECHO_FILE,*) 'Unable to open the zonal map: ',app_path(1:len_trim(app_path))//Sim_Parm%zone_map
   end if
   close(SHORT_USE_FILE)

   ! Read in the LAND versus ocean map.
   write(*,*) 'Reading in the land and ocean map.'
   call Read_Map ('land ocean  ', Sim_Parm%land_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (int(Globe(ix,iy)%temporary) .ne. 0) then
         Globe(ix,iy)%land = .true.
       else
         Globe(ix,iy)%land = .false.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The Land/Ocean map:'
     do iy=1,y_dim
       write(ECHO_FILE,1001) (Globe(ix,iy)%land,ix=1,x_dim)
     end do
   end if
1001  format(5000L1)                                    ! This value is hardwired

1002  format(5000I1)                                    ! This value is hardwired.


   ! Read in the LAND COVER CLASSIFICATION map.
   write(*,*) 'Reading in the land cover classification map.'
   call Read_Map ('land cover  ', Sim_Parm%class_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary .lt. 0) then
         Globe(ix,iy)%class = 0
       else
         Globe(ix,iy)%class = Globe(ix,iy)%temporary
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The classification map, with classes divided by 10, yielding single integers:'
     do iy=1,y_dim
       do ix=1,x_dim
         if (Globe(ix,iy)%class .gt. 0) then
           rowarr(ix) = Globe(ix,iy)%class / 10
           if (rowarr(ix) > 9) rowarr(ix) = 9
         else
           rowarr(ix) = 0
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   !************************************************************************
   !!!!!! Special section - not reading in a new map here              !!!!!
   !!!!!! Instead, filling a RANGELAND map with boolean values         !!!!!
   !************************************************************************
   write(*,*) 'Storing a map of areas that are rangeland.'
   range_cells = 0
   ! First, need the legend file to be read in.
   open(SHORT_USE_FILE, FILE=parm_path(1:len_trim(parm_path))//Sim_Parm%class_legend, ACTION='READ', IOSTAT=ioerr)
   if (ioerr == 0) then
     read(SHORT_USE_FILE,*) land_classes
     do i=1,land_classes
       read(SHORT_USE_FILE,*) irange, iclass, big_string
       if (irange == 1) then
         rangeland_classes(iclass) = .TRUE.
       else
         rangeland_classes(iclass) = .FALSE.
       end if
     end do
     close(SHORT_USE_FILE)
   else
     write(ECHO_FILE,*) 'Unable to open the rangeland parameter file: ',parm_path(1:len_trim(parm_path))//Sim_Parm%class_legend
   end if
   do i=1,land_classes
     ! Not storing the class names, so a skeletal echo
     write(ECHO_FILE,*) 'Land class (may include values not classed): ', i, '  Rangeland flag: ', rangeland_classes(i)
   end do
   ! Now fill a global map with rangeland values
   do iy=1,y_dim
     do ix=1,x_dim
       i=Globe(ix,iy)%class
       Globe(ix,iy)%rangeland = .FALSE.
       if (i>0) then
         if (rangeland_classes(i) .eq. .TRUE.) then
           Globe(ix,iy)%rangeland = .TRUE.
           range_cells = range_cells + 1
           Rng(range_cells)%x = ix
           Rng(range_cells)%y = iy
         end if
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The rangeland map, with ocean 0, land 1, and rangelands 2:'
     do iy=1,y_dim
       do ix=1,x_dim
         if (Globe(ix,iy)%rangeland .eq. .TRUE.) then
           rowarr(ix) = 2
         else if (Globe(ix,iy)%land .eq. .TRUE.) then
           rowarr(ix) = 1
         else
           rowarr(ix) = 0
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if
   write(ECHO_FILE,*) 'The number of rangeland cells across the globe at this resolution is: ', range_cells

   ! Read in the map showing the LATITUDE of each cell.  This could have been done analytically, but this is clearer
   write(*,*) 'Reading in the latitude map.'
   call Read_Map ('latitude    ', Sim_Parm%latitude_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       Globe(ix,iy)%latitude = Globe(ix,iy)%temporary
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The latitude map, with classes divided by 10, yielding two-character integers, counting sign:'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = Globe(ix,iy)%latitude / 10
       end do
       write(ECHO_FILE,1004) (rowarr(ix),ix=1,x_dim)
     end do
   end if
1004  format(5000I2)                                  ! This value is hardwired.

   ! Read in the topsoil sand percentage
   write(*,*) 'Reading in the topsoil sand percentage.'
   call Read_Map ('Top sand    ', Sim_Parm%top_sand_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%top_sand = Globe(ix,iy)%temporary              ! Percent
       else
         Globe(ix,iy)%top_sand = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The topsoil sand percent map, with values divided by 10'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%top_sand / 10. )
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the topsoil silt percentage
   write(*,*) 'Reading in the topsoil silt percentage.'
   call Read_Map ('Top silt    ', Sim_Parm%top_silt_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%top_silt = Globe(ix,iy)%temporary                     ! Percent
       else
         Globe(ix,iy)%top_silt = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The topsoil silt percent map, with values divided by 10'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%top_silt / 10. )
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the topsoil clay percentage
   write(*,*) 'Reading in the topsoil clay percentage.'
   call Read_Map ('Top clay    ', Sim_Parm%top_clay_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%top_clay = Globe(ix,iy)%temporary               ! Percent
       else
         Globe(ix,iy)%top_clay = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The topsoil clay percent map, with values divided by 10'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%top_clay / 10. )
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the topsoil gravel percent volume
   write(*,*) 'Reading in the topsoil gravel percent volume.'
   call Read_Map ('Top gravel    ', Sim_Parm%top_gravel_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%top_gravel = Globe(ix,iy)%temporary                 ! Percent
       else
         Globe(ix,iy)%top_gravel = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The topsoil gravel percent volume map, with values divided by 10.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%top_gravel / 10. )
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the topsoil bulk density
   write(*,*) 'Reading in the topsoil bulk density.'
   call Read_Map ('Top bulk    ', Sim_Parm%top_bulk_density_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%top_bulk_density = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%top_bulk_density = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The topsoil bulk density map, with values truncated to one-digit integer.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%top_bulk_density )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the topsoil organic carbon content
   write(*,*) 'Reading in the topsoil organic carbon content.'
   call Read_Map ('Top carbon  ', Sim_Parm%top_organic_carbon_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%top_organic_carbon = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%top_organic_carbon = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The topsoil organic carbon content map, with values truncated to one-digit integer.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%top_organic_carbon )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the sub-soil sand percentage
   write(*,*) 'Reading in the sub-soil sand percentage.'
   call Read_Map ('Sub sand    ', Sim_Parm%sub_sand_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%sub_sand = Globe(ix,iy)%temporary                   ! Percent
       else
         Globe(ix,iy)%sub_sand = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The sub-soil sand percent map, with values multiplied divided by 10'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%sub_sand / 10. )
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the sub-soil silt percentage
   write(*,*) 'Reading in the sub-soil silt percentage.'
   call Read_Map ('Sub silt    ', Sim_Parm%sub_silt_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%sub_silt = Globe(ix,iy)%temporary               ! Percent
       else
         Globe(ix,iy)%sub_silt = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The sub-soil silt percent map, with values multiplied by 10.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%sub_silt / 10. )
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the sub-soil clay percentage
   write(*,*) 'Reading in the sub-soil clay percentage.'
   call Read_Map ('Sub clay    ', Sim_Parm%sub_clay_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%sub_clay = Globe(ix,iy)%temporary                    ! Percent
       else
         Globe(ix,iy)%sub_clay = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The sub-soil clay percent map, with values multiplied divided 10.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%sub_clay / 10. )
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the sub-soil gravel percent volume
   write(*,*) 'Reading in the sub-soil gravel percent volume.'
   call Read_Map ('Sub gravel  ', Sim_Parm%sub_gravel_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%sub_gravel = Globe(ix,iy)%temporary                ! Percent
       else
         Globe(ix,iy)%sub_gravel = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The sub-soil gravel percent volume map, with values divided by 10.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%sub_gravel / 10. )
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the sub-soil bulk density
   write(*,*) 'Reading in the sub-soil bulk density.'
   call Read_Map ('Sub bulk    ', Sim_Parm%sub_bulk_density_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%sub_bulk_density = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%sub_bulk_density = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The sub-soil bulk density map, with values truncated to one-digit integer.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%sub_bulk_density )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the sub-soil organic carbon content
   write(*,*) 'Reading in the sub-soil organic carbon content.'
   call Read_Map ('Sub carbon  ', Sim_Parm%sub_organic_carbon_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%sub_organic_carbon = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%sub_organic_carbon = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The sub-soil organic carbon content map, with values truncated to one-digit integer.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%sub_organic_carbon )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if


   ! Read in the landscape type map, the main map identifying the landscape units for which parameters are provided.
   write(*,*) 'Reading in the landscape type map.'
   call Read_Map ('LScape type', Sim_Parm%landscape_type_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%landscape_type = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%landscape_type = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The landscape type map, with values divided by 1.5 to yield a one-digit integer.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%landscape_type / 1.5 )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if


   ! Read in the average annual precipitation
   write(*,*) 'Reading in the average annual precipitation map.'
   call Read_Map ('Avg precip  ', Sim_Parm%precip_average_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%precip_average = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%precip_average = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The average annual precipitation map, with values divided by 1000 to yield a one-digit integer.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%precip_average / 1000. )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if


   ! Read in the average annual temperature
   write(*,*) 'Reading in the average annual temperature map.'
   call Read_Map ('Avg temp.   ', Sim_Parm%temperature_average_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%temperature_average = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%temperature_average = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The average annual temperature map, with values divided by 70 to yield a one-digit integer.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%temperature_average / 70. )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if


   ! Read in the deciduous tree cover map
   write(*,*) 'Reading in the deciduous tree cover map.'
   call Read_Map ('D Tree cover.', Sim_Parm%deciduous_tree_cover_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%decid_tree_cover = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%decid_tree_cover = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The deciduous tree cover map can have values up to 100, so divide by 10 for 0 to 9 response.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%decid_tree_cover / 10. )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if

   ! Read in the evergreen tree cover map
   write(*,*) 'Reading in the evergreen tree cover map.'
   call Read_Map ('E Tree cover.', Sim_Parm%evergreen_tree_cover_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%egreen_tree_cover = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%egreen_tree_cover = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The evergreen tree cover map can have values up to 100, so divide by 10 for 0 to 9 response.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%egreen_tree_cover / 10. )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if


   ! Read in the shrub cover map  (MOD44B as modified in some way to isolate shrubs)
   write(*,*) 'Reading in the shrub cover map.'
   call Read_Map ('Shrb cover ', Sim_Parm%shrub_cover_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%shrub_cover = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%shrub_cover = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The shrub cover map can have values up to 100, so divide by 10 for 0 to 9 response.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%shrub_cover / 10. )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if


   ! Read in the herbaceous/shrub map  (MOD44B)
   write(*,*) 'Reading in the herbaceous/shrub cover map.'
   call Read_Map ('Herb cover ', Sim_Parm%herb_cover_map)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary > 0.) then
         Globe(ix,iy)%herb_cover = Globe(ix,iy)%temporary
       else
         Globe(ix,iy)%herb_cover = 0.
       end if
     end do
   end do
   ! Echo the map
   if (Sim_Parm%echo_maps .eq. .true.) then
     write(ECHO_FILE,*) 'The herbaceous cover map can have values up to 100, so divide by 10 for 0 to 9 response.'
     do iy=1,y_dim
       do ix=1,x_dim
         rowarr(ix) = int( Globe(ix,iy)%herb_cover / 10. )
         if (rowarr(ix) > 9) then
           rowarr(ix) = 9
         end if
       end do
       write(ECHO_FILE,1002) (rowarr(ix),ix=1,x_dim)
     end do
   end if


end subroutine


subroutine Read_Map (map_name, map_file)
!***** Read a map into the TEMPORARY variable that is part of the globe structure.  That will be transfered
!***** to the appropriate location in the calling routine.  Trying to read everything in as REAL from the ASCII
!***** files. The type will be converted in the calling routine.
!***** R. Boone    Last modified: September 24, 2010
   use Parameter_Vars
   use Structures
   implicit none
   character(12)  map_name                    ! A short name to identify the current map in the echo file
   character(100) map_file                    ! The full name to the map, appended to the application path
   integer ix, iy, i
   character(12) a_string
   character(40) big_string

   ! Read in the map.  Set the dimensions of the area as well, given cell size.
   write(ECHO_FILE,*) 'Reading in map: ', map_name, '  file: ', app_path(1:len_trim(app_path))//trim(adjustl(map_file))
   if (Sim_Parm%echo_level .gt. 0) write(*,*) 'Reading in map: ', map_name, '  file: ', &
                                                app_path(1:len_trim(app_path))//trim(adjustl(map_file))
   open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//trim(adjustl(map_file)), ACTION='READ', IOSTAT=ioerr)
   if (ioerr == 0) then
     ! Read just the first two lines of the header, as a double-check
     read(SHORT_USE_FILE,*) a_string, ix
     read(SHORT_USE_FILE,*) a_string, iy
     if (ix .ne. x_dim .or. iy .ne. y_dim) then
       write(ECHO_FILE,*) 'The dimension of the map does not match the zonal map: ', ix, iy
       stop
     end if
     do i=1,4                          ! Discard the remaining header lines
       read(SHORT_USE_FILE,*) big_string
     end do
     ! Read the data
     do iy=1,y_dim
       read(SHORT_USE_FILE,*) (Globe(ix,iy)%temporary, ix=1,x_dim)
     end do
     close(SHORT_USE_FILE)
   else
     write(*,*) 'Unable to open the map: ',app_path(1:len_trim(app_path))//trim(adjustl(map_file))
     write(ECHO_FILE,*) 'Unable to open the map: ',app_path(1:len_trim(app_path))//trim(adjustl(map_file))
     pause
   end if
   ! map_name = '            '
   ! map_file = '                                                                                                                  '

end subroutine


subroutine Initialize_Rangelands
!***** Take additional steps to initialize rangelands.  For example, wilting point and field capacity must
!***** be determined for the rangeland sites.
!***** (Facets were confirmed being reasonably populated through echoed statements)
!*****
!***** R. Boone    Last modified: Oct 10, 2014.

   use Parameter_Vars
   use Structures
   implicit none
   integer icell, ix, iy, ilayer, ifacet, iunit, plant_count, plant_count2
   real precip, temper, silt, clay, organic_carbon(SOIL_LAYERS), bulk_density(SOIL_LAYERS), gravel(SOIL_LAYERS)
   real temp_sum, remaining_c, fraction_metab

   write(*,*) 'Initializing the rangeland cells, which number: ',range_cells

   ! Calculate the soil parameters for the four soil layers, interpolating from the two Harmonized World Soils Database layers.  If new data are used, this will need modified.

   do icell=1,range_cells
     ix = Rng(icell)%x
     iy = Rng(icell)%y
     Rng(icell)%sand(1) = Globe(ix,iy)%top_sand           ! The top and bottom layers get their values directly from the two HWSD layers
     Rng(icell)%sand(4) = Globe(ix,iy)%sub_sand           ! The top and bottom layers get their values directly from the two HWSD layers
     Rng(icell)%sand(2) = ( Globe(ix,iy)%top_sand * 0.6667 ) + ( Globe(ix,iy)%sub_sand * 0.3333 )   ! The other layers get weighted values.
     Rng(icell)%sand(3) = ( Globe(ix,iy)%top_sand * 0.3333 ) + ( Globe(ix,iy)%sub_sand * 0.6667 )   ! The other layers get weighted values.
      Rng(icell)%silt(1) = Globe(ix,iy)%top_silt           ! The top and bottom layers get their values directly from the two HWSD layers
      Rng(icell)%silt(4) = Globe(ix,iy)%sub_silt           ! The top and bottom layers get their values directly from the two HWSD layers
      Rng(icell)%silt(2) = ( Globe(ix,iy)%top_silt * 0.6667 ) + ( Globe(ix,iy)%sub_silt * 0.3333 )   ! The other layers get weighted values.
      Rng(icell)%silt(3) = ( Globe(ix,iy)%top_silt * 0.3333 ) + ( Globe(ix,iy)%sub_silt * 0.6667 )   ! The other layers get weighted values.
     Rng(icell)%clay(1) = Globe(ix,iy)%top_clay           ! The top and bottom layers get their values directly from the two HWSD layers
     Rng(icell)%clay(4) = Globe(ix,iy)%sub_clay           ! The top and bottom layers get their values directly from the two HWSD layers
     Rng(icell)%clay(2) = ( Globe(ix,iy)%top_clay * 0.6667 ) + ( Globe(ix,iy)%sub_clay * 0.3333 )   ! The other layers get weighted values.
     Rng(icell)%clay(3) = ( Globe(ix,iy)%top_clay * 0.3333 ) + ( Globe(ix,iy)%sub_clay * 0.6667 )   ! The other layers get weighted values.
     ! The following are only used a few lines down, so not storing in Rng, for space-saving purposes.
      gravel(1) = Globe(ix,iy)%top_gravel           ! The top and bottom layers get their values directly from the two HWSD layers
      gravel(4) = Globe(ix,iy)%sub_gravel           ! The top and bottom layers get their values directly from the two HWSD layers
      gravel(2) = ( Globe(ix,iy)%top_gravel * 0.6667 ) + ( Globe(ix,iy)%sub_gravel * 0.3333 )   ! The other layers get weighted values.
      gravel(3) = ( Globe(ix,iy)%top_gravel * 0.3333 ) + ( Globe(ix,iy)%sub_gravel * 0.6667 )   ! The other layers get weighted values.
     bulk_density(1) = Globe(ix,iy)%top_bulk_density           ! The top and bottom layers get their values directly from the two HWSD layers
     bulk_density(4) = Globe(ix,iy)%sub_bulk_density           ! The top and bottom layers get their values directly from the two HWSD layers
     bulk_density(2) = ( Globe(ix,iy)%top_bulk_density * 0.6667 ) + ( Globe(ix,iy)%sub_bulk_density * 0.3333 )   ! The other layers get weighted values.
     bulk_density(3) = ( Globe(ix,iy)%top_bulk_density * 0.3333 ) + ( Globe(ix,iy)%sub_bulk_density * 0.6667 )   ! The other layers get weighted values.
      organic_carbon(1) = Globe(ix,iy)%top_organic_carbon           ! The top and bottom layers get their values directly from the two HWSD layers
      organic_carbon(4) = Globe(ix,iy)%sub_organic_carbon           ! The top and bottom layers get their values directly from the two HWSD layers
      organic_carbon(2) = ( Globe(ix,iy)%top_organic_carbon * 0.6667 ) + ( Globe(ix,iy)%sub_organic_carbon * 0.3333 )   ! The other layers get weighted values.
      organic_carbon(3) = ( Globe(ix,iy)%top_organic_carbon * 0.3333 ) + ( Globe(ix,iy)%sub_organic_carbon * 0.6667 )   ! The other layers get weighted values.

     ! Century uses these soil parameters from 0-1, so ...
     do ilayer = 1, SOIL_LAYERS
       Rng(icell)%sand(ilayer) = Rng(icell)%sand(ilayer) / 100.0
       Rng(icell)%silt(ilayer) = Rng(icell)%silt(ilayer) / 100.0
       Rng(icell)%clay(ilayer) = Rng(icell)%clay(ilayer) / 100.0
       gravel(ilayer) = gravel(ilayer) / 100.0
     end do

     ! Calculate the field capacity and wilting point for the rangeland cells.
     ! This process comes from Century 4.5, where they cite Gupta and Larson (1979).
     ! NB: The kg/dm3 for bulk density in the soils database is equal to g/cm3 in Gupta and Larson.
     ! Field capacity is done at a Matric potential of -0.33, as in Century, and includes only option SWFLAG=1, where both wilting point and field capacity are calculated.
     ! Wilting point is done at Matric potential of -15.0.
     ! (Century includes extra components, but the coefficients on those for SWFLAG=1 are 0, so following Gupta and Larson (1979) is correct.)

     ix = Rng(icell)%x
     iy = Rng(icell)%y
     Rng(icell)%range_type = Globe(ix,iy)%landscape_type
     if (Rng(icell)%range_type .lt. 1) then
        write(ECHO_FILE,*) 'A range cell has a landscape type 0.  Make sure GIS layers agree for X and Y: ', ix, iy
        Rng(icell)%range_type = 1
     end if
     iunit = Rng(icell)%range_type
     ! Calculating initial plant populations.  These are based on a 1 km^2 area, and the coverage maps.
     ! Three facets, plus bare ground.
     ! The potential populations of the plants are higher than the aerial coverage of the facets, at least
     ! for herbs and shrubs.  Tree cover and population are the same.  This is due to herbs being in the understory
     ! of shrubs and trees, and shrubs being in the understory of trees.
     do ifacet = 1, FACETS
       ! The following is the total population possible for the facet, if entirely dominated by that facet.
       Parms(iunit)%indiv_plant_area(ifacet) = Parms(iunit)%plant_dimension(ifacet) * Parms(iunit)%plant_dimension(ifacet)             ! m x m = m^2
       Parms(iunit)%pot_population(ifacet) = REF_AREA / Parms(iunit)%indiv_plant_area(ifacet)                                                       ! (m x m) / m^2 = #
     end do

     do ilayer=1,SOIL_LAYERS
       if ( Rng(icell)%sand(ilayer) + Rng(icell)%silt(ilayer) + Rng(icell)%clay(ilayer) + gravel(ilayer) > 0.01 ) then
         Rng(icell)%field_capacity(ilayer) = ( Rng(icell)%sand(ilayer) * 0.3075 ) + ( Rng(icell)%silt(ilayer) * 0.5886 ) + &
                                           ( Rng(icell)%clay(ilayer) * 0.8039 ) + ( organic_carbon(ilayer) * 0.002208 ) + &
                                           ( bulk_density(ilayer) * (-0.14340) )
         Rng(icell)%wilting_point(ilayer) =  ( Rng(icell)%sand(ilayer) * (-0.0059) ) + ( Rng(icell)%silt(ilayer) * 0.1142 ) + &
                                           ( Rng(icell)%clay(ilayer) * 0.5766 ) + ( organic_carbon(ilayer) * 0.002228 ) + &
                                           ( bulk_density(ilayer) * 0.02671 )
       else
         Rng(icell)%field_capacity(ilayer) = 0.03
         Rng(icell)%wilting_point(ilayer) = 0.01
         write(ECHO_FILE,*) 'Warning, check GIS: soil information is not defined for cell: ',icell,' and layer: ',ilayer
!         The following is commented out, to avoid distracting warnings with minor effects on outcomes.  But the error to ECHO.GOF is retained.
!         write(*,*) 'Warning, check GIS: soil information is not defined for cell: ',icell,' and layer: ',ilayer
       end if
       ! Correcting field capacity and wilting point based on gravel volume.
       Rng(icell)%field_capacity(ilayer) = Rng(icell)%field_capacity(ilayer) * ( 1. - gravel(ilayer) )
       Rng(icell)%wilting_point(ilayer) =  Rng(icell)%wilting_point(ilayer) * ( 1. - gravel(ilayer) )

       Rng(icell)%relative_water_content(ilayer) = 0.50
       ! Initialize asmos to the range between capacity and wilting, plus the bottom value, wilting.
       ! Then multiply that by the relative water content, and finally soil depth.  The other measures are for 1 cm deep soil, essentially.
       Rng(icell)%asmos(ilayer) = ( ( Rng(icell)%field_capacity(ilayer) - Rng(icell)%wilting_point(ilayer) ) * &
              Rng(icell)%relative_water_content(ilayer) + Rng(icell)%wilting_point(ilayer) ) * Rng(icell)%soil_depth(ilayer)
     end do
   end do

   ! Calculate total carbon in the soil.  This uses average temperature and precipitation, which Century
   ! truncates to fairly low values, and so here we do the same. NOTE:  How incorrect it is to initialize
   ! forested soils using the grassland initialization is a question.  But we won't be simulating forests per sey.
   do icell=1,range_cells
     ix = Rng(icell)%x
     iy = Rng(icell)%y
     iunit = Rng(icell)%range_type
     temper = Globe(ix,iy)%temperature_average
     precip = Globe(ix,iy)%precip_average
     silt = ( Rng(icell)%silt(1) + Rng(icell)%silt(2) + Rng(icell)%silt(3) + Rng(icell)%silt(4) ) / 4.0
     clay = ( Rng(icell)%clay(1) + Rng(icell)%clay(2) + Rng(icell)%clay(3) + Rng(icell)%clay(4) ) / 4.0
     if (temper .gt. 23.) temper = 23.
     if (precip .gt. 120.) precip = 120.
     ! Initialize total soil carbon in grams using the formula in Century, which combines som1c, som2c, and som3c
     Rng(icell)%soil_total_carbon = (-8.27E-01 * temper + 2.24E-02 * temper * temper + precip * 1.27E-01 - 9.38E-04 &
            * precip * precip + precip * silt * 8.99E-02 + precip * clay * 6.00E-02 + 4.09) * 1000.
     ! Truncated as in Century.  Not allowed to go below 500 g/m^2
     if ( Rng(icell)%soil_total_carbon .lt. 500. ) Rng(icell)%soil_total_carbon = 500.
     Rng(icell)%carbon_nitrogen_ratio(SURFACE_INDEX) = Parms(iunit)%init_soil_c_n_ratio
     Rng(icell)%carbon_nitrogen_ratio(SOIL_INDEX) = Parms(iunit)%init_soil_c_n_ratio

     ! Century cites equations by Burke to initialize carbon pool compartments
!     Rng(icell)%fast_soil_carbon(SURFACE_INDEX) = &
!                                   Rng(icell)%soil_total_carbon * 0.02 + ( ( Rng(icell)%soil_total_carbon * 0.02 ) * 0.011 )
     Rng(icell)%fast_soil_carbon(SURFACE_INDEX) = 10.0 + ( 10.0 * 0.011 )
     Rng(icell)%fast_soil_carbon(SOIL_INDEX) = ( Rng(icell)%soil_total_carbon * 0.02 ) + &
                                               ( ( Rng(icell)%soil_total_carbon * 0.02 ) * 0.011 )
     Rng(icell)%fast_soil_nitrogen(SURFACE_INDEX) = Rng(icell)%fast_soil_carbon(SURFACE_INDEX) * &
                                                    ( 1.0 / Rng(icell)%carbon_nitrogen_ratio(SURFACE_INDEX) )
     Rng(icell)%fast_soil_nitrogen(SOIL_INDEX) =  Rng(icell)%fast_soil_carbon(SOIL_INDEX) * &
                                                    ( 1.0 / Rng(icell)%carbon_nitrogen_ratio(SOIL_INDEX) )
     Rng(icell)%intermediate_soil_carbon = Rng(icell)%soil_total_carbon * 0.64 + &
                                           ( ( Rng(icell)%soil_total_carbon * 0.64 ) * 0.011)
     Rng(icell)%passive_soil_carbon = Rng(icell)%soil_total_carbon * 0.34 + ( ( Rng(icell)%soil_total_carbon * 0.34 ) * 0.011 )
     Rng(icell)%intermediate_soil_nitrogen = Rng(icell)%intermediate_soil_carbon * ( 1.0 / Parms(iunit)%init_soil_c_n_ratio )
     Rng(icell)%passive_soil_nitrogen = Rng(icell)%passive_soil_carbon * ( 1.0 / Parms(iunit)%init_soil_c_n_ratio )

     ! Surface fast soil carbon was removed from the following, since it is assigned 10.11 by default.
     remaining_c = Rng(icell)%soil_total_carbon - Rng(icell)%fast_soil_carbon(SOIL_INDEX) - &
                   Rng(icell)%intermediate_soil_carbon - Rng(icell)%passive_soil_carbon
     ! See multiple cites for the following, including Parten et al. (1993)
     fraction_metab = 0.85 - ( 0.018 * Parms(iunit)%init_lignin_n_ratio )
     fraction_metab = max(0.2, fraction_metab)
     ! Assigning initial carbon and nitrogen concentrations
     do ilayer=SURFACE_INDEX, SOIL_INDEX
       ! Values differ in Century parameter files.   100 appears typical, and spin-up should customize responses
       Rng(icell)%litter_structural_carbon(ilayer) = 100.
       Rng(icell)%litter_metabolic_carbon(ilayer) = 100. * fraction_metab
       Rng(icell)%litter_structural_nitrogen(ilayer) = Rng(icell)%litter_structural_carbon(ilayer) * &
                                                       Parms(iunit)%init_lignin_n_ratio
       Rng(icell)%litter_metabolic_nitrogen(ilayer) = Rng(icell)%litter_structural_nitrogen(ilayer) * fraction_metab
     end do

     do ifacet = 1, FACETS
       Rng(icell)%leaf_carbon(ifacet) = 200. + ( 200. * 0.011 )
       Rng(icell)%leaf_nitrogen(ifacet) = 3.0
       Rng(icell)%dead_standing_carbon(ifacet) = 80. + ( 80. * 0.011 )
       Rng(icell)%dead_standing_nitrogen(ifacet) = 1.6
       Rng(icell)%fine_root_carbon(ifacet) = 200. + ( 200. * 0.011 )
       Rng(icell)%fine_root_nitrogen(ifacet) = 3.0
       Rng(icell)%dead_fine_branch_carbon(ifacet) = 0.0                                      ! Setting to 0, but really just for herbs.
       Rng(icell)%dead_coarse_branch_carbon(ifacet) = 0.0
       Rng(icell)%dead_coarse_root_carbon(ifacet) = 0.0
       Rng(icell)%root_shoot_ratio(ifacet) = Rng(icell)%leaf_carbon(ifacet) / Rng(icell)%fine_root_carbon(ifacet)
     end do
     ! Standardize the three surfaces in case they sum to greater than 1.0  (they are allowed to be less than 1.0, with the remainder being bare ground)
     temp_sum = Globe(Rng(icell)%x, Rng(icell)%y)%herb_cover + Globe(Rng(icell)%x, Rng(icell)%y)%shrub_cover + &
                Globe(Rng(icell)%x, Rng(icell)%y)%decid_tree_cover + Globe(Rng(icell)%x, Rng(icell)%y)%egreen_tree_cover
     if (temp_sum .gt. 100.0) then
       Globe(Rng(icell)%x, Rng(icell)%y)%herb_cover = Globe(Rng(icell)%x, Rng(icell)%y)%herb_cover * (100.0 / temp_sum )
       Globe(Rng(icell)%x, Rng(icell)%y)%shrub_cover = Globe(Rng(icell)%x, Rng(icell)%y)%shrub_cover * (100.0 / temp_sum )
       Globe(Rng(icell)%x, Rng(icell)%y)%decid_tree_cover = Globe(Rng(icell)%x, Rng(icell)%y)%decid_tree_cover * (100.0/temp_sum)
       Globe(Rng(icell)%x, Rng(icell)%y)%egreen_tree_cover = Globe(Rng(icell)%x, Rng(icell)%y)%egreen_tree_cover * (100.0/temp_sum)
     end if

     ! Facet_cover is the straight proportion of each facet on the 1 km^2.  Facet_population includes understory plants.
     Rng(icell)%facet_cover(T_FACET) = ( Globe(Rng(icell)%x, Rng(icell)%y)%decid_tree_cover + &
                                         Globe(Rng(icell)%x, Rng(icell)%y)%egreen_tree_cover ) / 100.
     if (Rng(icell)%facet_cover(T_FACET) .gt. 0.0001) then
       Rng(icell)%prop_annual_decid(T_FACET) = Globe(Rng(icell)%x, Rng(icell)%y)%decid_tree_cover / &
            ( Globe(Rng(icell)%x, Rng(icell)%y)%decid_tree_cover + Globe(Rng(icell)%x, Rng(icell)%y)%egreen_tree_cover )
     else
       Rng(icell)%prop_annual_decid(T_FACET) = 0.0
     end if
     ! Shrub cover, which has no good surface to define it (confirmed by Dr. Hansen himself)
     Rng(icell)%facet_cover(S_FACET) = Globe(Rng(icell)%x, Rng(icell)%y)%shrub_cover / 100.
     if (Rng(icell)%facet_cover(S_FACET) .gt. 0.99) then
       Rng(icell)%facet_cover(S_FACET) = 0.99                              ! Trim any cell that is 100% shrubs to allow some herbs
     end if
     ! THE FOLLOWING COULD BE A PARAMETER.  For now, setting shrub deciduous proporation equal to tree deciduous portion, which should capture large-scale biome variation.
     Rng(icell)%prop_annual_decid(S_FACET) = Rng(icell)%prop_annual_decid(T_FACET)
     ! NOTE: Putting 1% cover into each cell, as an initial value only.
     Rng(icell)%facet_cover(H_FACET) = Globe(Rng(icell)%x, Rng(icell)%y)%herb_cover / 100.
     if (Rng(icell)%facet_cover(H_FACET) .lt. 0.01) then
       Rng(icell)%facet_cover(H_FACET) = 0.01
     end if
     ! The following is a parameter ...
     Rng(icell)%prop_annual_decid(H_FACET) = Parms(iunit)%prop_annuals
     ! Shrubs are the largest unknown, so if there is a problem, subtract from shrubs
     if ((Rng(icell)%facet_cover(T_FACET) + Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(H_FACET) .gt. 1.0)) then
       Rng(icell)%facet_cover(S_FACET) = 1.0 - (Rng(icell)%facet_cover(T_FACET) + Rng(icell)%facet_cover(H_FACET))
     end if
     Rng(icell)%bare_cover = (1.0 - ( Rng(icell)%facet_cover(T_FACET) + Rng(icell)%facet_cover(S_FACET) + &
                             Rng(icell)%facet_cover(H_FACET) ) )
!     write(ECHO_FILE,'(A10,I6,5(F7.4,2X))') 'FACETS: ',icell,Rng(icell)%facet_cover(T_FACET),Rng(icell)%facet_cover(S_FACET), &
!                         Rng(icell)%facet_cover(H_FACET), Rng(icell)%bare_cover, (Rng(icell)%facet_cover(T_FACET) + &
!                         Rng(icell)%facet_cover(S_FACET) + Rng(icell)%facet_cover(H_FACET) + Rng(icell)%bare_cover)

     ! Calculate facet populations.   ** Initializing herbs in 1/3 understory of shrubs.  Shrubs in 1/3 understory of trees, and
     !                                   herbs in 1/6 understory of trees and shrubs **    CONSIDER THIS AND ITS REPERCUSSIONS
     !                                   Ok for initialization, but in the simulation, woody cover.
     ! TREES
     Rng(icell)%total_population(T_LYR) = Rng(icell)%facet_cover(T_FACET) * Parms(iunit)%pot_population(T_FACET)   ! Tree population and cover are directly related.
     ! SHRUBS
     ! Calculate number of shrubs under trees, if 1/3 what could be fitted with full packing (trunks etc. are ignored here)
     plant_count = ( Rng(icell)%total_population(T_LYR) * Parms(iunit)%indiv_plant_area(T_FACET) ) / &
                     Parms(iunit)%indiv_plant_area(S_FACET)
     Rng(icell)%total_population(S_T_LYR) = plant_count * 0.3334
     ! Calculate total shrubs on shrub facet
     Rng(icell)%total_population(S_LYR) = Rng(icell)%facet_cover(S_FACET) * Parms(iunit)%pot_population(S_FACET)
     ! HERBS
     ! Calculate number of herbs under trees, if 1/6 what could be fitted with full packing (trunks etc. are ignored here)
     plant_count = ( Rng(icell)%total_population(T_LYR) * Parms(iunit)%indiv_plant_area(T_FACET) ) / &
                     Parms(iunit)%indiv_plant_area(H_FACET)
     Rng(icell)%total_population(H_T_LYR) = plant_count * 0.16667
     ! Calculate number of herbs under shrubs, if 1/3 what could be fitted with full packing (trunks etc. are ignored here)
     plant_count2 = ( Rng(icell)%total_population(S_LYR) * Parms(iunit)%indiv_plant_area(S_FACET) ) / &
                      Parms(iunit)%indiv_plant_area(H_FACET)
     Rng(icell)%total_population(H_S_LYR) = plant_count2 * 0.3334
     ! Calculate total herbs on herb facet, and sum those on tree and shrub facet to yield the total number of herbs
     Rng(icell)%total_population(H_LYR) = Rng(icell)%facet_cover(H_FACET) * Parms(iunit)%pot_population(H_FACET)

     ! Initialize lignin structural residue
     do ifacet = 1, FACETS
       Rng(icell)%plant_lignin_fraction(ifacet, SURFACE_INDEX) = 0.25
       Rng(icell)%plant_lignin_fraction(ifacet, SOIL_INDEX) = 0.25                                      ! Alter what these should be initialized to, as needed
     end do
   end do

   !***** Handle State Variable Request
   if (Sim_Parm%state_var_flag .eq. 2 .or. Sim_Parm%state_var_flag .eq. 3) then
     write(*,*) ' '
     write(*,*) 'Reading the state of the model from: ', app_path(1:len_trim(app_path))//Sim_Parm%state_var_file_in
     call Read_State
   end if

end subroutine


subroutine Initialize_Landscape_Parms
!**** Initialize_Landscape_Parms reads in the parameters that are associated with the landscape units to be simulated.
!**** It fills the Parms structure, which is dimensioned to be quite large (2000 units).
!**** The parameters will be stored in one large file, but if awkward, that may be split, or more likely, its creation automated,
!**** so that a separate program is used to set the values for a landscape unit, then the results are written to one file.
!****
!**** R. Boone    Last modified: July 14, 2014
   use Parameter_Vars
   use Structures
   implicit none
   integer iunit, ifacet
   integer i, j

   ! NOTE:  **** THIS IS NOT YET PROTECTED FROM ERRORS - ASSUMES A GOOD FILE STRUCTURE ****

   ! -------------------------------------------
   ! ----  PROCESS THE LANDSCAPE UNIT PARAMETER FILE  ----
   ! ----  A negative or zero entry for landscape unit number will stop reading of the file
   ! ----  Select reads will be commented out, where the parameter values are better as fixed values
   ! ----  This routine must be kept in sync with the Structures Parms type, and with the parameter file (entries, order, etc.)
   open(SHORT_USE_FILE, FILE=app_path(1:len_trim(app_path))//'/Parms/'//Sim_Parm%parms_file_name, ACTION='READ', IOSTAT=ioerr)
   if (ioerr == 0) then
     iunit = 1
     do while (iunit .gt. 0)
       read(SHORT_USE_FILE,*) iunit
       write(ECHO_FILE,*) 'READING IN IUNIT: ',iunit
       if (iunit .gt. 0) then
         Parms(iunit)%melting_temp = 0.00000  ! read(SHORT_USE_FILE,*) Parms(iunit)%melting_temp
         write(ECHO_FILE,*) 'Melting temperature: ',  Parms(iunit)%melting_temp
         Parms(iunit)%melting_slope = 0.00200 ! read(SHORT_USE_FILE,*) Parms(iunit)%melting_slope
         write(ECHO_FILE,*) 'Slope of line in melting relationship: ',  Parms(iunit)%melting_slope
         read(SHORT_USE_FILE,*) Parms(iunit)%prcp_threshold
         write(ECHO_FILE,*) 'Precipitation threshold: ',  Parms(iunit)%prcp_threshold
         read(SHORT_USE_FILE,*) Parms(iunit)%prcp_threshold_fraction
         write(ECHO_FILE,*) 'Precipitation threshold fraction: ',  Parms(iunit)%prcp_threshold_fraction
         read(SHORT_USE_FILE,*) Parms(iunit)%base_flow_fraction
         write(ECHO_FILE,*) 'Base flow fraction: ',  Parms(iunit)%base_flow_fraction
         read(SHORT_USE_FILE,*) (Parms(iunit)%soil_transpiration_fraction(i),i=1,SOIL_LAYERS)
         write(ECHO_FILE,*) 'Soil transpiration fraction, by layer: ',(Parms(iunit)%soil_transpiration_fraction(i),i=1,SOIL_LAYERS)
         read(SHORT_USE_FILE,*) Parms(iunit)%init_soil_c_n_ratio
         write(ECHO_FILE,*) 'Initial soil carbon nitrogen ratio: ',  Parms(iunit)%init_soil_c_n_ratio
         read(SHORT_USE_FILE,*) Parms(iunit)%init_lignin_n_ratio
         write(ECHO_FILE,*) 'Initial lignin to nitrogen ratio: ',  Parms(iunit)%init_lignin_n_ratio
         read(SHORT_USE_FILE,*) ((Parms(iunit)%shrub_carbon(i,j),i=1,WOODY_PARTS),j=1,2)
         write(ECHO_FILE,*) 'Shrub carbon, by part, by alive and dead: ',  ((Parms(iunit)%shrub_carbon(i,j),j=1,2),i=1,WOODY_PARTS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%tree_carbon(i,j),i=1,WOODY_PARTS),j=1,2)
         write(ECHO_FILE,*) 'Tree carbon, by part, by alive and dead: ',  ((Parms(iunit)%tree_carbon(i,j),j=1,2),i=1,WOODY_PARTS)
         ! read(SHORT_USE_FILE,*) (Parms(iunit)%plant_dimension(i),i=1,FACETS)
         Parms(iunit)%plant_dimension(1) = 0.5
         Parms(iunit)%plant_dimension(2) = 2.0
         Parms(iunit)%plant_dimension(3) = 8.0
         write(ECHO_FILE,*) 'Plant dimension, by facet: ',  (Parms(iunit)%plant_dimension(i),i=1,FACETS)
         Parms(iunit)%litter_effect_on_soil_temp = 0.4   ! read(SHORT_USE_FILE,*) Parms(iunit)%litter_effect_on_soil_temp
         write(ECHO_FILE,*) 'Litter effect on soil temperature: ',  Parms(iunit)%litter_effect_on_soil_temp
         Parms(iunit)%biomass_effect_on_min_soil_temp = 0.004    ! read(SHORT_USE_FILE,*) Parms(iunit)%biomass_effect_on_min_soil_temp
         write(ECHO_FILE,*) 'Biomass effect on minimum soil temperature: ',  Parms(iunit)%biomass_effect_on_min_soil_temp
         Parms(iunit)%maximum_biomass_soil_temp = 600.0   ! read(SHORT_USE_FILE,*) Parms(iunit)%maximum_biomass_soil_temp
         write(ECHO_FILE,*) 'Maximum biomass affecting soil temperature: ',  Parms(iunit)%maximum_biomass_soil_temp
         Parms(iunit)%biomass_effect_on_max_soil_temp = -0.00350   ! read(SHORT_USE_FILE,*) Parms(iunit)%biomass_effect_on_max_soil_temp
         write(ECHO_FILE,*) 'Biomass effect on maximum soil temperature: ',  Parms(iunit)%biomass_effect_on_max_soil_temp
         Parms(iunit)%ppt_regression_points(1) = 0.0   ! read(SHORT_USE_FILE,*) (Parms(iunit)%ppt_regression_points(i),i=1,3)
         Parms(iunit)%ppt_regression_points(2) = 1.0
         Parms(iunit)%ppt_regression_points(3) = 0.8
         write(ECHO_FILE,*) 'Precipitation regression points, of three: ',  (Parms(iunit)%ppt_regression_points(i),i=1,3)
         read(SHORT_USE_FILE,*) (Parms(iunit)%temperature_production(i),i=1,4)
         write(ECHO_FILE,*) 'Temperature regression points, of four: ',  (Parms(iunit)%temperature_production(i),i=1,4)
         read(SHORT_USE_FILE,*) Parms(iunit)%standing_dead_production_halved
         write(ECHO_FILE,*) 'Standing dead production halved: ',  Parms(iunit)%standing_dead_production_halved
         read(SHORT_USE_FILE,*) Parms(iunit)%radiation_production_coefficient
         write(ECHO_FILE,*) 'Radiation production coefficient: ',  Parms(iunit)%radiation_production_coefficient
         read(SHORT_USE_FILE,*) (Parms(iunit)%fraction_carbon_to_roots(i),i=1,3)
         write(ECHO_FILE,*) 'Fraction carbon going to roots: ',  (Parms(iunit)%fraction_carbon_to_roots(i),i=1,3)
         read(SHORT_USE_FILE,*) Parms(iunit)%grazing_effect
         write(ECHO_FILE,*) 'Grazing effect integer, 1 through 6: ',  Parms(iunit)%grazing_effect
         Parms(iunit)%grazing_effect_multiplier = 0.0   ! read(SHORT_USE_FILE,*) Parms(iunit)%grazing_effect_multiplier
         write(ECHO_FILE,*) 'Grazing effect multiplier: ',  Parms(iunit)%grazing_effect_multiplier
         Parms(iunit)%temperature_effect_decomposition(1) =  15.40  ! Century values
         Parms(iunit)%temperature_effect_decomposition(2) =  11.75
         Parms(iunit)%temperature_effect_decomposition(3) =  29.70
         Parms(iunit)%temperature_effect_decomposition(4) =  0.031
         write(ECHO_FILE,*) 'Temperature effect on decomposition, four: ',  (Parms(iunit)%temperature_effect_decomposition(i),i=1,4)
         Parms(iunit)%anerobic_effect_decomposition(1) = 1.5       ! read(SHORT_USE_FILE,*) (Parms(iunit)%anerobic_effect_decomposition(i),i=1,3)
         Parms(iunit)%anerobic_effect_decomposition(2) = 3.0
         Parms(iunit)%anerobic_effect_decomposition(3) = 0.5
         write(ECHO_FILE,*) 'Anerobic effect on decomposition, of three: ',  (Parms(iunit)%anerobic_effect_decomposition(i),i=1,3)
         Parms(iunit)%decomp_rate_structural_litter(1) = 3.9   !  read(SHORT_USE_FILE,*) (Parms(iunit)%decomp_rate_structural_litter(i),i=1,2)
         Parms(iunit)%decomp_rate_structural_litter(2) = 4.9
         write(ECHO_FILE,*) 'Decomposition rate of structural litter, surface and soil: ', &
                                                        (Parms(iunit)%decomp_rate_structural_litter(i),i=1,2)
         Parms(iunit)%decomp_rate_metabolic_litter(1) = 14.8    !read(SHORT_USE_FILE,*) (Parms(iunit)%decomp_rate_metabolic_litter(i),i=1,2)
         Parms(iunit)%decomp_rate_metabolic_litter(2) = 18.5
         write(ECHO_FILE,*) 'Decomposition rate of metabolic litter, surface and soil: ',  &
                                                        (Parms(iunit)%decomp_rate_metabolic_litter(i),i=1,2)
         Parms(iunit)%decomp_rate_fast_som(1) = 6.0     ! read(SHORT_USE_FILE,*) (Parms(iunit)%decomp_rate_fast_som(i),i=1,2)
         Parms(iunit)%decomp_rate_fast_som(2) = 7.3
         write(ECHO_FILE,*) 'Decomposition rate of fast soil organic matter, surface and soil: ', &
                                                        (Parms(iunit)%decomp_rate_fast_som(i),i=1,2)
         Parms(iunit)%decomp_rate_slow_som = 0.00450    ! read(SHORT_USE_FILE,*) Parms(iunit)%decomp_rate_slow_som
         write(ECHO_FILE,*) 'Decomposition rate of slow soil organic matter: ',  Parms(iunit)%decomp_rate_slow_som
         Parms(iunit)%decomp_rate_inter_som = 0.2       ! read(SHORT_USE_FILE,*) Parms(iunit)%decomp_rate_inter_som
         write(ECHO_FILE,*) 'Decomposition rate of intermediate soil organic matter: ',  Parms(iunit)%decomp_rate_inter_som
         Parms(iunit)%decomp_rate_fine_branch = 1.5     ! read(SHORT_USE_FILE,*) Parms(iunit)%decomp_rate_fine_branch
         write(ECHO_FILE,*) 'Decomposition rate of fine branches: ',  Parms(iunit)%decomp_rate_fine_branch
         Parms(iunit)%decomp_rate_coarse_branch = 0.5   ! read(SHORT_USE_FILE,*) Parms(iunit)%decomp_rate_coarse_branch
         write(ECHO_FILE,*) 'Decomposition rate of coarse branches: ',  Parms(iunit)%decomp_rate_coarse_branch
         Parms(iunit)%decomp_rate_coarse_root = 0.6     ! read(SHORT_USE_FILE,*) Parms(iunit)%decomp_rate_coarse_root
         write(ECHO_FILE,*) 'Decomposition rate of coarse roots: ',  Parms(iunit)%decomp_rate_coarse_root
         read(SHORT_USE_FILE,*) (Parms(iunit)%decomp_rate_structural_litter_inverts(i),i=1,2)
         write(ECHO_FILE,*) 'Decomposition rate of structural litter by invertebrates, surface and soil: ',  &
                                                        (Parms(iunit)%decomp_rate_structural_litter_inverts(i),i=1,2)
         read(SHORT_USE_FILE,*) Parms(iunit)%drainage_affecting_anaerobic_decomp
         write(ECHO_FILE,*) 'Drainage affecting anaerobic decomposition: ', Parms(iunit)%drainage_affecting_anaerobic_decomp
         read(SHORT_USE_FILE,*) Parms(iunit)%feces_lignin
         write(ECHO_FILE,*) 'Lignin concentration in feces: ',  Parms(iunit)%feces_lignin
         read(SHORT_USE_FILE,*) ((Parms(iunit)%lignin_content_fraction_and_precip(i,j),j=1,2),i=1,2)
         write(ECHO_FILE,*) 'Lignin content fraction and precipitation, two by two: ', &
                                                       ((Parms(iunit)%lignin_content_fraction_and_precip(i,j),j=1,2),i=1,2)
         read(SHORT_USE_FILE,*) Parms(iunit)%fraction_urine_volatized
         write(ECHO_FILE,*) 'Fraction of urine volatized: ',  Parms(iunit)%fraction_urine_volatized
         Parms(iunit)%fraction_gross_n_mineral_volatized = 0.05  ! read(SHORT_USE_FILE,*) Parms(iunit)%fraction_gross_n_mineral_volatized
         write(ECHO_FILE,*) 'Fraction gross nitrogen mineral volatized: ',  Parms(iunit)%fraction_gross_n_mineral_volatized
         Parms(iunit)%rate_volatization_mineral_n = 0.02   ! read(SHORT_USE_FILE,*) Parms(iunit)%rate_volatization_mineral_n
         write(ECHO_FILE,*) 'Rate of volatization of mineral nitrogen: ',  Parms(iunit)%rate_volatization_mineral_n
         read(SHORT_USE_FILE,*) (Parms(iunit)%precip_n_deposition(i),i=1,2)
         write(ECHO_FILE,*) 'Precipitation nitrogen from deposition, of two: ',  (Parms(iunit)%precip_n_deposition(i),i=1,2)
!         read(SHORT_USE_FILE,*) (Parms(iunit)%precip_n_symbiotic(i),i=1,2)
!         write(ECHO_FILE,*) 'Precipitation nitrogen from symbiotics, of two: ',  (Parms(iunit)%precip_n_symbiotic(i),i=1,2)
         read(SHORT_USE_FILE,*) Parms(iunit)%decomp_litter_mix_facets
         write(ECHO_FILE,*) 'Litter mix between facets for decomposition: ',  Parms(iunit)%decomp_litter_mix_facets
         ! New approach to read in degree_days_phen.  No change in the science, just the formatting, to leave the parameter file more structured
         ! The echoed portion must match with past simulations, to ensure it is being read in correctly.  After some struggles they do in fact agree.
         do i=1,FACETS
           Parms(iunit)%degree_days_phen(i,1)  = 0.0
           Parms(iunit)%degree_days_phen(i,2)  = 0.0
           Parms(iunit)%degree_days_phen(i,4)  = 1.0  ! 3 controls phen stage 1
           Parms(iunit)%degree_days_phen(i,6)  = 2.0  ! 5 controls phen stage 2
           Parms(iunit)%degree_days_phen(i,8)  = 3.0  ! 7 controls phen stage 3
           Parms(iunit)%degree_days_phen(i,10) = 4.0  ! 9 controls phen stage 4
         end do
         read(SHORT_USE_FILE,*) (Parms(iunit)%degree_days_phen(H_FACET,j),j=3,10,2),&
                                (Parms(iunit)%degree_days_phen(S_FACET,j),j=3,10,2),&
                                (Parms(iunit)%degree_days_phen(T_FACET,j),j=3,10,2)
!         read(SHORT_USE_FILE,*) ((Parms(iunit)%degree_days_phen(i,j),j=1,10),i=1,FACETS)
         write(ECHO_FILE,*) 'Degree days and phenology, per facet: ',  ((Parms(iunit)%degree_days_phen(i,j),j=1,10),i=1,FACETS)
         read(SHORT_USE_FILE,*) (Parms(iunit)%degree_days_reset(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Degree days to reset phenology, per facet: ',  (Parms(iunit)%degree_days_reset(i),i=1,FACETS)
         Parms(iunit)%root_effect_on_nutrients = 0.015    ! read(SHORT_USE_FILE,*) Parms(iunit)%root_effect_on_nutrients
         write(ECHO_FILE,*) 'Root effect on nutrients: ',  Parms(iunit)%root_effect_on_nutrients
         Parms(iunit)%root_intercept_on_nutrients = 0.8   ! read(SHORT_USE_FILE,*) Parms(iunit)%root_intercept_on_nutrients
         write(ECHO_FILE,*) 'Root intercept on nutrients: ',  Parms(iunit)%root_intercept_on_nutrients
         read(SHORT_USE_FILE,*) Parms(iunit)%tree_site_potential
         write(ECHO_FILE,*) 'Tree site potential: ',  Parms(iunit)%tree_site_potential
         Parms(iunit)%tree_basal_area_to_grass_nitrogen = 1.0   ! read(SHORT_USE_FILE,*) Parms(iunit)%tree_basal_area_to_grass_nitrogen
         write(ECHO_FILE,*) 'Tree basal area versus grass nitrogen: ',  Parms(iunit)%tree_basal_area_to_grass_nitrogen
         Parms(iunit)%tree_basal_area_to_wood_biomass = 400.0    ! read(SHORT_USE_FILE,*) Parms(iunit)%tree_basal_area_to_wood_biomass
         write(ECHO_FILE,*) 'Tree basal area versus woody biomass: ',  Parms(iunit)%tree_basal_area_to_wood_biomass
         read(SHORT_USE_FILE,*) Parms(iunit)%max_symbiotic_n_fixation_ratio
         write(ECHO_FILE,*) 'Maximum symbiotic nitrogen fixation ratio: ',  Parms(iunit)%max_symbiotic_n_fixation_ratio
         Parms(iunit)%fraction_nitrogen_available = 0.9   ! read(SHORT_USE_FILE,*) Parms(iunit)%fraction_nitrogen_available
         write(ECHO_FILE,*) 'Fraction of nitrogen available: ',  Parms(iunit)%fraction_nitrogen_available
         read(SHORT_USE_FILE,*) ((Parms(iunit)%minimum_c_n_ratio(i,j),j=1,WOODY_PARTS),i=1,FACETS)
         write(ECHO_FILE,*) 'Minimum carbon to nitrogen ratio, by facet, by woody part: ', &
                                                         ((Parms(iunit)%minimum_c_n_ratio(i,j),j=1,WOODY_PARTS),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%maximum_c_n_ratio(i,j),j=1,WOODY_PARTS),i=1,FACETS)
         write(ECHO_FILE,*) 'Maximum carbon to nitrogen ratio, by facet, by woody part: ', &
                                                         ((Parms(iunit)%maximum_c_n_ratio(i,j),j=1,WOODY_PARTS),i=1,FACETS)
         do i = 1, FACETS
           Parms(iunit)%fraction_npp_to_respiration(i) = 1.0  !  read(SHORT_USE_FILE,*) (Parms(iunit)%fraction_npp_to_respiration(i),i=1,FACETS)
         end do
         write(ECHO_FILE,*) 'Fraction of net primary production going to respiration, by facet: ', &
                                                         (Parms(iunit)%fraction_npp_to_respiration(i),i=1,FACETS)
         Parms(iunit)%herb_max_fraction_npp_to_respiration(ABOVE) = 0.26   !  read(SHORT_USE_FILE,*) (Parms(iunit)%herb_max_fraction_npp_to_respiration(i),i=1,2)
         Parms(iunit)%herb_max_fraction_npp_to_respiration(BELOW) = 0.26
         write(ECHO_FILE,*) 'Herb fraction of net primary production going to respiration, above and below: ', &
                                                         (Parms(iunit)%herb_max_fraction_npp_to_respiration(i),i=1,2)
         Parms(iunit)%woody_max_fraction_npp_to_respiration(1) = 0.4   ! read(SHORT_USE_FILE,*) (Parms(iunit)%woody_max_fraction_npp_to_respiration(i),i=1,5)
         Parms(iunit)%woody_max_fraction_npp_to_respiration(2) = 0.4
         Parms(iunit)%woody_max_fraction_npp_to_respiration(3) = 0.4
         Parms(iunit)%woody_max_fraction_npp_to_respiration(4) = 0.4
         Parms(iunit)%woody_max_fraction_npp_to_respiration(5) = 0.4
         write(ECHO_FILE,*) 'Woody fraction of net primary production going to respiration, by woody part: ', &
                                                         (Parms(iunit)%woody_max_fraction_npp_to_respiration(i),i=1,5)
         read(SHORT_USE_FILE,*) Parms(iunit)%maximum_leaf_area_index
         write(ECHO_FILE,*) 'Maximum leaf area index: ',  Parms(iunit)%maximum_leaf_area_index
         read(SHORT_USE_FILE,*) Parms(iunit)%k_leaf_area_index
         write(ECHO_FILE,*) 'Biomass where leaf area index is half of maximum: ',  Parms(iunit)%k_leaf_area_index
         read(SHORT_USE_FILE,*) Parms(iunit)%biomass_to_leaf_area_index_factor
         write(ECHO_FILE,*) 'Biomass to leaf area index factor: ',  Parms(iunit)%biomass_to_leaf_area_index_factor
         read(SHORT_USE_FILE,*) Parms(iunit)%annual_fraction_volatilized_n
         write(ECHO_FILE,*) 'Annual fraction volatized nitrogen: ',  Parms(iunit)%annual_fraction_volatilized_n
         read(SHORT_USE_FILE,*) Parms(iunit)%max_herb_root_death_rate
         write(ECHO_FILE,*) 'Maximum herb root death rate, by facet: ',  Parms(iunit)%max_herb_root_death_rate
         Parms(iunit)%fraction_n_absorbed_by_residue(1) = 0.0000    ! read(SHORT_USE_FILE,*) (Parms(iunit)%fraction_n_absorbed_by_residue(i),i=1,2)
         Parms(iunit)%fraction_n_absorbed_by_residue(2) = 0.0200
         write(ECHO_FILE,*) 'Fraction nitrogen absorbed by residue, two: ',(Parms(iunit)%fraction_n_absorbed_by_residue(i),i=1,2)
         read(SHORT_USE_FILE,*) (Parms(iunit)%shoot_death_rate(i),i=1,4)
         write(ECHO_FILE,*) 'Shoot death rate, in four catagories: ',  (Parms(iunit)%shoot_death_rate(i),i=1,4)
         read(SHORT_USE_FILE,*) Parms(iunit)%prop_annuals
         write(ECHO_FILE,*) 'Proportion annuals in herb facet: ',  Parms(iunit)%prop_annuals
         read(SHORT_USE_FILE,*) Parms(iunit)%month_to_remove_annuals
         write(ECHO_FILE,*) 'Month to remove remaining annuals: ',  Parms(iunit)%month_to_remove_annuals
         read(SHORT_USE_FILE,*) (Parms(iunit)%relative_seed_production(i),i=1,FACETS)
         do ifacet = 1, FACETS
           Parms(iunit)%relative_seed_production(ifacet) = Parms(iunit)%relative_seed_production(ifacet) / 10000.                  ! To avoid overflows.  Perhaps include a check of the summed values.
         end do
         write(ECHO_FILE,*) 'Relative seed production, by facet: ',  (Parms(iunit)%relative_seed_production(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%water_effect_on_establish(i,j),j=1,4),i=1,FACETS)
         ! The following may need to be made a variable.  Values are HARDWIRED but not well known.
         Parms(iunit)%fraction_aground_npp_to_seeds(H_FACET) = 0.05
         Parms(iunit)%fraction_aground_npp_to_seeds(S_FACET) = 0.03
         Parms(iunit)%fraction_aground_npp_to_seeds(T_FACET) = 0.01 ! read(SHORT_USE_FILE,*) ((Parms(iunit)%fraction_aground_npp_to_seeds(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Fraction of aboveground npp to go to seeds, by facet: ', &
                                                     (Parms(iunit)%fraction_aground_npp_to_seeds(i),i=1,FACETS)
         ! The following may need to be made a variable.  Values are HARDWIRED but not well known.
         Parms(iunit)%fraction_seeds_not_germinated(H_FACET) = 0.50
         Parms(iunit)%fraction_seeds_not_germinated(S_FACET) = 0.80
         Parms(iunit)%fraction_seeds_not_germinated(T_FACET) = 0.90 ! read(SHORT_USE_FILE,*) ((Parms(iunit)%fraction_seeds_not_germinated(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Fraction of seeds that are not germinated, by facet: ', &
                                                     (Parms(iunit)%fraction_seeds_not_germinated(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Water effect on establishment, by facet, by four: ', &
                                                          ((Parms(iunit)%water_effect_on_establish(i,j),j=1,4),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%herb_root_effect_on_establish(i,j),j=1,4),i=1,FACETS)
         write(ECHO_FILE,*) 'Herbaceous root effect on establishment, by facet, by four: ', &
                                                         ((Parms(iunit)%herb_root_effect_on_establish(i,j),j=1,4),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%litter_effect_on_establish(i,j),j=1,4),i=1,FACETS)
         write(ECHO_FILE,*) 'Litter effect on establishment, by facet, by four: ', &
                                                         ((Parms(iunit)%litter_effect_on_establish(i,j),j=1,4),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%woody_cover_effect_on_establish(i,j),j=1,4),i=1,FACETS)
         write(ECHO_FILE,*) 'Woody cover effect on establishment, by facet, by four: ', &
                                                         ((Parms(iunit)%woody_cover_effect_on_establish(i,j),j=1,4),i=1,FACETS)
         read(SHORT_USE_FILE,*) (Parms(iunit)%nominal_plant_death_rate(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Nominal plant death rate, by facet: ',  (Parms(iunit)%nominal_plant_death_rate(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%water_effect_on_death_rate(i,j),j=1,4),i=1,FACETS)
         write(ECHO_FILE,*) 'Water effect on death rate, by facet, by four: ', &
                                                         ((Parms(iunit)%water_effect_on_death_rate(i,j),j=1,4),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%grazing_effect_on_death_rate(i,j),j=1,4),i=1,FACETS)
         write(ECHO_FILE,*) 'Grazing effect on death rate, by facet, by four: ', &
                                                         ((Parms(iunit)%grazing_effect_on_death_rate(i,j),j=1,4),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%shading_effect_on_death_rate(i,j),j=1,4),i=1,FACETS)
         write(ECHO_FILE,*) 'Shading effect on death rate, by facet, by four: ', &
                                                         ((Parms(iunit)%shading_effect_on_death_rate(i,j),j=1,4),i=1,FACETS)
         read(SHORT_USE_FILE,*) (Parms(iunit)%fall_rate_of_standing_dead(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Fall rate of standing dead to litter, by facet: ', &
                                                         (Parms(iunit)%fall_rate_of_standing_dead(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) Parms(iunit)%death_rate_of_deciduous_leaves
         write(ECHO_FILE,*) 'Death rate of deciduous leaves: ',  Parms(iunit)%death_rate_of_deciduous_leaves
         Parms(iunit)%temperature_leaf_out_and_fall(1) = 10.0   ! read(SHORT_USE_FILE,*) (Parms(iunit)%temperature_leaf_out_and_fall(i),i=1,2)
         Parms(iunit)%temperature_leaf_out_and_fall(2) = 7.0
         write(ECHO_FILE,*) 'Temperature leaf-out and leaf-fall: ',  (Parms(iunit)%temperature_leaf_out_and_fall(i),i=1,2)
         read(SHORT_USE_FILE,*) (Parms(iunit)%drought_deciduous(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Proportion that are drought deciduous, by facet: ',  (Parms(iunit)%drought_deciduous(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) Parms(iunit)%fraction_woody_leaf_n_translocated
         write(ECHO_FILE,*) 'Fraction of woody leaf nitrogen translocated prior to scenecense: ', &
                                                           Parms(iunit)%fraction_woody_leaf_n_translocated
         read(SHORT_USE_FILE,*) (Parms(iunit)%leaf_death_rate(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Leaf death rate, by facet: ',  (Parms(iunit)%leaf_death_rate(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) (Parms(iunit)%fine_root_death_rate(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Fine root death rate, by facet: ',  (Parms(iunit)%fine_root_death_rate(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) (Parms(iunit)%fine_branch_death_rate(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Fine branch death rate, by facet: ',  (Parms(iunit)%fine_branch_death_rate(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) (Parms(iunit)%coarse_branch_death_rate(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Coarse branch death rate, by facet: ',  (Parms(iunit)%coarse_branch_death_rate(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) (Parms(iunit)%coarse_root_death_rate(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Coarse root death rate, by facet: ',  (Parms(iunit)%coarse_root_death_rate(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) Parms(iunit)%fraction_carbon_grazed_returned
         write(ECHO_FILE,*) 'Fraction of carbon grazed returned through feces: ',  Parms(iunit)%fraction_carbon_grazed_returned
         read(SHORT_USE_FILE,*) Parms(iunit)%fraction_excreted_nitrogen_in_feces
         write(ECHO_FILE,*) 'Fraction excreted nitrogen that is in feces: ',  Parms(iunit)%fraction_excreted_nitrogen_in_feces
         read(SHORT_USE_FILE,*) (Parms(iunit)%fraction_grazed_by_facet(i),i=1,FACETS)
         write(ECHO_FILE,*) 'Fraction grazed, by facet: ',  (Parms(iunit)%fraction_grazed_by_facet(i),i=1,FACETS)
         read(SHORT_USE_FILE,*) Parms(iunit)%fraction_grazed
         write(ECHO_FILE,*) 'Fraction grazed: ',  Parms(iunit)%fraction_grazed
         read(SHORT_USE_FILE,*) Parms(iunit)%frequency_of_fire
         write(ECHO_FILE,*) 'Frequency of fire: ', Parms(iunit)%frequency_of_fire
         read(SHORT_USE_FILE,*) Parms(iunit)%fraction_burned
         write(ECHO_FILE,*) 'Fraction burned: ', Parms(iunit)%fraction_burned
         read(SHORT_USE_FILE,*) Parms(iunit)%burn_month
         write(ECHO_FILE,*) 'Burn month: ', Parms(iunit)%burn_month
         read(SHORT_USE_FILE,*) (Parms(iunit)%fuel_vs_intensity(i),i=1,FIRE_SEVERITIES)
         write(ECHO_FILE,*) 'Fuel versus intensity: ', (Parms(iunit)%fuel_vs_intensity(i),i=1,FIRE_SEVERITIES)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%green_vs_intensity(i,j),j=1,FIRE_SEVERITIES),i=1,2)
         write(ECHO_FILE,*) 'Green versus intensity: ', ((Parms(iunit)%green_vs_intensity(i,j),j=1,FIRE_SEVERITIES),i=1,2)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%fraction_shoots_burned(i,j),j=1,FIRE_SEVERITIES),i=1,FACETS)
         write(ECHO_FILE,*) 'Fraction shoots burned: ', ((Parms(iunit)%fraction_shoots_burned(i,j),j=1,FIRE_SEVERITIES),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%fraction_standing_dead_burned(i,j),j=1,FIRE_SEVERITIES),i=1,FACETS)
         write(ECHO_FILE,*) 'Fraction standing dead burned: ', &
                            ((Parms(iunit)%fraction_standing_dead_burned(i,j),j=1,FIRE_SEVERITIES),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%fraction_plants_burned_dead(i,j),j=1,FIRE_SEVERITIES),i=1,FACETS)
         write(ECHO_FILE,*) 'Fraction plants burned dead: ', &
                            ((Parms(iunit)%fraction_plants_burned_dead(i,j),j=1,FIRE_SEVERITIES),i=1,FACETS)
         read(SHORT_USE_FILE,*) ((Parms(iunit)%fraction_litter_burned(i,j),j=1,FIRE_SEVERITIES),i=1,FACETS)
         write(ECHO_FILE,*) 'Fraction litter burned: ', ((Parms(iunit)%fraction_litter_burned(i,j),j=1,FIRE_SEVERITIES),i=1,FACETS)
         read(SHORT_USE_FILE,*) Parms(iunit)%fraction_burned_carbon_as_ash
         write(ECHO_FILE,*) 'Fraction burned carbon that is ash: ', Parms(iunit)%fraction_burned_carbon_as_ash
         read(SHORT_USE_FILE,*) Parms(iunit)%fraction_burned_nitrogen_as_ash
         write(ECHO_FILE,*) 'Fraction burned nitrogen that is ash: ', Parms(iunit)%fraction_burned_nitrogen_as_ash
         read(SHORT_USE_FILE,*) Parms(iunit)%frequency_of_fertilization
         write(ECHO_FILE,*) 'Frequency of fertilization: ', Parms(iunit)%frequency_of_fertilization
         read(SHORT_USE_FILE,*) Parms(iunit)%fraction_fertilized
         write(ECHO_FILE,*) 'Fraction fertilized: ', Parms(iunit)%fraction_fertilized
         read(SHORT_USE_FILE,*) Parms(iunit)%fertilize_month
         write(ECHO_FILE,*) 'Fertilize month: ', Parms(iunit)%fertilize_month
         read(SHORT_USE_FILE,*) Parms(iunit)%fertilize_nitrogen_added
         write(ECHO_FILE,*) 'Fertilize nitrogen added: ', Parms(iunit)%fertilize_nitrogen_added
         read(SHORT_USE_FILE,*) Parms(iunit)%fertilize_carbon_added
         write(ECHO_FILE,*) 'Fertilize carbon added: ', Parms(iunit)%fertilize_carbon_added
!         read(SHORT_USE_FILE,*) Parms(iunit)%fertilize_carbon_nitrogen_ratio
!         write(ECHO_FILE,*) 'Fertilize carbon nitrogen ratio: ', Parms(iunit)%fertilize_carbon_nitrogen_ratio
       end if
     end do
     write(*,*) 'Landscape unit parameters have been read in.'
   else
     write(ECHO_FILE,*) 'The landscape unit parameter file cannot be opened: ', &
                        app_path(1:len_trim(app_path))//'/Parms/'//Sim_Parm%parms_file_name
   end if
   close(SHORT_USE_FILE)

end subroutine


subroutine Set_Defaults
!*****  Set initial values to 0 or other values as a first step in the simulation.
!*****  This will need reviewed.
!*****
!*****  R Boone    Last Modified:  Oct 2, 2010
!*****
  use Parameter_Vars
  use Structures
  implicit none

  integer icell, iunit, ilayer, ifacet

  do icell=1,MAX_RANGE_CELLS

    iunit = Rng(icell)%range_type
    Rng(icell)%snow = 0.0
    Rng(icell)%old_snow = 0.0
    Rng(icell)%melt = 0.0
    Rng(icell)%snow_liquid = 0.0
    Rng(icell)%old_snow_liquid = 0.0
    Rng(icell)%pet_remaining = 0.0
    Rng(icell)%ppt_soil = 0.0
    Rng(icell)%evaporation = 0.0
    Rng(icell)%water_available(1) = 0.0
    Rng(icell)%water_available(2) = 0.0
    Rng(icell)%water_available(3) = 0.0
    do ifacet = 1, FACETS
      Rng(icell)%facet_cover(ifacet) = 0.0
    end do
    do ilayer = 1, V_LYRS
      Rng(icell)%total_population(ilayer) = 0.0
    end do
    Rng(icell)%bare_cover = 0.0
    do ilayer = 1, SOIL_LAYERS
      Rng(icell)%mineral_nitrogen(ilayer) = 0.0
    end do
    Rng(icell)%large_error_count = 0
    Rng(icell)%neg_error_count = 0
  end do


end subroutine

