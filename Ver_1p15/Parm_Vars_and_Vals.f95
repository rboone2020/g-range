module Parameter_Vars
  ! The Parameter_Vars module stores parameters used in the G-Range model, and constants that are needed.  
   
  ! Constants
  integer, parameter  :: MAX_X_DIM=1440, MAX_Y_DIM=720              ! MAX_X_DIM=720, MAX_Y_DIM=360               ! MAX_X_DIM=3600, MAX_Y_DIM=1800             ! Appropraite for a 0.1 degree grid of the globe.  World arrays should be made dynamic at some point, instead of hardwired. 
  integer, parameter  :: MAX_RANGE_CELLS=80000                      ! MAX_RANGE_CELLS=60000   ! MAX_RANGE_CELLS=1200000
  character(11)       :: INI_FILE='\GRange.ini'                     ! The one-line init file in the program directory, pointing to the application directory.
  character(19)       :: SIMPARM_NAME='\Parms\Sim_Parm.grg'         ! The main parameter file
  character(22)       :: OUT_LIST_NAME='\Parms\Output_List.grg'     ! The list of output choices, showing what the user wishes to include in output
  character(16)       :: ECHO_NAME='\Output\Echo.gof'               ! The file storing echoed input or other program data
  character(20)       :: OUT_DATA_NAME='\Output\Rng_Data.gof'       ! The list of output choices, showing what the user wishes to include in output
  character(17)       :: OUT_MAP_NAME= '\Output\Globe.gof'          ! The globe, storing 0s and 1s showing land and ocean
  character(20)       :: RUN_TIME_NAME='\Output\Run_Time.gof'       ! The file that stores simulation time
  character(18)       :: EXCEED_NAME='\Output\Exceed.gof'           ! The name of the file storing the number of times cell contents were reset because they exceeded some limit (i.e., V_LARGE or 0.0)
  integer, parameter  :: ECHO_FILE=99                               ! The unit number for the echo file
  integer, parameter  :: SHORT_USE_FILE=98                          ! The unit number for short-term file use
  ! Vegetation layer constants
  integer, parameter  :: H_LYR=1                                    ! Herb layer index
  integer, parameter  :: H_S_LYR=2                                  ! Herbs under shrubs layer index
  integer, parameter  :: H_T_LYR=3                                  ! Herbs under tree layer index
  integer, parameter  :: S_LYR=4                                    ! Shrub layer index
  integer, parameter  :: S_T_LYR=5                                  ! Shrub under tree layer index
  integer, parameter  :: T_LYR=6                                    ! Tree layer index
  integer, parameter  :: V_LYRS=6                                   ! Total number of layers (v used to distinguish from other layers)
  ! Vegetation facet constants.  Facets are used for unit input, so that users provide 3 values, rather than 6
  integer, parameter  :: H_FACET=1                                  ! Herb layer index
  integer, parameter  :: S_FACET=2                                  ! Shrub facet index
  integer, parameter  :: T_FACET=3                                  ! Herbs under tree layer index  
  integer, parameter  :: FACETS=3
  ! Woody part constants
  integer, parameter  :: WOODY_PARTS=5                             ! Woody parts, leaf, fine root, fine branch, large wood, and coarse root    
  integer, parameter  :: LEAF_INDEX=1
  integer, parameter  :: FINE_ROOT_INDEX=2
  integer, parameter  :: FINE_BRANCH_INDEX=3
  integer, parameter  :: COARSE_BRANCH_INDEX=4
  integer, parameter  :: COARSE_ROOT_INDEX=5                        ! Indices for woody parts.
  ! Soil and layer constants
  integer, parameter  :: SOIL_LAYERS=4                              ! The number of soil layers
  integer, parameter  :: SURFACE_INDEX=1                            ! Surface index in litter array and perhaps elsewhere
  integer, parameter  :: SOIL_INDEX=2                               ! Soil index in litter array and perhaps elsewhere                                                                                                   
  integer, parameter  :: N_DECOMP_LOOPS=4                           ! Number of times per month decomposition is called, defaults to 4 and a warning not to change it in CENTURY docs.
  integer, parameter  :: N_STORE=1                                  ! Element to storage, using in Growth, etc.
  integer, parameter  :: N_SOIL=2                                   ! Element to soil, used in Growth, etc.
  integer, parameter  :: N_FIX=3                                    ! Element to fix, used in Growth, etc.
  integer, parameter  :: ABOVE=1                                    ! Aboveground index
  integer, parameter  :: BELOW=2                                    ! Belowground index
  integer, parameter  :: ALIVE=1                                    ! Plant material that is alive, index
  integer, parameter  :: DEAD=2                                     ! Plant material that is dead, index
  ! Other constants
  real, parameter     :: PI=3.141592653589793                       ! Pi ...
  real, parameter     :: PI2=6.283185307179586                      ! Pi doubled
  real, parameter     :: WRAD=0.0174532925                          ! Radians
  real, parameter     :: BASE_TEMP=4.4                              ! Base temperature for heat accumulation
  integer, parameter  :: JULIAN_DAY_MID(12)=(/ 16,46,75,106,136,167,197,228,259,289,320,350 /)  ! The Julian day at the middle of each month.
  integer, parameter  :: JULIAN_DAY_START(12)=(/ 1,32,61,92,122,153,183,214,245,275,306,337 /)  ! The Julian day at the start of each month.
  integer, parameter  :: MONTH_DAYS(12)=(/31,28,31,30,31,30,31,31,30,31,30,31 /)                !The days in each month
  integer, parameter  :: REF_AREA=1000000.                          ! A 1 km x 1 km square area, in meters.
  integer, parameter  :: WHOLE_MAP_COUNT=127                        ! Number of potential output variables ... UPDATE AS NEEDED
  real, parameter     :: V_LARGE=1000000000.                        ! A very large number; entries will not be allowed to exceed this value
  integer, parameter  :: FIRE_SEVERITIES=2                          ! The number of fire severities used
  ! Pathways
  character(100)      :: bin_path                                   ! The path to the G-RANGE program, filled using a Unix library call
  character(100)      :: app_path                                   ! The path to the application of interest
  character(100)      :: parm_path                                  ! The path to the parameter directory
  character(100)      :: out_path                                   ! The path to the output directory
  ! Basic globals
  integer             :: ioerr                                      ! Receives errors when opening files
  integer             :: x_dim, y_dim                               ! Actual maximum X and Y dimensions of the globe 
  integer             :: lower_x, lower_y, upper_x, upper_y         ! Lower and upper dimensions of the world
  real                :: cellsize                                   ! The resolution of cells, in decimal degrees            integer             :: 
  integer             :: range_cells                                ! The number of cells identified as rangeland
  integer             :: year                                       ! Year of simulation
  integer             :: month                                      ! Month of simulation
  real                :: start_time                                 ! The starting time of the simulation
  integer             :: clock_start_time                           ! The start time in clock units, rather than CPU units
  ! Other global variables
  logical             :: rangeland_classes(255)                     ! Rangeland classes corresponding to the land cover classes.  Just dimensioned quite large - many classes.
  ! Model control variables
  logical             :: check_nan_flag                             ! A flag indicating whether -NaN should be checked for.   This will slow the model, but help in debugging.
  logical             :: stop_on_nan_flag                           ! A flag indicating whether or not -NaN should stop the model.  If the model does not stop, error messages will scroll and the NAN will remain.
  
end module