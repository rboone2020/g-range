program GRange
  ! G-Range is a global model that simulates plant and animal dynamics on rangelands.
  ! The following uses calls per rangeland cell.  That may be less than efficient, but
  ! most easily parallizable.  Details may be seen in supplemental material from Boone et al. (Global Change Biology, 2018).
  !
  ! R.B. Boone and R.T. Conant  Last modified: April 9, 2019. Last change of logic: July 1, 2014
  use Structures
  use Parameter_Vars
  implicit none

  integer icell

  call Initialize_Parms                          ! Read the initial parameters for the simulation.
  call Initialize_Landscape_Parms                ! Read the parameters that correspond to the different landscape units simulated.
  call Initialize_Globe                          ! Do the spatial initialization, setting up the global data.
  call Initialize_Rangelands                     ! Do additional initialization of rangelands, including checking state variable file.
  call Initialize_Outputs                        ! Produces positional data for rangeland cells in the output file.

  do year = Sim_Parm%start_yr, Sim_Parm%end_Yr   ! Loop through each year ...
    call Each_Year                               ! Miscellaneous steps that need to be done each year
    do month = 1, 12                             ! ... and each month.  Simulations must end in December
      call Progress                              ! Report the progress of the simulation
      call Read_Weather                          ! Get the month's weather data
      call Read_Other                            ! Get any maps associated with fire or management

      do icell=1,range_cells                     ! Process all of the cells classed as rangeland only
        call Update_Vegetation (icell)           ! Update metrics for vegetation.
        call Update_Weather (icell)              ! Calculate snowfall, evapotranspiration, etc.  Also updates heat accumulation.
        call Potential_Production (icell)        ! Calculate potential production, and plant allometrics adjusted by grazing fraction
        call Herb_Growth (icell)                 ! Calculate herbaceous growth
        call Woody_Growth (icell)                ! Calculate woody plant growth
        call Grazing (icell)                     ! Remove material grazed by livestock
        call Plant_Part_Death (icell)            ! Plant part death
        call Whole_Plant_Death (icell)           ! Whole plant death
        call Management (icell)                  ! Fertilization and other management
        call Plant_Reproduction (icell)          ! Seed-based reproduction by plants
        call Update_Vegetation (icell)           ! Update metrics for vegetation
        call Water_Loss (icell)                  ! Calculate water loss
        call Decomposition (icell)               ! Decomposition
        call Nitrogen_Losses (icell)             ! Leaching and volatilization of nitrogen
      end do

      call Each_Month                            ! Miscellaneous steps that need to be done each month
      call Output_Surfaces                       ! Produce output surfaces
      call Zero_Accumulators                     ! Zero-out the accumulators storing dead materials

    end do                                       ! End monthly loop
  end do                                         ! End yearly loop

  call Wrap_Up                                   ! Do simulation-ending tasks, including some output
  write(*,*) ' '
  write(*,*) 'G-Range is complete.'

end program
