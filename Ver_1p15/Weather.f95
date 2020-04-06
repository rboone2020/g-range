subroutine Read_Weather 
!**** Read in the PRECIPITATION, MAXIMUM TEMPERATURE, and MINIMUM TEMPERATURE maps.
!**** !! Recall to correct for the multipliers used in the GRIDASCII files !!
!**** 
!**** R. Boone    Last modified: January 16, 2014
   use Parameter_Vars
   use Structures
   implicit none
   integer ix, iy, w_year
   character(100) big_string
   character(6) mid_string
   logical file_exists 

   ! Create a ** DRAFT ** date string that is a substring of the file names.
   w_year = year
   select case (month)
     case (1) ; write(mid_string, fmt='(I4,A2)') w_year, '01'
     case (2) ; write(mid_string, fmt='(I4,A2)') w_year, '02'
     case (3) ; write(mid_string, fmt='(I4,A2)') w_year, '03'
     case (4) ; write(mid_string, fmt='(I4,A2)') w_year, '04'
     case (5) ; write(mid_string, fmt='(I4,A2)') w_year, '05'
     case (6) ; write(mid_string, fmt='(I4,A2)') w_year, '06'
     case (7) ; write(mid_string, fmt='(I4,A2)') w_year, '07'
     case (8) ; write(mid_string, fmt='(I4,A2)') w_year, '08'
     case (9) ; write(mid_string, fmt='(I4,A2)') w_year, '09'
     case (10) ; write(mid_string, fmt='(I4,A2)') w_year, '10'
     case (11) ; write(mid_string, fmt='(I4,A2)') w_year, '11'
     case (12) ; write(mid_string, fmt='(I4,A2)') w_year, '12'
   end select

   ! ********************************************
   ! Check to see if there is a precipitation map for that date.  We will assume the other two maps are available if the PPT map is.
   ! First, make the string name
   write(big_string,'(a85,a6,a9)') trim(Sim_Parm%precip_path_prefix),mid_string,adjustl(Sim_Parm%precip_temp_suffix)
   inquire(FILE=app_path(1:len_trim(app_path))//trim(adjustl(big_string)), EXIST=file_exists)   

   if (file_exists .eq. .FALSE.) then
     ! Now decrement w_year until the precipitation file is there.  Assumes there are 100 (or more) years of layers available.
     ! (I suppose the if then above and this do while could be combined, but I like that the do is not touched at all if the file is present)
     do while (file_exists .eq. .FALSE.)
       w_year = w_year - 100
       select case (month)
         case (1) ; write(mid_string, fmt='(I4,A2)') w_year, '01'
         case (2) ; write(mid_string, fmt='(I4,A2)') w_year, '02'
         case (3) ; write(mid_string, fmt='(I4,A2)') w_year, '03'
         case (4) ; write(mid_string, fmt='(I4,A2)') w_year, '04'
         case (5) ; write(mid_string, fmt='(I4,A2)') w_year, '05'
         case (6) ; write(mid_string, fmt='(I4,A2)') w_year, '06'
         case (7) ; write(mid_string, fmt='(I4,A2)') w_year, '07'
         case (8) ; write(mid_string, fmt='(I4,A2)') w_year, '08'
         case (9) ; write(mid_string, fmt='(I4,A2)') w_year, '09'
         case (10) ; write(mid_string, fmt='(I4,A2)') w_year, '10'
         case (11) ; write(mid_string, fmt='(I4,A2)') w_year, '11'
         case (12) ; write(mid_string, fmt='(I4,A2)') w_year, '12'
       end select

       ! ********************************************
       ! Check to see if there is a precipitation map for that date.  We will assume the other two maps are available.
       ! First, make the string name
       write(big_string,'(a85,a6,a9)') trim(Sim_Parm%precip_path_prefix),mid_string,adjustl(Sim_Parm%precip_temp_suffix)
       inquire(FILE=app_path(1:len_trim(app_path))//trim(adjustl(big_string)), EXIST=file_exists)   
     end do
   end if     ! End logical if.  If the file needed rewound by multiples of 100 years, the program has done so, so just proceed.
              ! If it wasn't required (i.e., if files were present), no action has occurred and the program may proceed normally.
   write(ECHO_FILE,*) 'To model precip in month ', month, ' of ', year, ' the program is using ', trim(adjustl(big_string))
   
   ! ********************************************
   ! Read in the PRECIPITATION map for the month.
   write(big_string,'(a85,a6,a9)') trim(Sim_Parm%precip_path_prefix),mid_string,adjustl(Sim_Parm%precip_temp_suffix)
   if (Sim_Parm%echo_level .gt. 0) write(*,*) 'Reading in the precipitation map and converting to cm: ',adjustl(big_string)
   call Read_Map ('precip     ', big_string)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary .lt. 0) then
         Globe(ix,iy)%precip = 0                                   ! Will assign null data 0 precip, but that is fine.
       else
         Globe(ix,iy)%precip = Globe(ix,iy)%temporary / 100.0      ! Stored as * 10 in the GRIDS, so / 100 takes that to cm.
       end if    
     end do
   end do
   ! These maps are not echoed, due to their large size.  Specific maps may be printed as output, rather than echoed, as a check of their form.

   ! **************************************************
   ! Read in the MAXIMUM TEMPERATURE map for the month.
   write(big_string,'(a85,a6,a9)') trim(Sim_Parm%max_temp_path_prefix),mid_string,adjustl(Sim_Parm%precip_temp_suffix)
   if (Sim_Parm%echo_level .gt. 0) write(*,*) 'Reading in the maximum temperature map and converting to degrees: ', &
                                               adjustl(big_string)
   call Read_Map ('max_temp      ', big_string)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary .lt. -400.0) then
         Globe(ix,iy)%max_temp = -40.
       else
         Globe(ix,iy)%max_temp = Globe(ix,iy)%temporary / 10.0
       end if    
     end do
   end do
   ! These maps are not echoed, due to their large size.  Specific maps may be printed as output, rather than echoed, as a check of their form.

   ! **************************************************
   ! Read in the MINIMUM TEMPERATURE map for the month.
   write(big_string,'(a85,a6,a9)') trim(Sim_Parm%min_temp_path_prefix),mid_string,adjustl(Sim_Parm%precip_temp_suffix)
   if (Sim_Parm%echo_level .gt. 0) write(*,*) 'Reading in the minimum temperature map and converting to degrees: ', &
                                               adjustl(big_string)
   call Read_Map ('min_temp      ', big_string)
   ! Transfer the map
   do iy=1,y_dim
     do ix=1,x_dim
       if (Globe(ix,iy)%temporary .lt. -400.0) then
         Globe(ix,iy)%min_temp = -40.                                             
       else
         Globe(ix,iy)%min_temp = Globe(ix,iy)%temporary / 10.0
       end if    
     end do
   end do
   ! These maps are not echoed, due to their large size.  Specific maps may be printed as output, rather than echoed, as a check of their form.

end subroutine



subroutine Update_Weather (icell)
!**** Main weather calling routine.
!**** 
!**** 
!**** R. Boone    Last modified: March 5, 2011
  use Parameter_Vars
  use Structures
  implicit none

  integer icell, iunit
  real temp_mean
  
  iunit = Rng(icell)%range_type
  Rng(icell)%storm_flow = Parms(iunit)%prcp_threshold_fraction * Globe(Rng(icell)%x, Rng(icell)%y)%precip
  temp_mean  = ( Globe(Rng(icell)%x,Rng(icell)%y)%max_temp + Globe(Rng(icell)%x,Rng(icell)%y)%min_temp ) / 2.0
  if (temp_mean .gt. BASE_TEMP) then
    Rng(icell)%heat_accumulation = Rng(icell)%heat_accumulation + ( ( temp_mean - BASE_TEMP ) * MONTH_DAYS(month) )
  end if

  call P_Evap (icell)
  call Day_Length (icell)
  call Snow_Dynamics (icell)

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'WEATHER')

end subroutine
   


subroutine Snow_Dynamics (icell)
!**** Calculate for all the rangeland cells any weather-related attributes, such as snowfall, snowmelt, 
!**** and evapotranspiration.
!**** This section draws heavily on the CENTURY code, as much of this material does.  Especially SNOWCENT
!**** R. Boone    Last modified: March 14, 2011
  use Parameter_Vars
  use Structures
  implicit none

  real shortwave           ! Declare function
  integer icell, iunit
  real accum_snow, add_to_soil, sublimated, snow_total, temp_avg

  iunit = Rng(icell)%range_type
  accum_snow = 0.0
  add_to_soil = 0.0
  sublimated = 0.0
  snow_total = 0.0
  Rng(icell)%ppt_soil = Globe(Rng(icell)%x, Rng(icell)%y)%precip
  temp_avg = ( Globe(Rng(icell)%x, Rng(icell)%y)%max_temp + Globe(Rng(icell)%x, Rng(icell)%y)%min_temp ) / 2.0      ! The method used in Century
 
  ! Judge whether precipitation is snow or liquid
  if (temp_avg .le. 0.0) then
    Rng(icell)%snow = Rng(icell)%snow + Globe(Rng(icell)%x, Rng(icell)%y)%precip           ! Recall snowpack is water equivalent
    accum_snow = Globe(Rng(icell)%x, Rng(icell)%y)%precip                                  ! Track snow accumulation
    Rng(icell)%ppt_soil = 0.0                                                              ! No water left to move into soil
  end if
  ! Add rain-on-snow to snowpack liquid 
  if (Rng(icell)%snow .gt. 0.0) then
    Rng(icell)%snow_liquid = Rng(icell)%snow_liquid + Rng(icell)%ppt_soil
    Rng(icell)%ppt_soil = 0.0
  end if
 
  ! Evaporate water from the snowpack                           
  if (Rng(icell)%snow .gt. 0.0) then
    ! Calculate cm of snow that remaining PET energy can evaporate
    sublimated = Rng(icell)%pet_remaining * 0.87                    ! 0.87 relates to heat of fusion for ice versus liquid water
    ! Calculate total snowpack water, ice and liquid
    snow_total = Rng(icell)%snow + Rng(icell)%snow_liquid
    if (sublimated .gt. snow_total) then
      sublimated = snow_total                                       ! Don't sublimate more than is present
    end if
    if (sublimated .lt. 0.0) then
      sublimated = 0.0
    end if
    ! Take sublimation from snow and snow liquid in proportion
    Rng(icell)%snow = Rng(icell)%snow - ( sublimated * (Rng(icell)%snow / snow_total ) )
    Rng(icell)%snow_liquid = Rng(icell)%snow_liquid - ( sublimated * (Rng(icell)%snow_liquid / snow_total) )       ! Snow_total cannot be zero, but may be very small.  A problem?)
    Rng(icell)%evaporation = Rng(icell)%evaporation + sublimated    ! Accumulate sublimated snow
    ! Decrement remaining PET by the energy that was used to evaporate snow
    Rng(icell)%pet_remaining = Rng(icell)%pet_remaining - ( sublimated / 0.87 )
    if (Rng(icell)%pet_remaining .lt. 0.0) then
      Rng(icell)%pet_remaining = 0.0
    end if
  end if

  ! Melt snow if the temperature is high enough
  if (Rng(icell)%snow .gt. 0.0 .and. temp_avg .ge. Parms(iunit)%melting_temp) then
    Rng(icell)%melt = Parms(iunit)%melting_slope * (temp_avg - &
                      Parms(iunit)%melting_temp) * shortwave(icell)
    if (Rng(icell)%melt .lt. 0.0) then
      Rng(icell)%melt = 0.0
    end if
    if ((Rng(icell)%snow - Rng(icell)%melt) .gt. 0.0) then
      Rng(icell)%snow = Rng(icell)%snow - Rng(icell)%melt
    else
      Rng(icell)%melt = Rng(icell)%snow 
      Rng(icell)%snow = 0.0
    end if
    ! Melting snow is liquid and drains excess
    Rng(icell)%snow_liquid = Rng(icell)%snow_liquid + Rng(icell)%melt
    ! Drain snowpack to 50% liquid content, and excess to soil.
    if (Rng(icell)%snow_liquid .gt. (0.5 * Rng(icell)%snow)) then
      add_to_soil = Rng(icell)%snow_liquid - ( 0.5 * Rng(icell)%snow )
      Rng(icell)%snow_liquid = Rng(icell)%snow_liquid - add_to_soil
      Rng(icell)%ppt_soil = Rng(icell)%ppt_soil + add_to_soil
      ! Return drained water into the soil
      Rng(icell)%melt = Rng(icell)%ppt_soil
    end if
  end if

end subroutine


real function shortwave(icell)
!**** Calculates the short wave radiation outside the atmosphere using Pennman's equation (1948).
!**** R. Boone, almost directly from CENTURY.    Last modified: September 25, 2010
  use Parameter_Vars
  use Structures
  implicit none

  real ahou, declination, par1, par2, radians_lat, solar_radiation
  real temp
  integer icell
  ! In CENTURY, a transmission coefficient is set to 0.8 for each month throughout the year.  
  ! Given that it has probably been that for years, I will hardwire the value.
  
  ! Convert latitude of site to radians
  radians_lat = Globe(Rng(icell)%x, Rng(icell)%y)%latitude * ( PI / 180.0 )
  
  ! Calculate short wave solar radiation on a clear day using equation in Sellers (1965)
  declination = 0.401426 * sin(6.283185 * (real(JULIAN_DAY_MID(month))-77.0)/365.0)
  temp = 1.0 - (-tan(radians_lat) * tan(declination))**2
  if (temp .lt. 0.0) then
    temp = 0.0
  end if
  par1 = sqrt(temp)
  par2 = (-tan(radians_lat)*tan(declination))
  
  ahou = atan2(par1, par2)
  if (ahou .lt. 0.0) then
    ahou = 0.0
  end if
  
  solar_radiation = 917.0 * 0.8 * (ahou * sin(radians_lat) * sin(declination) + cos(radians_lat) * cos(declination) * sin(ahou))
  
  shortwave = solar_radiation / 0.8  

end function


subroutine P_Evap (icell)
!**** Calculate Penmon-Monteith potential evapotraspiration for all the rangeland cells
!**** NOTE:  Using a subroutine here, rather than function call, given the structures used in GRange.
!**** 
!**** R. Boone, drawing from Century 4.5 PEvap.f    Last modified: March 14, 2011
   use Parameter_Vars
   use Structures
   implicit none

   integer  :: icell, nx, ny
   real     :: shortwave
   real     :: day_pet, month_pet, temp_range, temp_mean, const1, const2, langley2watts, site_latitude, fwloss_4

   const1 = 0.0023
   const2 = 17.8
   langley2watts = 54.0
   fwloss_4 = 0.8               ! A variable in Century, but defaults to 0.8

   ! Calculate PET for the Julian day in the middle of the current month
   nx = Rng(icell)%x
   ny = Rng(icell)%y
   site_latitude = Globe(nx,ny)%latitude
   temp_range = Globe(nx,ny)%max_temp - Globe(nx,ny)%min_temp
   temp_mean  = ( Globe(nx,ny)%max_temp + Globe(nx,ny)%min_temp ) / 2.0
   day_pet = ( const1 * ( temp_mean + const2 ) * sqrt(temp_range) * ( shortwave(icell) / langley2watts ) )
   ! Calculate monthly PET and convert to cm
   month_pet = ( day_pet * 30. ) / 10.
   if (month_pet .lt. V_LARGE .and. month_pet .ge. 0.5) then
     ! Do nothing.   Trying to avoid -NaN by using else contents
   else
     month_pet = 0.5
   end if
     
   ! Modified by FWLoss_4, the scaling factor for potential evapotranspiration
   Rng(icell)%pot_evap = month_pet * fwloss_4

end subroutine


subroutine Day_Length (icell)
!*** Calculate day length.  Rather than passing month, site latitude, and daylength (as a function), the information 
!*** is stored in structures.   Original was in C.  Modified as needed (e.g., different array base).  See DAYLEN.C
!***
!*** NOTE: This function also resets heat accumulation at the appropriate month, based on when the day length is at a minimum. 
!***
!***
!*** R. Boone   Last updated: January 2, 2014 
   use Parameter_Vars
   use Structures
   implicit none

   real     :: radians_lat, temp_1, adelt, temp_2, ahou
!   real     :: alat, dayl, dec, sinld, cosld, aob
   integer  :: icell
    
   ! Convert latitude of site to radians
   radians_lat = Globe(Rng(icell)%x, Rng(icell)%y)%latitude * ( PI / 180.0 )

   temp_1 = 2.0 * PI * (JULIAN_DAY_START(month) - 77.0) / 365.0
   adelt = 0.4014 * sin(temp_1)
   temp_1 = 1.0 - ( -tan(radians_lat) * adelt )**2.0
   if (temp_1 < 0.0) then
     temp_1 = 0.0
   end if
   temp_1 = sqrt(temp_1)
   temp_2 = -tan(radians_lat) * tan(adelt)
   ahou = atan2(temp_1, temp_2)
   Rng(icell)%day_length = ( ahou / PI ) * 24.0

   ! Set day length for this and the previous month, to be able to judge if seasons are changing.
   if (Rng(icell)%day_length .gt. Rng(icell)%last_month_day_length) then
     if (Rng(icell)%day_length_increasing .eq. .FALSE.) then  
       ! Starting a new ecological year, so to speak.  The days started getting longer this month.
       Rng(icell)%heat_accumulation = 0.0
     end if
     Rng(icell)%day_length_increasing = .TRUE.
   else
     Rng(icell)%day_length_increasing = .FALSE.
     ! The season is changing, it is the middle of the local summer.   Hopefully not too early to shift snow and snow liquid to the long-term storage
     Rng(icell)%old_snow = Rng(icell)%old_snow + Rng(icell)%snow
     Rng(icell)%snow = 0.0
     Rng(icell)%old_snow_liquid = Rng(icell)%old_snow_liquid + Rng(icell)%snow_liquid
     Rng(icell)%snow_liquid = 0.0
   end if
   Rng(icell)%last_month_day_length = Rng(icell)%day_length
  
end subroutine


