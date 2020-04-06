subroutine Water_Loss(icell)
!**** The primary water submodel, based on H2OLOSS in CENTURY.
!****
!**** R. Boone   Last modified: July 15, 2014

  use Parameter_Vars
  use Structures
  implicit none

  real     sd, avg_int, avg_bare_soil_evap, evl, add_water, trap, temp_avg, pttr
  real     base_flow, afl, strm, base, asimx, avail_water, avhsm
  real     rwc1, tot, tot2, awwt(SOIL_LAYERS), avinj, trl, fwlos, evmt, evlos
  real     litterd, avg_live_biomass, avg_dead_biomass, avg_litter_biomass
  real     fwloss_1, fwloss_2
  integer  icell, soil_layer, iunit

  iunit = Rng(icell)%range_type
  fwloss_1 = 1.0                           ! A variable in CENTURY, but defaults to 1.0
  fwloss_2 = 1.0                           ! A variable in CENTURY, but defaults to 1.0.  Include as needed.
  Rng(icell)%water_available(1) = 0.0      ! Zero-ed out at the beginning of each pass, as in Century
  Rng(icell)%water_available(2) = 0.0
  Rng(icell)%water_available(3) = 0.0

  ! The following is in Cycle, prior to the H2OLoss call, so moved to here.  The co2 production value (CO2CTR) is calculated in co2eff.f in CENTURY.
  ! Rng(icell)%co2_value(ifacet) = Parms(iunit)%effect_of_co2_on_production(ifacet) ... wrong ... need CO2 coming in, not this.
  ! CO2 is not used directly.   Instead, the CO2 EFFECT ON PRODUCTION is used, and only that.   So if production increases, that increases.   I.E.  It is at the landscape unit level.
  if (Rng(icell)%melt .gt. 0.0) then
    Rng(icell)%ratio_water_pet = Rng(icell)%melt / Rng(icell)%pot_evap
  else
    Rng(icell)%ratio_water_pet = ( Rng(icell)%water_available(3) + Globe(Rng(icell)%x, Rng(icell)%y)%precip ) / &
                                   Rng(icell)%pot_evap    ! Excludes irrigation
  endif
  ! Concludes the part from Cycle prior to the call to H2OLoss.

  trap = 0.01
  add_water = 0.0
  base_flow = 0.0
  base = 0
  strm = 0.0
  asimx = 0.0
  Rng(icell)%transpiration = 0.0
  Rng(icell)%evaporation = 0.0
  temp_avg = ( Globe(Rng(icell)%x, Rng(icell)%y)%max_temp + Globe(Rng(icell)%x, Rng(icell)%y)%min_temp ) / 2.0      ! The method used in Century

  ! Skipping irrigation, so "inputs" in h2oloss is just precipitation

  ! PET remaining stores the energy available to evaporate after each step that saps energy
  Rng(icell)%pet_remaining = Rng(icell)%pot_evap

  ! Set the snowpack, snow melt, and sublimation
  ! call Snow_Dynamics(icell) ... NO, already handled in Weather Update

  ! Calculate runoff using Probert et al. (1995).  See CENTURY H2OLoss for full citation.
  Rng(icell)%runoff = MAX(0.0, Parms(iunit)%prcp_threshold_fraction * (Rng(icell)%ppt_soil - Parms(iunit)%prcp_threshold))
  Rng(icell)%ppt_soil = Rng(icell)%ppt_soil - Rng(icell)%runoff                           ! PPT_SOIL from SNOWCENT is WINPUTS in H2OLOSS.  It is the precipitation that reaches the soil, versus that which is converted to snow.

  ! Compute bare soil water loss and interception, when there is no snow.
  ! (The following line was changed from CENTURY to avoid an equality test on 0.0)
  if (Rng(icell)%snow .le. 0.00001) then
    avg_live_biomass = Rng(icell)%total_aground_live_biomass
    avg_dead_biomass = Rng(icell)%total_bground_live_biomass

    avg_litter_biomass = ( Rng(icell)%litter_structural_carbon(SURFACE_INDEX) + &
                           Rng(icell)%litter_metabolic_carbon(SURFACE_INDEX) ) * 2.5

    sd = avg_live_biomass + avg_dead_biomass
    ! The following were modified to use modern functions, and the 0 test added to avoid infinite results in avg_bare_soil_evap.  In short, biomass should not be negative.  Going negative suggests a problem elsewhere, but one debug at a time.
    !    if (sd .gt. 800.0) sd = 800.0
    sd = min(800.0, sd)
    sd = max(sd, 0.0)
    if (avg_litter_biomass .gt. 400.0) then
      litterd = 400.0
    else
      litterd = avg_litter_biomass                                                   ! Avg_litter_biomass is not used after this point
    end if

    ! Calculate canopy interception
    avg_int = (0.0003 * litterd + 0.0006 * sd) * fwloss_1     ! Not stored long-term until shown as needed
    ! Calculate bare soil evaporation
    avg_bare_soil_evap = 0.5 * exp((-0.002 * litterd) - (0.004 * sd)) * fwloss_2     ! Not stored long-term so far.
    ! Calculate total surface evaporation losses.  The maximum allowable is 0.4 * PET
    evl = MIN(((avg_bare_soil_evap + avg_int) * Rng(icell)%ppt_soil), (0.4 * Rng(icell)%pet_remaining))
    Rng(icell)%evaporation = Rng(icell)%evaporation + evl
    ! Calculate remaining water to add to soil and potential transpiration as remaining pet
    add_water = Rng(icell)%ppt_soil - evl
    add_water = max(0.0,add_water)
    trap = Rng(icell)%pet_remaining - evl
  else
    avg_bare_soil_evap = 0.0                                                         ! Initialized in H2OLOS, here in if-else, but functionally the same.
    trap = 0.0
    add_water = 0.0
    evl = 0.0
  end if

  ! Determine potential transpiration water loss (cm/mon) as a function of precipitation and live biomass.
  ! If the temperature is less than 2 degrees C, transpiration will be turned off.
  if (temp_avg .lt. 2.0) then
    pttr = 0.0
  else
    pttr = Rng(icell)%pet_remaining * 0.65 * (1.0 - exp(-0.020 * avg_live_biomass)) * &
      ( ( Rng(icell)%co2_effect_on_production(H_FACET) +  Rng(icell)%co2_effect_on_production(S_FACET) + &
          Rng(icell)%co2_effect_on_production(T_FACET) ) / 3 )      ! NOT facet-based.  Water loss will be based on the average of the
                                                                   ! three facets.
  end if
  if (pttr .le. trap) trap = pttr
  if (trap .le. 0.0) trap = 0.01

  ! CENTURY maintains pttr for work on harvest.  I don't think we need that here now.  REVISIT

  ! Calculate potential evapotranspiration rate from the top soil layer (cm/day).  This is not subtracted until
  ! after transpiration losses are calculated.
  Rng(icell)%pet_top_soil = Rng(icell)%pet_remaining - trap - evl
  if (Rng(icell)%pet_top_soil .lt. 0.0) Rng(icell)%pet_top_soil = 0.0

  ! Transpire water from added water, then pass it on to soil.
  Rng(icell)%transpiration = MIN((trap - 0.01), add_water)
  Rng(icell)%transpiration = MAX(0.0, Rng(icell)%transpiration)
  trap = trap - Rng(icell)%transpiration
  add_water = add_water - Rng(icell)%transpiration

  ! Add water to the soil, including base flow and storm flow
  ! stream_water = 0.0
  base_flow = 0.0
  ! Rng(icell)%stream(1) = 0.0

  do soil_layer=1,SOIL_LAYERS
    Rng(icell)%asmos(soil_layer) = Rng(icell)%asmos(soil_layer) + add_water          ! Add water to layer
    ! Calculate field capacity of soil, drain, and pass excess on to amov
    afl = Rng(icell)%soil_depth(soil_layer) * Rng(icell)%field_capacity(soil_layer)
    if (Rng(icell)%asmos(soil_layer) .gt. afl) then
      Rng(icell)%amov(soil_layer) = Rng(icell)%asmos(soil_layer) - afl
      Rng(icell)%asmos(soil_layer) = afl
      ! If the bottom layer, the remainder is storm flow
      if (soil_layer .eq. SOIL_LAYERS) strm = Rng(icell)%amov(soil_layer) * Rng(icell)%storm_flow
    else
      Rng(icell)%amov(soil_layer) = 0.0
    end if
    add_water = Rng(icell)%amov(soil_layer)
  end do

  ! Compute base flow and stream flow
  Rng(icell)%holding_tank = Rng(icell)%holding_tank + add_water - strm
  ! Drain base flow fraction from holding tank
  base = Rng(icell)%holding_tank * Parms(iunit)%base_flow_fraction
  Rng(icell)%holding_tank = Rng(icell)%holding_tank - base
  ! Add runoff to stream flow
  ! Rng(icell)%stream(1) = strm + base + Rng(icell)%runoff          Commenting out stream flow modeling for now.  Not needed for now.
  ! Save asmos(1) before transpiration
  asimx = Rng(icell)%asmos(1)

  ! Calculate transpiration water loss
  rwc1 = 0.0
  tot = 0.0
  tot2 = 0.0
  do soil_layer=1,SOIL_LAYERS
    avail_water = Rng(icell)%asmos(soil_layer) - Rng(icell)%wilting_point(soil_layer) * Rng(icell)%soil_depth(soil_layer)
    if (avail_water .lt. 0.0) avail_water = 0.0
    ! Calculate available water weighted by transpiration depth distribution factors
    awwt(soil_layer) = avail_water * Parms(iunit)%soil_transpiration_fraction(soil_layer)
    tot = tot + avail_water
    tot2 = tot2 + awwt(soil_layer)
    ! Moving up a copy of this update, which will do no harm, and fix a problem if the following clause is skipped (i.e.,
    !  if tot2 .le. 0
    Rng(icell)%relative_water_content(soil_layer) = (Rng(icell)%asmos(soil_layer) / Rng(icell)%soil_depth(soil_layer) - &
        Rng(icell)%wilting_point(soil_layer)) / ( Rng(icell)%field_capacity(soil_layer) - Rng(icell)%wilting_point(soil_layer))
  end do

  ! Calculate transpiration water loss (cm/month)
  trap = MIN(tot, trap)

  ! CENTURY calculates the layers providing water, based on crops, trees, and total leaf area.  Another
  ! simplification here.  All layers will provide water, given that they go down to 60 cm.  In short,
  ! I'm not simulating crops, and trees are a secondary interest.  REVISIT as necessary.
  ! But ... this was the way it was done in Century until 2003.

  if (tot2 .gt. 0.0) then
    do soil_layer=1,SOIL_LAYERS
      avinj = Rng(icell)%asmos(soil_layer) - Rng(icell)%wilting_point(soil_layer) * Rng(icell)%soil_depth(soil_layer)
      if (avinj .lt. 0.0) avinj = 0.0
      ! Calculate transpiration loss from soil layer, using weighted availabilities
      trl = (trap * awwt(soil_layer)) / tot2
      if (trl .gt. avinj) trl = avinj
      Rng(icell)%asmos(soil_layer) = Rng(icell)%asmos(soil_layer) - trl
      avinj = avinj - trl
      Rng(icell)%transpiration = Rng(icell)%transpiration + trl
      Rng(icell)%relative_water_content(soil_layer) = (Rng(icell)%asmos(soil_layer) / Rng(icell)%soil_depth(soil_layer) - &
           Rng(icell)%wilting_point(soil_layer)) / ( Rng(icell)%field_capacity(soil_layer) - Rng(icell)%wilting_point(soil_layer))
      if (Rng(icell)%relative_water_content(soil_layer) .gt. V_LARGE) then
        Rng(icell)%relative_water_content(soil_layer) = V_LARGE
      end if
      if (Rng(icell)%relative_water_content(soil_layer) .lt.   0.0  ) then
        Rng(icell)%relative_water_content(soil_layer) = 0.0
      end if

      ! Calculate water available to plants for growth
      Rng(icell)%water_available(1) = Rng(icell)%water_available(1) + avinj
      if (soil_layer .le. 2) Rng(icell)%water_available(3) = Rng(icell)%water_available(3) + avinj
    end do
  end if

  ! Sum water available for plants to survive.
  do soil_layer=1,SOIL_LAYERS
    Rng(icell)%water_available(2) = Rng(icell)%water_available(2) + avinj
  end do
  ! Ignoring entry regarding harvesting of crops ... REVISIT?

  ! Minimum relative water content for top layer to evaporate  (MANY OF THESE COMMENTS COME DIRECTLY FROM CENTURY)
  fwlos = 0.25

  ! Fraction of water content between FWLOS and field capacity
  evmt = (Rng(icell)%relative_water_content(1) - fwlos) / (1. - fwlos)
  if (evmt .lt. 0.01) evmt = 0.01

  ! Evaporation loss from layer 1
  evlos = evmt * Rng(icell)%pet_top_soil * avg_bare_soil_evap * 0.10
  avinj = Rng(icell)%asmos(1) - Rng(icell)%wilting_point(1) * Rng(icell)%soil_depth(1)
  if (avinj .lt. 0.0) avinj = 0.0
  if (evlos .gt. avinj) evlos = avinj
  Rng(icell)%asmos(1) = Rng(icell)%asmos(1) - evlos
  Rng(icell)%evaporation = Rng(icell)%evaporation + evlos

  ! Recalculate RWCF(1) to estimate mid-month water content
  avhsm = ( Rng(icell)%asmos(1) + Rng(icell)%relative_water_content(1) * asimx ) / ( 1. + Rng(icell)%relative_water_content(1) )
  Rng(icell)%relative_water_content(1) = ( avhsm / Rng(icell)%soil_depth(1) - &
       Rng(icell)%wilting_point(1)) / ( Rng(icell)%field_capacity(1) - Rng(icell)%wilting_point(1))
  if (Rng(icell)%relative_water_content(1) .gt. V_LARGE) then
    Rng(icell)%relative_water_content(1) = V_LARGE
    write(ECHO_FILE,*) 'Relative water content reset to very large in Water_Loss, layer 1: ',icell
  end if
  if (Rng(icell)%relative_water_content(1) .lt.   0.0  ) then
    Rng(icell)%relative_water_content(1) = 0.0
    write(ECHO_FILE,*) 'Relative water content reset to 0.0 in Water_Loss, layer 1: ',icell
  end if

  ! Update water available pools minus evaporation from top layer
  Rng(icell)%water_available(1) = Rng(icell)%water_available(1) - evlos
  Rng(icell)%water_available(1) = max(0.0, Rng(icell)%water_available(1))
  Rng(icell)%water_available(2) = Rng(icell)%water_available(2) - evlos
  Rng(icell)%water_available(2) = max(0.0, Rng(icell)%water_available(2))
  Rng(icell)%water_available(3) = Rng(icell)%water_available(3) - evlos
  Rng(icell)%water_available(3) = max(0.0, Rng(icell)%water_available(3))

  ! Compute annual actual evapotranspiration
  Rng(icell)%annual_evapotranspiration = Rng(icell)%annual_evapotranspiration + Rng(icell)%evaporation + Rng(icell)%transpiration

  ! From CYCLE, following the call to H2OLoss, a single if-then, which will be incorporated here.
  if (Rng(icell)%snow .gt. 0.) then
    Rng(icell)%soil_surface_temperature = 0.0
  endif
  ! Concludes the material from CYCLE that relates (sort of) to water loss.

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'H2O_LSS')

end subroutine



real function line(x, x1, y1, x2, y2)
!**** The function for a line, yielding a Y value for a given X and the parameters for a line.
!****
!**** R Boone from CENTURY LINE function.   Last modified: September 25, 2010
  real x, x1, y1, x2, y2

  line = (y2 - y1) / (x2 - x1) * (x - x2) + y2

  return

end function
