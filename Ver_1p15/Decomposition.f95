subroutine Decomposition(icell)
!**** The main decomposition model, which is a simplified (or really, very similar version of) version of CENTURY that M. Coughenour adapted for use in Savanna.
!****
!**** R. Boone   Last modified: April 6, 2020
  use Parameter_Vars
  use Structures
  implicit none

  real decinv(2), prf(3), immobil, ps2s3(2), ps1s3(2)
  real peftxa, peftxb, fxmca, fxmcb, fxmxs, eftext, decart, tndec, tcdec, co2los, fps1s3, tosom1net, tosom2net, tosom3net
  real p1co2, demnsom2, demnsom3, demn, somtomin, smintosom, rces1, rces2, rces3, fr, grmin
  real dec2rt, fps2s3, dec3rt, decc, decn, volex, temper, decmrt, demnsom1, fxbiom, fwdfx
  real fixa
  real tbio, biof, fixs
  real decodt
  real dthc, dthn, frac_lignin

  integer icell, iunit, ifacet, ilayer

  decodt = 1./12.        ! 48.            ! 12 months, 4 calls to decomposition each month, in CENTURY and in SAVANNA.  I am going to experiment with a single call, until it proves ineffective.  Multiple calls in a monthly model is not clear.  In Savanna's weekly model it is clearer.
  fxmca = -0.125
  fxmcb = 0.005
  fxmxs = 0.35
  ps1s3(1) = 0.003
  ps1s3(2) = 0.032
  ps2s3(1) = 0.003
  ps2s3(2) = 0.009
  peftxa = 0.25
  peftxb = 0.75
  rces1 = 13.0               ! Default values for rces in Savanna
  rces2 = 18.0
  rces3 = 8.0

  iunit = Rng(icell)%range_type

  ! nsoiltp is soil type in Savanna.  I have access to that through the structures
  !**** NOTE that in Savanna, the tree facet is 2, the shrub facet is 3.  The opposite here.  *****

  ! Partition some of litter fall among facets.  Calculated elsewhere ... so commented out.
  ! Rng(icell)%herb_cover = 1. - Rng(icell)%woody_cover - Rng(icell)%shrub_cover

  do ifacet = 1, FACETS
    ! decomp_litter_mix_facets determines mixing, litter flows to facets in proportion to their cover. 1=in proportion to their cover 0=only to the current facet.
    select case (ifacet)
      case (H_FACET)
        prf(T_FACET) = Rng(icell)%facet_cover(T_FACET) * Parms(iunit)%decomp_litter_mix_facets
        prf(S_FACET) = Rng(icell)%facet_cover(S_FACET) * Parms(iunit)%decomp_litter_mix_facets
        prf(H_FACET) = 1. - prf(S_FACET) - prf(T_FACET)
      case (S_FACET)
        prf(H_FACET) = Rng(icell)%facet_cover(H_FACET) * Parms(iunit)%decomp_litter_mix_facets
        prf(T_FACET) = Rng(icell)%facet_cover(T_FACET) * Parms(iunit)%decomp_litter_mix_facets
        prf(S_FACET) = 1. - prf(H_FACET) - prf(T_FACET)
      case (T_FACET)
        prf(H_FACET) = Rng(icell)%facet_cover(H_FACET) * Parms(iunit)%decomp_litter_mix_facets
        prf(S_FACET) = Rng(icell)%facet_cover(S_FACET) * Parms(iunit)%decomp_litter_mix_facets
        prf(T_FACET) = 1. - prf(H_FACET) - prf(S_FACET)
    end select

    ! Cut-off litter input into a facet with 0 cover
    if (Rng(icell)%facet_cover(T_FACET) .lt. 0.01 .and. ifacet .ne. T_FACET) then
      prf(T_FACET) = 0.
      prf(H_FACET) = prf(H_FACET) + Rng(icell)%facet_cover(T_FACET)
    endif
    if (Rng(icell)%facet_cover(S_FACET) .lt. 0.01 .and. ifacet .ne. S_FACET) then
      prf(S_FACET) = 0.
      prf(H_FACET) = prf(H_FACET) + Rng(icell)%facet_cover(S_FACET)
    endif

    ! Loop skipped that deals with species.

    ! Note that units in this section may be gC/m^2 or gB/m^2, so use caution.
    ! The array index here is:  1 - Phenological death, 2 - Incremental death, 3 - Herbivory, 4 - Fire  (NO IT DOESN'T, RIGHT NOW)
    ! Fire is not transfered here
    if (prf(ifacet) .gt. 0.) then
      ! do itype=1,3    (not used presently)
      ! Fine roots
      if (Rng(icell)%dead_fine_root_carbon(ifacet) .gt. 0.0) then
        dthc = Rng(icell)%dead_fine_root_carbon(ifacet)    ! * prf(ifacet)
      else
        dthc = 0.0
      end if
      if (Rng(icell)%dead_fine_root_nitrogen(ifacet) .gt. 0.0) then
        dthn = Rng(icell)%dead_fine_root_nitrogen(ifacet)  ! * prf(ifacet)
      else
        dthn = 0.0
      end if
      frac_lignin = Rng(icell)%lignin_fine_root(ifacet) / Rng(icell)%fine_root_carbon(ifacet)
      frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
      frac_lignin = min(0.50, frac_lignin)
      call Partition_Litter (icell, SOIL_INDEX, dthc, dthn, frac_lignin, 'froots')
      ! Standing dead, leaves and stems
      !    Standing dead is already partitioned in PLANT_DEATH.   That includes LEAVES (plus shoots)
      ! Seed
      if (Rng(icell)%dead_seed_carbon(ifacet) .gt. 0.0) then
        dthc = Rng(icell)%dead_seed_carbon(ifacet)     ! * prf(ifacet)
      else
        dthc = 0.0
      end if
      if (Rng(icell)%dead_seed_nitrogen(ifacet) .gt. 0.0) then
        dthn = Rng(icell)%dead_seed_nitrogen(ifacet)   ! * prf(ifacet)
      else
        dthn = 0.0
      end if
        ! Using leaf lignin for simplicity
      frac_lignin = Rng(icell)%lignin_leaf(ifacet) / Rng(icell)%leaf_carbon(ifacet)
      frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
      frac_lignin = min(0.50, frac_lignin)
      call Partition_Litter (icell, SURFACE_INDEX, dthc, dthn, frac_lignin, 'seeds')
      ! Leaf                                          AN ACCUMULATOR ONLY - DEAD LEAF GOES TO STANDING DEAD THEN TO LITTER THEN DECOMPOSITION
      ! Fine branch and stem wood
      dthc = Rng(icell)%dead_fine_branch_carbon(ifacet)      ! * prf(ifacet)
      dthn = Rng(icell)%dead_fine_branch_nitrogen(ifacet)    ! * prf(ifacet)
      Rng(icell)%dead_total_fine_branch_carbon = Rng(icell)%dead_total_fine_branch_carbon + dthc          ! Accumulating values for shrubs and trees.
      Rng(icell)%dead_total_fine_branch_nitrogen = Rng(icell)%dead_total_fine_branch_nitrogen + dthn
      frac_lignin = Rng(icell)%lignin_fine_branch(ifacet) / Rng(icell)%fine_branch_carbon(ifacet)
      frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
      frac_lignin = min(0.50, frac_lignin)
      call Partition_Litter (icell, SURFACE_INDEX, dthc, dthn, frac_lignin, 'cbrnch')
      ! Coarse roots
      dthc = Rng(icell)%dead_coarse_root_carbon(ifacet)      ! * prf(ifacet)
      dthn = Rng(icell)%dead_coarse_root_nitrogen(ifacet)    ! * prf(ifacet)
      Rng(icell)%dead_total_coarse_root_carbon = Rng(icell)%dead_total_coarse_root_carbon + dthc          ! Accumulating values for shrubs and trees.  Be sure they are partitioned when done.
      Rng(icell)%dead_total_coarse_root_nitrogen = Rng(icell)%dead_total_coarse_root_nitrogen + dthn
      ! Do flows to litter, keeping track of structural and metabolic components
      frac_lignin = Rng(icell)%lignin_coarse_root(ifacet) / Rng(icell)%coarse_root_carbon(ifacet)
      frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
      frac_lignin = min(0.50, frac_lignin)
      call Partition_Litter (icell, SOIL_INDEX, dthc, dthn, frac_lignin, 'croots')
      ! Coarse branches
      dthc = Rng(icell)%dead_coarse_branch_carbon(ifacet)    ! * prf(ifacet)
      dthn = Rng(icell)%dead_coarse_branch_nitrogen(ifacet)  ! * prf(ifacet)
      Rng(icell)%dead_total_coarse_branch_carbon = Rng(icell)%dead_total_coarse_branch_carbon + dthc
      Rng(icell)%dead_total_coarse_branch_nitrogen = Rng(icell)%dead_total_coarse_branch_nitrogen + dthn
      ! Do flows to litter, keeping track of structural and metabolic components
      ! Coarse branches cannot be part of dead standing carbon. It could be, but functionally coarse branches are not useful for herbivores, for example
      frac_lignin = Rng(icell)%lignin_coarse_branch(ifacet) / Rng(icell)%coarse_branch_carbon(ifacet)
      frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
      frac_lignin = min(0.50, frac_lignin)
      call Partition_Litter (icell, SURFACE_INDEX, dthc, dthn, frac_lignin, 'cbranch')

      ! end do                    ! End type of mortality (e.g., phenology, fire, herbivory)
    end if                    ! End if any dead carbon flow
  end do                  ! End facet

  ! Calculate effect of attributes of cell on decomposition.  These values are in DECOMP in Savanna directly, but I had programmed them drawing from
  ! CENTURY, and so will use those values.  They should be essentially directly replacable.
  call Effects_on_Decomposition(icell)

  do ilayer=SURFACE_INDEX, SOIL_INDEX               ! Confirmed two layers
    Rng(icell)%tnetmin(ilayer) = 0.
    Rng(icell)%tminup(ilayer) = 0.
    Rng(icell)%grossmin(ilayer) = 0.

    ! SOM1 (fast microbial carbon) decomposition goes to som2 and som3, but some is ignored for now.
    ! Note:  Som1 n/c ratio is always higher than som2 and some3, always mineralize n
    eftext = peftxa + peftxb * Rng(icell)%sand(ilayer)
    decart = Parms(iunit)%decomp_rate_fast_som(ilayer) * Rng(icell)%all_effects_on_decomp * eftext * decodt
    tcdec = Rng(icell)%fast_soil_carbon(ilayer) * decart

    ! Respiration loss
    select case (ilayer)
      case (SURFACE_INDEX)
        p1co2 = 0.6
      case (SOIL_INDEX)
        p1co2 = 0.17 + 0.68 * Rng(icell)%sand(2)                                               ! Second layer of sand used.  Sand is calculated in Savanna.  We have that directly, if I am interpreting it correctly.
    end select
    co2los = tcdec * p1co2                                                                     ! Fixed - June 2014

    ! Net flow to SOM3 (passive_soil_carbon)
    fps1s3 = ps1s3(1) + ps1s3(2) * Rng(icell)%clay(2)                                                ! Second layer of clay used.
    tosom3net = tcdec * fps1s3 * ( 1.0 + 5.0 * ( 1. - Rng(icell)%anerobic_effect_on_decomp ) )

    ! Net flow to SOM2 (intermediate_soil_carbon)
    tosom2net = tcdec - co2los - tosom3net
    tosom2net = max(0., tosom2net)

    if (tcdec .gt. 0.) then
      if (Rng(icell)%fast_soil_carbon(ilayer) .gt. 0.0001) then
        tndec = tcdec * Rng(icell)%fast_soil_nitrogen(ilayer) / Rng(icell)%fast_soil_carbon(ilayer)
      else
        tndec = 0.0
      end if
      demnsom2 = tosom2net / rces2
      demnsom3 = tosom3net / rces3

      demn = demnsom2 + demnsom3

      if (tndec .gt. demn) then
        somtomin = tndec - demn
        smintosom = 0.
      else
        smintosom = demn - tndec
        somtomin = 0.
        if (smintosom .gt. 0.00001) then
          if (smintosom .gt. Rng(icell)%mineral_nitrogen(ilayer)) then
            fr = Rng(icell)%mineral_nitrogen(ilayer) / smintosom                       ! Division by zero prevented two lines above
            fr = max(0., fr)
            smintosom = Rng(icell)%mineral_nitrogen(ilayer)
            smintosom = max(0., smintosom)
            tcdec = tcdec * fr
            tndec = tndec * fr
            tosom2net = tosom2net * fr
            tosom3net = tosom3net * fr
            demnsom2 = demnsom2 * fr
            demnsom3 = demnsom3 * fr
          end if
        end if
      end if

      ! Doing summaries that change the values stored in structures
      ! Note that only fast SOM and SOME have two layers
      Rng(icell)%fast_soil_carbon(ilayer) = Rng(icell)%fast_soil_carbon(ilayer) - tcdec
      Rng(icell)%intermediate_soil_carbon = Rng(icell)%intermediate_soil_carbon + tosom2net
      Rng(icell)%passive_soil_carbon = Rng(icell)%passive_soil_carbon + tosom3net
      Rng(icell)%fast_soil_nitrogen(ilayer) = Rng(icell)%fast_soil_nitrogen(ilayer) - tndec
      Rng(icell)%intermediate_soil_nitrogen = Rng(icell)%intermediate_soil_nitrogen + demnsom2
      Rng(icell)%passive_soil_nitrogen = Rng(icell)%passive_soil_nitrogen + demnsom3
      Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) + somtomin - smintosom
      Rng(icell)%grossmin(ilayer) = Rng(icell)%grossmin(ilayer) + tndec
      Rng(icell)%tnetmin(ilayer) = Rng(icell)%tnetmin(ilayer) + somtomin - smintosom
      if (Rng(icell)%tnetmin(ilayer) .lt. 0.0) Rng(icell)%tnetmin(ilayer) = 0.0
    end if

    ! Metabolic decomposition
    decmrt = Parms(iunit)%decomp_rate_metabolic_litter(ilayer) * Rng(icell)%all_effects_on_decomp * decodt
    tcdec = Rng(icell)%litter_metabolic_carbon(ilayer) * decmrt
    tcdec = min(tcdec, Rng(icell)%litter_metabolic_carbon(ilayer))
    tcdec = max(tcdec, 0.)
    if (Rng(icell)%litter_metabolic_carbon(ilayer) .gt. 0.0001) then
      tndec = tcdec * Rng(icell)%litter_metabolic_nitrogen(ilayer) / Rng(icell)%litter_metabolic_carbon(ilayer)
      tndec = min(tndec, Rng(icell)%litter_metabolic_nitrogen(ilayer))
      tndec = max(tndec, 0.)
    else
      tndec = 0.0
    end if

    ! 0.55 is the fraction respired, PMCO2 from Century
    co2los = tcdec * 0.55

    tosom1net = tcdec - co2los
    demnsom1 = tosom1net / rces1

    if (tndec .gt. demnsom1) then
      somtomin = tndec - demnsom1
      smintosom = 0.
    else                                                       ! Changed in Sav5 ... a multi-line addition
      somtomin = 0.
      smintosom = demnsom1 - tndec
      ! Reduce decomp if not enough mineral N
      if (smintosom .gt. Rng(icell)%mineral_nitrogen(ilayer) ) then
        fr = Rng(icell)%mineral_nitrogen(ilayer) / smintosom
        fr = max(0., fr)
        smintosom = Rng(icell)%mineral_nitrogen(ilayer)
        smintosom = max(0., smintosom)
        tcdec = tcdec * fr
        tndec = tndec * fr
        tosom1net = tosom1net * fr
        demnsom1 = demnsom1 * fr
      end if
    end if

    ! Doing summaries that change the values stored in structures
    ! Note that only fast SOM and SOME have two layers
    Rng(icell)%litter_metabolic_carbon(ilayer) = Rng(icell)%litter_metabolic_carbon(ilayer) - tcdec
    Rng(icell)%fast_soil_carbon(ilayer) = Rng(icell)%fast_soil_carbon(ilayer) + tosom1net
    Rng(icell)%litter_metabolic_nitrogen(ilayer) = Rng(icell)%litter_metabolic_nitrogen(ilayer) - tndec
    Rng(icell)%fast_soil_nitrogen(ilayer) = Rng(icell)%fast_soil_nitrogen(ilayer) + demnsom1
    Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) + somtomin - smintosom

    Rng(icell)%grossmin(ilayer) = Rng(icell)%grossmin(ilayer) + tndec
    Rng(icell)%tnetmin(ilayer) = Rng(icell)%tnetmin(ilayer) + somtomin - smintosom
    if (Rng(icell)%tnetmin(ilayer) .lt. 0.0) Rng(icell)%tnetmin(ilayer) = 0.0
    ! Structural decomposition, goes to SOM1 (fast) and SOM2 (intermediate)
    frac_lignin = ( Rng(icell)%plant_lignin_fraction(H_FACET,ilayer) + Rng(icell)%plant_lignin_fraction(S_FACET,ilayer) + &
                    Rng(icell)%plant_lignin_fraction(T_FACET,ilayer) ) / 3.0
    frac_lignin = max(0.02, frac_lignin)                                              ! From Century CmpLig.f
    frac_lignin = min(0.50, frac_lignin)
    call Track_Lignin(icell, ilayer, 0, Rng(icell)%litter_structural_carbon(ilayer), &
      Rng(icell)%litter_structural_nitrogen(ilayer), Rng(icell)%all_effects_on_decomp, &
      decodt, frac_lignin, grmin, immobil)
    Rng(icell)%grossmin(ilayer) = Rng(icell)%grossmin(ilayer) + grmin
    Rng(icell)%tnetmin(ilayer) = Rng(icell)%tnetmin(ilayer) - immobil
    Rng(icell)%tminup(ilayer) = Rng(icell)%tminup(ilayer) + immobil

    ! Intermediate pool (SOM2) decomposition, which goes to SOM1 (fast) and SOM3 (passive)
    if (ilayer .eq. SOIL_INDEX) then
      dec2rt = Parms(iunit)%decomp_rate_inter_som * Rng(icell)%all_effects_on_decomp * decodt
      tcdec = Rng(icell)%intermediate_soil_carbon * dec2rt
      if (tcdec .gt. 0.) then
        co2los = tcdec * 0.55 * Rng(icell)%anerobic_effect_on_decomp                  ! 0.55 is fraction respired (p2co2 in Century)

        ! Net flow to SOM3 (passive)
        fps2s3 = ps2s3(1) + ps2s3(2) * Rng(icell)%clay(ilayer)                     ! Using top two soil layers.  Note confusion of ilayer.
        tosom3net = tcdec * fps2s3 * (1.0 + 5.0 * (1.0 - Rng(icell)%anerobic_effect_on_decomp))
        ! Net flow to SOM1 (fast)
        tosom1net = tcdec - co2los - tosom3net
        ! N flows
        demnsom1 = tosom1net / rces1
        demnsom3 = tosom3net / rces3
        if (Rng(icell)%intermediate_soil_carbon .gt. 0.0) then
          tndec = tcdec * Rng(icell)%intermediate_soil_nitrogen / Rng(icell)%intermediate_soil_carbon
        else
          tndec = 0.0
        end if

        demn = demnsom1 + demnsom3

        if (tndec .gt. demn) then
          somtomin = tndec - demn
          smintosom = 0.
        else
          smintosom = demn - tndec
          somtomin = 0.
        end if

        ! Reduce decomposition if not enough mineral nitrogen to support it via immoblization
        if (smintosom .gt. 0.) then
          if (smintosom .gt. Rng(icell)%mineral_nitrogen(ilayer)) then                                     ! Check the use of ILAYER.  Caution not cited anymore.
            fr = Rng(icell)%mineral_nitrogen(ilayer) / smintosom                                           ! Checked for 0 above.
            fr = max(0., fr)
            smintosom = Rng(icell)%mineral_nitrogen(ilayer)
            smintosom = max(0., smintosom)
            tcdec = tcdec * fr
            tndec = tndec * fr
            tosom1net = tosom1net * fr
            tosom3net = tosom3net * fr
            demnsom1 = demnsom1 * fr
            demnsom3 = demnsom3 * fr
          end if
        end if

        ! Do the updates to the information in the structure
        Rng(icell)%fast_soil_carbon(ilayer) = Rng(icell)%fast_soil_carbon(ilayer) + tosom1net
        Rng(icell)%intermediate_soil_carbon = Rng(icell)%intermediate_soil_carbon - tcdec
        Rng(icell)%passive_soil_carbon = Rng(icell)%passive_soil_carbon + tosom3net

        Rng(icell)%intermediate_soil_nitrogen = Rng(icell)%intermediate_soil_nitrogen - tndec
        Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) - smintosom + somtomin
        Rng(icell)%fast_soil_nitrogen(ilayer) = Rng(icell)%fast_soil_nitrogen(ilayer) + demnsom1
        Rng(icell)%passive_soil_nitrogen = Rng(icell)%passive_soil_nitrogen + demnsom3

        Rng(icell)%grossmin(ilayer) = Rng(icell)%grossmin(ilayer) + tndec
        Rng(icell)%tnetmin(ilayer) = Rng(icell)%tnetmin(ilayer) + somtomin - smintosom
        if (Rng(icell)%tnetmin(ilayer) .lt. 0.0) Rng(icell)%tnetmin(ilayer) = 0.0
        Rng(icell)%tminup(ilayer) = Rng(icell)%tminup(ilayer) + smintosom

      end if

      ! SOM3 (passive component) decomposition
      dec3rt = Parms(iunit)%decomp_rate_slow_som * Rng(icell)%all_effects_on_decomp * decodt
      tcdec = Rng(icell)%passive_soil_carbon * dec3rt

      if (tcdec .gt. 0.) then
        co2los = tcdec * 0.55 * Rng(icell)%anerobic_effect_on_decomp                              ! 0.55 is fraction respired (p3co2 in Century)
        ! Net flow to SOM1
        tosom1net = tcdec - co2los
        ! N flows, mineralization, because N/C of SOM3 (1/rces3=7) is less than SOM1 (microb nc = 10)  (comment from Savanna)
        if (Rng(icell)%passive_soil_carbon .gt. 0.0) then
          tndec = tcdec * Rng(icell)%passive_soil_nitrogen / Rng(icell)%passive_soil_carbon
        else
          tndec = 0.0
        end if
        demnsom1 = tosom1net / rces1
        somtomin = tndec - demnsom1

        ! Do the updates to the main material in the strucutre.
        Rng(icell)%passive_soil_carbon = Rng(icell)%passive_soil_carbon - tcdec
        Rng(icell)%fast_soil_carbon(ilayer) = Rng(icell)%fast_soil_carbon(ilayer) + tosom1net
        Rng(icell)%fast_soil_nitrogen(ilayer) = Rng(icell)%fast_soil_nitrogen(ilayer) + somtomin
        Rng(icell)%passive_soil_nitrogen = Rng(icell)%passive_soil_nitrogen - tndec
        Rng(icell)%fast_soil_nitrogen(ilayer) = Rng(icell)%fast_soil_nitrogen(ilayer) + demnsom1                 ! Note two summations for fast soil nitrogen.  Odd, but ok.

        Rng(icell)%grossmin(ilayer) = Rng(icell)%grossmin(ilayer) + tndec
        Rng(icell)%tnetmin(ilayer) = Rng(icell)%tnetmin(ilayer) + somtomin - smintosom
        if (Rng(icell)%tnetmin(ilayer) .lt. 0.0) Rng(icell)%tnetmin(ilayer) = 0.0
        Rng(icell)%tminup(ilayer) = Rng(icell)%tminup(ilayer) + smintosom

      end if
    end if                  ! End if layer = 2

    ! Invertebrate decomposition or herbivory of structural litter.   C is respired, N is recycled
    decinv(ilayer) = Parms(iunit)%decomp_rate_structural_litter_inverts(ilayer) * &
                     Rng(icell)%temp_effect_on_decomp * decodt
    decc = Rng(icell)%litter_structural_carbon(ilayer) * decinv(ilayer)
    decc = min(decc, Rng(icell)%litter_structural_carbon(ilayer))
    if (Rng(icell)%litter_structural_carbon(ilayer) .gt. 0.0) then
      decn = decc * Rng(icell)%litter_structural_nitrogen(ilayer) / Rng(icell)%litter_structural_carbon(ilayer)
    else
      decn = 0.0
    end if
    Rng(icell)%litter_structural_nitrogen(ilayer) = Rng(icell)%litter_structural_nitrogen(ilayer) - decn
    Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) + decn
    Rng(icell)%tnetmin(ilayer) = Rng(icell)%tnetmin(ilayer) + decn
    Rng(icell)%litter_structural_carbon(ilayer) = Rng(icell)%litter_structural_carbon(ilayer) - decc
    ! An entry in Savanna storing a temporary accumulator appears never to be used, and was skipped.

    ! Fine branch  (these are ordered, first do fine then do coarse)
    call Track_Lignin(icell, SURFACE_INDEX, 1, Rng(icell)%dead_total_fine_branch_carbon, &
                     Rng(icell)%dead_total_fine_branch_nitrogen, Rng(icell)%all_effects_on_decomp, decodt, 0.25, grmin, immobil)
    Rng(icell)%grossmin(ilayer) = Rng(icell)%grossmin(ilayer) + grmin
    Rng(icell)%tnetmin(ilayer) = Rng(icell)%tnetmin(ilayer) - immobil
    if (Rng(icell)%tnetmin(ilayer) .lt. 0.0) Rng(icell)%tnetmin(ilayer) = 0.0
    Rng(icell)%tminup(ilayer) = Rng(icell)%tminup(ilayer) + immobil
    ! Coarse branch
    call Track_Lignin(icell, SURFACE_INDEX, 2, Rng(icell)%dead_total_coarse_branch_carbon, &
                     Rng(icell)%dead_coarse_total_branch_nitrogen, Rng(icell)%all_effects_on_decomp, decodt, 0.25, grmin, immobil)
    Rng(icell)%grossmin(ilayer) = Rng(icell)%grossmin(ilayer) + grmin
    Rng(icell)%tnetmin(ilayer) = Rng(icell)%tnetmin(ilayer) - immobil
    if (Rng(icell)%tnetmin(ilayer) .lt. 0.0) Rng(icell)%tnetmin(ilayer) = 0.0
    Rng(icell)%tminup(ilayer) = Rng(icell)%tminup(ilayer) + immobil
    ! Coarse root
    call Track_Lignin(icell, SOIL_INDEX, 3, Rng(icell)%dead_total_coarse_root_carbon, &
                     Rng(icell)%dead_total_coarse_root_nitrogen, Rng(icell)%all_effects_on_decomp, decodt, 0.25, grmin, immobil)
    Rng(icell)%grossmin(ilayer) = Rng(icell)%grossmin(ilayer) + grmin
    Rng(icell)%tnetmin(ilayer) = Rng(icell)%tnetmin(ilayer) - immobil
    if (Rng(icell)%tnetmin(ilayer) .lt. 0.0) Rng(icell)%tnetmin(ilayer) = 0.0
    Rng(icell)%tminup(ilayer) = Rng(icell)%tminup(ilayer) + immobil

    ! N gains and losses from system
    if (ilayer .eq. 1) then
      ! Fraction could include a clay component, being a function of soil texture
      Rng(icell)%volitn(ilayer) = Parms(iunit)%fraction_gross_n_mineral_volatized * Rng(icell)%grossmin(ilayer)
      volex = Parms(iunit)%rate_volatization_mineral_n * Rng(icell)%mineral_nitrogen(ilayer) * decodt
      Rng(icell)%volitn(ilayer) = Rng(icell)%volitn(ilayer) + volex
      Rng(icell)%volitn(ilayer) = min(Rng(icell)%volitn(ilayer), Rng(icell)%mineral_nitrogen(ilayer))
      Rng(icell)%volitn(ilayer) = max(0., Rng(icell)%volitn(ilayer))
      Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) - Rng(icell)%volitn(ilayer)
    end if
  end do

  ! Fixation and deposition
  ! Using the calculation for base_n_deposition done with Century, already completed and more clear than in Savanna.  That is the one that used EPNFA in Century, and calculated here in Misc_Materials, once a year.
  ! REMOVEING the symbiotic parameter in the LAND UNITS set.  Symbiotic fixation is handled elsewhere.
  fixa = Parms(iunit)%precip_n_deposition(1) * ( 12. / 365. ) + ( Parms(iunit)%precip_n_deposition(2) * &
     Globe(Rng(icell)%x, Rng(icell)%y)%precip )
  fixa = max(0., fixa)

  ! tbio is used to sum total green biomass.  "total_aground_live_biomass" will store that variable, so using that.
  tbio = Rng(icell)%total_aground_live_biomass * 0.4
  biof = fxmca + (fxmcb * tbio)        ! * N_DECOMP_LOOPS
  fxbiom = 1. - biof
  fxbiom = min(1., fxbiom)
  temper = Globe(Rng(icell)%x, Rng(icell)%y)%temperature_average
  if (fxbiom .lt. 0. .or. temper .lt. 7.5) then
    fwdfx = 0.
  else
    fwdfx = fxbiom
  end if
  fixs = fxmxs * fwdfx

  ! Total fixation
  Rng(icell)%fixnit = fixa + fixs

  ! Convert to g/m2 for the facet cover, and add to the soil
  ! Rng(icell)%fixnit = Rng(icell)%fixnit
  Rng(icell)%mineral_nitrogen(SOIL_INDEX) = Rng(icell)%mineral_nitrogen(SOIL_INDEX) + Rng(icell)%fixnit
  ! Incorporate runoff from the cell.  By rights, this should have a minor effect on our large cells.
  ! Note different units on runoff here than in SAVANNA.  CM in Century, most likely, MM in Savanna.
  Rng(icell)%runoffn = Parms(iunit)%precip_n_deposition(2) * Rng(icell)%runoff
  Rng(icell)%mineral_nitrogen(1) = Rng(icell)%mineral_nitrogen(1) - Rng(icell)%runoffn

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'DECOMPN')

end subroutine



subroutine Partition_Litter (icell, ilayer, dthc, dthn, lignin_conc, calling)
!**** Partition the litter into structural and metabolic components, based on the lignin C to N ratio (PARTLITR in Savanna, kept in DECOMP.F)
!****
!**** R. Boone   Last modified: April 9, 2011
  use Parameter_Vars
  use Structures
  implicit none

  real dthc, dthn, lignin_conc
  real new_lignin, frn_conc, rlnres, frmet, frstruc, flow_struc, flow_metab, old_lignin, flow_strucn, flow_metabn
  character(20) calling

  integer icell, ilayer, ifacet                               ! Ilayer is the surface (1) or soil (2)

  if (dthc .lt. 0.0001) return                            ! Nothing to partition

  ! A section of code is turned off in this version using dirabs=0, and commented out in Sav5b4, and so not included here. No direct absorbtion included here.

  ! N content, using a biomass basis and 2.5 conversion
  frn_conc = dthn / (dthc*2.5)                                                           ! Greater than 0 checked above.
  rlnres = lignin_conc / (frn_conc+1.e-6)                                                ! Addition of 1.e-6 prevents 0 division.
  rlnres = max(0.,rlnres)
  frmet = 0.85 - 0.013 * rlnres                       ! The values are parameters taken from Century, spl()
  frmet = max(0.,frmet)
  frstruc = 1. - frmet
  ! Ensure the structural fraction is not greater than the lignin fraction
  if (lignin_conc .gt. frstruc) then
    frstruc = lignin_conc
    frmet = 1. - frstruc
  endif

  ! Put a minimum fraction of materials to metabolic
  if (frmet .lt. 0.2) then
    frmet = 0.2
    frstruc = 0.8
  endif

  flow_struc = frstruc * dthc
  flow_metab = frmet * dthc

  ! Flows to metaboloic and structural components
  do ifacet = 1, FACETS
    ! Lignin content of the structural litter.  Lignin_structural_residue is a fraction.
    old_lignin = Rng(icell)%plant_lignin_fraction(ifacet, ilayer) * Rng(icell)%litter_structural_carbon(ilayer)
    new_lignin = lignin_conc * dthc

    if ((Rng(icell)%litter_structural_carbon(ilayer) + flow_struc) .gt. 0.0001) then
      Rng(icell)%plant_lignin_fraction(ifacet,ilayer) = &
                                     ( old_lignin + new_lignin ) / ( Rng(icell)%litter_structural_carbon(ilayer) + flow_struc )
    else
      Rng(icell)%plant_lignin_fraction(ifacet,ilayer) = 0.05                                                                ! Assigned a typical lignin concentration.  Should not occur, but to prevent errors ...  Adjust at will.
    end if
    Rng(icell)%plant_lignin_fraction(ifacet, ilayer) = max(0.02, Rng(icell)%plant_lignin_fraction(ifacet, ilayer) )
    Rng(icell)%plant_lignin_fraction(ifacet, ilayer) = min(0.50, Rng(icell)%plant_lignin_fraction(ifacet, ilayer) )
  end do
  Rng(icell)%litter_structural_carbon(ilayer) = Rng(icell)%litter_structural_carbon(ilayer) + flow_struc
  Rng(icell)%litter_metabolic_carbon(ilayer) = Rng(icell)%litter_metabolic_carbon(ilayer) + flow_metab

  flow_strucn = flow_struc / 200.                                                    ! 200 is C:N ratio in the "elements" structural litter C:E ratio array RCESTR() in Century and Savanna.  CENTURY uses 150 for nitrogen, and a note in Savanna says "200 is used everywhere" so for simplicity ...
  if (flow_strucn .gt. dthn) flow_strucn = dthn

  Rng(icell)%litter_structural_nitrogen(ilayer) = Rng(icell)%litter_structural_nitrogen(ilayer) + flow_strucn

  flow_metabn = dthn - flow_strucn
  flow_metabn = max(0.,flow_metabn)
  Rng(icell)%litter_metabolic_nitrogen(ilayer) = Rng(icell)%litter_metabolic_nitrogen(ilayer) + flow_metabn

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'PARTLIT')

end subroutine


subroutine Effects_on_Decomposition(icell)
!**** Calculate effects on decomposition, including temperature, water, and anerobic effects from precipitation.
!****
!**** R. Boone   Last modified: July 14, 2014
!****            The differences in decomposition are the source of very different soil total carbon seen,
!****            and it appears to be because of these coefficients.  Error in calculating minimum of temperature effect.
!****            Replacing runoff estimate, which is not well tracked, with drain estimate, which is a variable in Century and now in G-Range.  That variable is a general coefficient, DRAIN=1 for sandy soils, DRAIN=0 for clay soils.
  use Parameter_Vars
  use Structures
  implicit none

  real normalizer, a, b, c, d, x, catanf, xh2o, newrat, slope
  integer icell, iunit

  ! CYCLE Contains calculations used in decomposition, after the H2OLoss call, and after plant production has been calculated.
  ! I will include those here.  STEMP deals with temperature, and is in productivity.  So here we start with TFUNC, line 151.

  ! Not using TCALC function call here, just incorporating the function here.
  ! An exponential function is used in CENTURY, dependent upon 4 parameters, stored in an array.
  iunit = Rng(icell)%range_type
  x = 30.0
  a = Parms(iunit)%temperature_effect_decomposition(1)
  b = Parms(iunit)%temperature_effect_decomposition(2)
  c = Parms(iunit)%temperature_effect_decomposition(3)
  d = Parms(iunit)%temperature_effect_decomposition(4)
  normalizer = b + (c / PI) * atan(PI * d * (x - a))                                 ! Note standardized to 30.
  x = Rng(icell)%soil_surface_temperature
  catanf = b + (c / PI) * atan(PI * d * (x - a))
  Rng(icell)%temp_effect_on_decomp = max( 0.01, ( catanf / normalizer ) )

  ! Two methods of calculating water effect on decomposition in CENTURY.  I am incorporating the second, as it appears to perform better in this setting.
  if (Rng(icell)%ratio_water_pet .gt. 9.) then
    Rng(icell)%water_effect_on_decomp = 1.0
  else
    Rng(icell)%water_effect_on_decomp = 1./(1.+30.0*exp(-8.5*Rng(icell)%ratio_water_pet))         ! 1+ demonator prevents division by 0.
  endif

  ! Anerobic effects on decomposition  (ANEROB.F)
  Rng(icell)%anerobic_effect_on_decomp = 1.
  if (Rng(icell)%ratio_water_pet .gt. Parms(iunit)%anerobic_effect_decomposition(1)) then
    xh2o = ( Rng(icell)%ratio_water_pet - Parms(iunit)%anerobic_effect_decomposition(1) ) * &
             Rng(icell)%pot_evap * ( 1.0 - Parms(iunit)%drainage_affecting_anaerobic_decomp )

    if (xh2o .gt. 0.) then
      if (Rng(icell)%pot_evap .gt. 0.0) then                                               ! Avoiding division error
        newrat = Parms(iunit)%anerobic_effect_decomposition(1) + (xh2o / Rng(icell)%pot_evap)
      else
        newrat = Parms(iunit)%anerobic_effect_decomposition(1)
      end if
      slope = (1.0 - Parms(iunit)%anerobic_effect_decomposition(3)) / &
              (Parms(iunit)%anerobic_effect_decomposition(1) - Parms(iunit)%anerobic_effect_decomposition(2))          ! Parameters, so unlikely to yield division by 0.
      Rng(icell)%anerobic_effect_on_decomp = 1.0 + slope * (newrat - Parms(iunit)%anerobic_effect_decomposition(1))
    endif
    if (Rng(icell)%anerobic_effect_on_decomp .lt. Parms(iunit)%anerobic_effect_decomposition(3)) then
      Rng(icell)%anerobic_effect_on_decomp = Parms(iunit)%anerobic_effect_decomposition(3)
    endif
  endif

  if (Rng(icell)%anerobic_effect_on_decomp .lt. 0.0) Rng(icell)%anerobic_effect_on_decomp = 0.0
  if (Rng(icell)%anerobic_effect_on_decomp .gt. 1.0) Rng(icell)%anerobic_effect_on_decomp = 1.0
  if (Rng(icell)%temp_effect_on_decomp .lt. 0.0) Rng(icell)%temp_effect_on_decomp = 0.0
  if (Rng(icell)%temp_effect_on_decomp .ge. 1.0) Rng(icell)%temp_effect_on_decomp = 1.0
  if (Rng(icell)%water_effect_on_decomp .lt. 0.0) Rng(icell)%water_effect_on_decomp = 0.0
  if (Rng(icell)%water_effect_on_decomp .ge. 1.0) Rng(icell)%water_effect_on_decomp = 1.0

  ! Combining effects of temperature and moisture  (CYCLE stores a monthly value ... I don't think I want or need this.  Only monthly cycles are represented here, and the past is the past, and the future unknown.
  Rng(icell)%all_effects_on_decomp =  Rng(icell)%temp_effect_on_decomp * Rng(icell)%water_effect_on_decomp * &
                                      Rng(icell)%anerobic_effect_on_decomp
  if (Rng(icell)%all_effects_on_decomp .lt. 0.0) Rng(icell)%all_effects_on_decomp = 0.0                ! DEFAC in Savanna DECOMP.F
  if (Rng(icell)%all_effects_on_decomp .gt. 1.0) Rng(icell)%all_effects_on_decomp = 1.0                ! DEFAC in Savanna DECOMP.F

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'EFF_LIG')

end subroutine



subroutine Track_Lignin(icell, ilayer, ic, woodc, woodn, defac, decodt, frlig, tndec, smintosom)
!**** A routine that tracks lignin, called DECLIG in Savanna.
!****
!**** R. Boone   Last modified: April 6, 2020
  use Parameter_Vars
  use Structures
  implicit none

  real defac, decodt, tndec, smintosom, woodc, woodn, frlig
  real ps1co2(2), rsplig, eflig, decrt, tcdec, tosom2gross, co2los2, tosom2net, tosom1net, tosom1gross, co2los1
  real demnsom1, demnsom2, demnt, rces1, rces2, rces3, fr, demn1, demn2
  integer icell, ilayer, ic
  integer iunit

  ps1co2(1) = 0.45
  ps1co2(2) = 0.55
  rsplig = 0.3
  rces1 = 13.0               ! Default values for rces in Savanna
  rces2 = 18.0
  rces3 = 8.0

  tndec = 0.0
  smintosom = 0.0

  iunit = Rng(icell)%range_type

  eflig = exp(-3. * frlig)
  select case (ic)
    case (0)
      decrt = Parms(iunit)%decomp_rate_structural_litter(ilayer) * defac * eflig * decodt
    case (1)
      decrt = Parms(iunit)%decomp_rate_fine_branch * defac * eflig * decodt
    case (2)
      decrt = Parms(iunit)%decomp_rate_coarse_branch * defac * eflig * decodt
    case (3)
      decrt = Parms(iunit)%decomp_rate_coarse_root * defac * eflig * decodt
  end select

  tcdec = woodc * decrt

  if (tcdec .gt. 0.) then
    tosom2gross = tcdec * frlig
    co2los2 = tosom2gross * rsplig
    tosom2net = tosom2gross - co2los2

    tosom1gross = tcdec - tosom2gross
    co2los1 = tosom1gross * ps1co2(ilayer)
    tosom1net = tosom1gross - co2los1

    if (woodc .gt. 0.0) then
      tndec = tcdec * woodn / woodc
    else
      tndec = 0.0
    end if

    demnsom1 = tosom1net / rces1
    demnsom2 = tosom2net / rces2
    demnt = demnsom1 + demnsom2

    if (demnt .gt. tndec) then
      smintosom = demnt - tndec
    else
      smintosom = 0.
    end if

    if (smintosom .gt. 0.) then
      if (smintosom .gt. Rng(icell)%mineral_nitrogen(ilayer)) then                             ! NOTE: Possible error, mineral nitrogen is soil layers, ilayer is surface versus soil.  Won't break program, but logic may not work.
        Rng(icell)%mineral_nitrogen(ilayer) = max(Rng(icell)%mineral_nitrogen(ilayer), 0.)
        fr = Rng(icell)%mineral_nitrogen(ilayer) / smintosom
        smintosom = Rng(icell)%mineral_nitrogen(ilayer)
        tcdec = tcdec * fr
        tndec = tndec * fr
        tosom1net = tosom1net * fr
        tosom2net = tosom2net * fr
      end if
    end if

    Rng(icell)%fast_soil_carbon(ilayer) = Rng(icell)%fast_soil_carbon(ilayer) + tosom1net
    Rng(icell)%intermediate_soil_carbon = Rng(icell)%intermediate_soil_carbon + tosom2net
    woodc = woodc - tcdec

    demn1 = tosom1net / rces1
    demn2 = tosom2net / rces2
    Rng(icell)%fast_soil_nitrogen(ilayer) = Rng(icell)%fast_soil_nitrogen(ilayer) + demn1
    Rng(icell)%intermediate_soil_nitrogen = Rng(icell)%intermediate_soil_nitrogen + demn2
    woodn = woodn - tndec
    Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) - smintosom
  end if

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'TRK_LIG')

end subroutine


subroutine Nitrogen_Losses (icell)
!**** Nitrogen volatization and leaching routine, building on Century's version.
!**** Century's version was simplified to only deal with nitrogen.
!****
!**** R. Boone   Last modified: October 5, 2013
  use Parameter_Vars
  use Structures
  implicit none

  real base, linten, strm, minlch, texture_effect
  real volex
  integer icell, ilayer, nxt, iunit

  iunit = Rng(icell)%range_type

  !-----------------------------
  ! LEACHING FOLLOWS
  minlch = 18.0                      ! Default value for Century, a parameter there.  Hardwired here.  Revisit if needed.

  ! Century includes stream flow.  These cells will be large, and so the importance of stream flow should be limited.
  ! Beyond that, it appears to be an accumulator just for output and inspection, so skipped here.
  base = 0.0

  do ilayer = 1, SOIL_LAYERS
    Rng(icell)%n_leached(ilayer) = 0.
    strm = 0.
    nxt = ilayer + 1                                                      ! Call for ilayer = 4 will cause the system to fail, given that nxt would equal 5, which is not defined.

    if ((Rng(icell)%amov(ilayer) .gt. 0.0) .and. (Rng(icell)%mineral_nitrogen(ilayer) .gt. 0.0)) then
      linten = min(1.0 - (minlch - Rng(icell)%amov(ilayer)) / minlch, 1.0)
      linten = max(linten, 0.0)
      texture_effect = 0.2 + 0.7 * Rng(icell)%sand(ilayer)                ! Using Savanna's version, rather than Century
      Rng(icell)%n_leached(ilayer) = texture_effect * Rng(icell)%mineral_nitrogen(ilayer) * linten

      ! If at the bottom of the stack of layers, compute storm flow
      if (ilayer .eq. SOIL_LAYERS) then
        strm = Rng(icell)%n_leached(ilayer) * Rng(icell)%storm_flow
      end if
      Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) - Rng(icell)%n_leached(ilayer)
      if (nxt .lt. SOIL_LAYERS) then
        Rng(icell)%mineral_nitrogen(nxt) = Rng(icell)%mineral_nitrogen(nxt) + (Rng(icell)%n_leached(ilayer) - strm)
      else
        base = Rng(icell)%mineral_nitrogen(ilayer) * Parms(Rng(icell)%range_type)%base_flow_fraction
        Rng(icell)%mineral_nitrogen(ilayer) = Rng(icell)%mineral_nitrogen(ilayer) - base
      end if
    end if
  end do

  ! Century uses extra soil layers (up to 10), with layer + 1 storing base flow.   I only have the four layers.
  ! Leaching is moving down.  Must store that using a unique approach, relfected in the else section immediately above.
  ! Streamflow not simulated.
  ! END OF LEACHING
  !---------------------------------

  !---------------------------------
  ! Given the relatedness and brevity, adding VOLATILIZATION here
  if (Rng(icell)%mineral_nitrogen(SURFACE_INDEX) .gt. 0.0) then
    ! Annual fraction of mineral n volatized
    volex = ( Parms(iunit)%annual_fraction_volatilized_n * (1./12.)) * Rng(icell)%mineral_nitrogen(SURFACE_INDEX)
    Rng(icell)%mineral_nitrogen(SURFACE_INDEX) = Rng(icell)%mineral_nitrogen(SURFACE_INDEX) - volex
    Rng(icell)%nitrogen_source_sink = Rng(icell)%nitrogen_source_sink + volex
  end if

  if (check_nan_flag .eq. .TRUE.) call check_for_nan (icell, 'NIT_LOSS')

end subroutine
