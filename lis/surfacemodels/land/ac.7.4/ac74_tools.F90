!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ac74_tools
  !BOP
  !
  ! !MODULE: ac74_tools
  !
  ! !DESCRIPTION:
  !  Tools of the LIS AquaCrop driver:
  !  \begin{description}
  !  \item[AC74\_ETo\_calc] reference evapotranspiration from the forcing
  !  \item[AC74\_plan\_run] plans the next run of a tile (fallow, then crop)
  !  \item[AC74\_plant\_search] planting criteria, day by day
  !  \item[AC74\_set\_run\_input] passes the planned run and the external
  !    inputs to AquaCrop, before the run is initialized
  !  \end{description}
  !
  !EOP

contains
  !BOP
  !
  ! !ROUTINE: AC74_ETo_calc
  ! \label{AC74_ETo_calc}
  !
  ! !REVISION HISTORY:
  !  04 NOV 2024, Louise Busschaert; initial implementation
  !

  ! !INTERFACE:
  subroutine AC74_ETo_calc(P, Tmax, Tmin, Tdew, ws, Rs, z, lat, eto)

    ! !USES:
    use LIS_constantsMod, only: LIS_CONST_PI, LIS_CONST_TKFRZ
    use LIS_coreMod, only: LIS_rc
    use LIS_tbotAdjustMod, only: LIS_tbotTimeUtil

    !
    ! !DESCRIPTION:
    !
    !  This routine computes the reference evapotranspiration (ETo)
    !  with the Penman-Monteith equation, following the guidelines
    !  of the FAO Irrigation and Drainage Paper No. 56.
    !
    !
    !EOP

    implicit none

    real, intent(in)      :: P ! kPa
    real, intent(in)      :: Tmax, Tmin, Tdew
    real, intent(in)      :: ws ! wind speed
    real, intent(in)      :: Rs ! radiation
    real, intent(in)      :: z ! elevation
    real, intent(in)      :: lat
    real, intent(inout)   :: eto ! returns eto [mm day-1]

    real                  :: Tmean
    real                  :: gamma, ea, es, slope
    real                  :: phi, dr, delta, omega, Ra, Rso
    real                  :: Rns, Rnl, Rn
    real                  :: julian_in, ratio

    ! Tmean (degC)
    Tmean = (Tmin + Tmax)/2.

    ! Psychrometric constant [kPa degC-1]
    gamma = 0.664742 * 0.001 * P

    ! Mean saturation vapour pressure (es) [kPa]
    es = (0.6108 * EXP((17.27 * Tmin) / (237.3 + Tmin)) &
         + 0.6108 * EXP((17.27 * Tmax) / (237.3 + Tmax))) / 2
    ! Actual vapor pressure (ea) [kPa]
    ea = 0.6108 * EXP((17.27 * Tdew) / (237.3 + Tdew))

    ! Slope of saturation vapour pressure curve [kPa degC-1]
    slope = (4098 * (0.6108 * EXP((17.27 * Tmean)/(Tmean + 237.3)))) &
         / (Tmean + 237.3)**2

    ! Extraterrestrial radiation (Ra) [MJ m-2 day-1]
    call LIS_tbotTimeUtil(julian_in,LIS_rc%yr) ! First get day of year
    phi = lat * LIS_CONST_PI/180 ! get latitude in rad
    dr = 1 + 0.033 * COS(2*LIS_CONST_PI/365 * julian_in) ! Inverse relative distance Earth-Sun
    delta = 0.409 * SIN(2*LIS_CONST_PI/365 * julian_in - 1.39) ! Solar declination
    omega = ACOS(-TAN(phi)*TAN(delta))
    Ra = (24 * 60/LIS_CONST_PI) * 0.082 * dr * (omega*SIN(phi)*SIN(delta) &
         + COS(phi)*COS(delta)*SIN(omega))
    if(Ra.le.0) then ! return 0 ETo for unrealisic numbers (e.g. polar night)
       eto = 0
       return
    endif

    ! Net radiation at crop surface (Rn) [MJ m-2 day-1]
    Rso = (0.75 + 2E-5 * z) * Ra
    Rns = (1 - 0.23) * Rs ! fixed value of 0.23 for albedo
    ratio = Rs/Rso
    if (ratio.gt.1) then
       ratio=1
    endif
    Rnl = 4.903E-9 * (((Tmax + LIS_CONST_TKFRZ)**4 + (Tmin + LIS_CONST_TKFRZ)**4) / 2.) &
         * (0.34 - 0.14 * SQRT(ea)) &
         * (1.35 * ratio - 0.35) ! Net longwave radiation [MJ m-2 day-1]
    Rn = Rns - Rnl

    ! Penman-Monteith --> reference evapotranspiration [mm day-1]
    eto = (0.408 * slope * Rn + (gamma * (900 / (Tmean + LIS_CONST_TKFRZ)) * ws * (es - ea))) &
         / (slope + gamma * (1 + 0.34 * ws))
    if(eto.le.0) then
       eto = 0 ! avoid negative values
    endif

  end subroutine AC74_ETo_calc

  !BOP
  !
  ! !ROUTINE: AC74_plan_run
  ! \label{AC74_plan_run}
  !
  ! !REVISION HISTORY:
  !  19 OCT 2026, Louise Busschaert; initial implementation
  !

  ! !INTERFACE:
  subroutine AC74_plan_run(n, t, FromDayNr)

    ! !USES:
    use AC74_lsmMod, only: AC74_struc
    use ac74_global, only: DetermineDate, DetermineDayNr

    !
    ! !DESCRIPTION:
    !
    !  Plans the next run of a tile, starting on FromDayNr: fallow until the
    !  planting day of the main season, then the crop until the end of its
    !  cycle, at the latest MaxCropLength days after planting.
    !  With a planting criterion, the planting day is the last day of the
    !  window, until the criterion is met (AC74_plant_search).
    !
    !EOP

    implicit none

    integer, intent(in) :: n, t
    integer, intent(in) :: FromDayNr

    integer :: Dayi, Monthi, Yeari, PlantDayNr

    ! first planting day on or after FromDayNr
    call DetermineDate(FromDayNr, Dayi, Monthi, Yeari)
    call DetermineDayNr(AC74_struc(n)%Crop_AnnualStartDay, &
         AC74_struc(n)%Crop_AnnualStartMonth, Yeari, PlantDayNr)
    if (PlantDayNr < FromDayNr) then
       call DetermineDayNr(AC74_struc(n)%Crop_AnnualStartDay, &
            AC74_struc(n)%Crop_AnnualStartMonth, Yeari+1, PlantDayNr)
    endif

    AC74_struc(n)%ac74(t)%PlantSearch = &
         AC74_struc(n)%Temp_crit .or. AC74_struc(n)%Rainfall_crit
    AC74_struc(n)%ac74(t)%WindowStart = PlantDayNr
    AC74_struc(n)%ac74(t)%TempCritDays = 0
    AC74_struc(n)%ac74(t)%TempCritOcc = 0
    AC74_struc(n)%ac74(t)%TempCritMet = .not. AC74_struc(n)%Temp_crit
    AC74_struc(n)%ac74(t)%RainCritN = 0
    AC74_struc(n)%ac74(t)%RainCritOcc = 0
    AC74_struc(n)%ac74(t)%RainCritMet = .not. AC74_struc(n)%Rainfall_crit
    if (AC74_struc(n)%ac74(t)%PlantSearch) then
       PlantDayNr = PlantDayNr + AC74_struc(n)%crit_window - 1
    endif

    AC74_struc(n)%ac74(t)%RunFromDayNr = FromDayNr
    AC74_struc(n)%ac74(t)%RunCropDay1 = PlantDayNr
    AC74_struc(n)%ac74(t)%RunLastDayNr = PlantDayNr + AC74_struc(n)%MaxCropLength - 1
  end subroutine AC74_plan_run

  !BOP
  !
  ! !ROUTINE: AC74_plant_search
  ! \label{AC74_plant_search}
  !
  ! !REVISION HISTORY:
  !  19 OCT 2026, Louise Busschaert; initial implementation
  !

  ! !INTERFACE:
  subroutine AC74_plant_search(n, t, DayNr, Tmin, Rain)

    ! !USES:
    use AC74_lsmMod, only: AC74_struc

    !
    ! !DESCRIPTION:
    !
    !  Updates the planting criteria of a tile with the forcing of the day
    !  DayNr, within the planting window. The crop is planted the first day
    !  all criteria are met: the tile then starts a new run on that day
    !  (InitializeRun), or on the last day of the window if not met before.
    !  Temperature: Temp_crit_occurrence times Temp_crit_days consecutive
    !  days with Tmin >= Temp_crit_tmin. Rainfall: Rainfall_crit_occurrence
    !  times more than Rainfall_crit_amount within Rainfall_crit_days.
    !  Occurrences do not overlap.
    !
    !EOP

    implicit none

    integer, intent(in) :: n, t
    integer, intent(in) :: DayNr
    real, intent(in)    :: Tmin, Rain

    integer :: nd

    if (.not. AC74_struc(n)%ac74(t)%PlantSearch) return
    if (DayNr < AC74_struc(n)%ac74(t)%WindowStart) return
    ! last day of the window: planted anyway
    if (DayNr >= AC74_struc(n)%ac74(t)%RunCropDay1) then
       AC74_struc(n)%ac74(t)%PlantSearch = .false.
       return
    endif

    ! temperature criterion
    if (.not. AC74_struc(n)%ac74(t)%TempCritMet) then
       if (Tmin >= AC74_struc(n)%Temp_crit_tmin) then
          AC74_struc(n)%ac74(t)%TempCritDays = AC74_struc(n)%ac74(t)%TempCritDays + 1
          if (AC74_struc(n)%ac74(t)%TempCritDays == AC74_struc(n)%Temp_crit_days) then
             AC74_struc(n)%ac74(t)%TempCritOcc = AC74_struc(n)%ac74(t)%TempCritOcc + 1
             AC74_struc(n)%ac74(t)%TempCritDays = 0
          endif
       else
          AC74_struc(n)%ac74(t)%TempCritDays = 0
       endif
       AC74_struc(n)%ac74(t)%TempCritMet = &
            (AC74_struc(n)%ac74(t)%TempCritOcc >= AC74_struc(n)%Temp_crit_occurrence)
    endif

    ! rainfall criterion: rain of the last days, since the last occurrence
    if (.not. AC74_struc(n)%ac74(t)%RainCritMet) then
       nd = AC74_struc(n)%Rainfall_crit_days
       if (AC74_struc(n)%ac74(t)%RainCritN < nd) then
          AC74_struc(n)%ac74(t)%RainCritN = AC74_struc(n)%ac74(t)%RainCritN + 1
       else
          AC74_struc(n)%ac74(t)%RainCritBuf(1:nd-1) = AC74_struc(n)%ac74(t)%RainCritBuf(2:nd)
       endif
       AC74_struc(n)%ac74(t)%RainCritBuf(AC74_struc(n)%ac74(t)%RainCritN) = Rain
       if (sum(AC74_struc(n)%ac74(t)%RainCritBuf(1:AC74_struc(n)%ac74(t)%RainCritN)) &
            > AC74_struc(n)%Rainfall_crit_amount) then
          AC74_struc(n)%ac74(t)%RainCritOcc = AC74_struc(n)%ac74(t)%RainCritOcc + 1
          AC74_struc(n)%ac74(t)%RainCritN = 0
       endif
       AC74_struc(n)%ac74(t)%RainCritMet = &
            (AC74_struc(n)%ac74(t)%RainCritOcc >= AC74_struc(n)%Rainfall_crit_occurrence)
    endif

    ! planted today: new run from today, without search
    if (AC74_struc(n)%ac74(t)%TempCritMet .and. AC74_struc(n)%ac74(t)%RainCritMet) then
       AC74_struc(n)%ac74(t)%PlantSearch = .false.
       AC74_struc(n)%ac74(t)%RunFromDayNr = DayNr
       AC74_struc(n)%ac74(t)%RunCropDay1 = DayNr
       AC74_struc(n)%ac74(t)%RunLastDayNr = DayNr + AC74_struc(n)%MaxCropLength - 1
       AC74_struc(n)%ac74(t)%irun = AC74_struc(n)%ac74(t)%irun + 1
       AC74_struc(n)%ac74(t)%InitializeRun = 1
    endif
  end subroutine AC74_plant_search

  !BOP
  !
  ! !ROUTINE: AC74_set_run_input
  ! \label{AC74_set_run_input}
  !
  ! !REVISION HISTORY:
  !  19 OCT 2026, Louise Busschaert; initial implementation
  !

  ! !INTERFACE:
  subroutine AC74_set_run_input(n, t)

    ! !USES:
    use AC74_lsmMod, only: AC74_struc
    use ac74_global, only: SetCrop_GDDaysToGermination, SetCrop_GDDaysToHarvest, &
         SetCrop_GDDaysToSenescence, SetCrop_GDDaysToMaxRooting, &
         SetCrop_GDDCGC, SetCrop_GDDCDC, SetManagement_Brel
    use ac74_project_input, only: set_project_input

    !
    ! !DESCRIPTION:
    !
    !  Writes the planned run of a tile in the (single) AquaCrop project run,
    !  before the run is initialized. Also sets the inputs that are driven
    !  externally (LIS param file): the crop GDD parameters, used where the
    !  crop file has 'External', and the relative biomass (soil fertility).
    !
    !EOP

    implicit none

    integer, intent(in) :: n, t

    if (AC74_struc(n)%GDDfromLDT) then
       call SetCrop_GDDaysToGermination(AC74_struc(n)%ac74(t)%GDDaysToGermination)
       call SetCrop_GDDaysToHarvest(AC74_struc(n)%ac74(t)%GDDaysToHarvest)
       call SetCrop_GDDaysToSenescence(AC74_struc(n)%ac74(t)%GDDaysToSenescence)
       call SetCrop_GDDaysToMaxRooting(AC74_struc(n)%ac74(t)%GDDaysToMaxRooting)
       call SetCrop_GDDCGC(AC74_struc(n)%ac74(t)%GDDCGC)
       call SetCrop_GDDCDC(AC74_struc(n)%ac74(t)%GDDCDC)
    endif
    ! Brel < 0: not supplied, the soil fertility of the management file is used
    call SetManagement_Brel(AC74_struc(n)%ac74(t)%Brel)

    call set_project_input(1, 'Simulation_DayNr1', AC74_struc(n)%ac74(t)%RunFromDayNr)
    call set_project_input(1, 'Simulation_DayNrN', AC74_struc(n)%ac74(t)%RunLastDayNr)
    call set_project_input(1, 'Crop_Day1', AC74_struc(n)%ac74(t)%RunCropDay1)
    call set_project_input(1, 'Crop_LastDayNr', AC74_struc(n)%ac74(t)%RunLastDayNr)
    call set_project_input(1, 'Crop_Filename', &
         trim(AC74_struc(n)%ac74(t)%cropt)//'.CRO')
    ! initial conditions from the config for the first run, then kept
    if (AC74_struc(n)%ac74(t)%irun == 1) then
       call set_project_input(1, 'SWCIni_Filename', '(None)')
    else
       call set_project_input(1, 'SWCIni_Filename', 'KeepSWC')
    endif
  end subroutine AC74_set_run_input

end module ac74_tools
