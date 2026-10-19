!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
#include "LIS_plugins.h"
!BOP
!
! !ROUTINE: AC74_setup
! \label{AC74_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial
!  specification of the subroutine is defined by Sujay Kumar.
!   04 NOV 2024; Louise Busschaert, Initial implementation
!
! !INTERFACE:
subroutine AC74_setup()

  ! !USES:
  use ac74_global, only: CheckForKeepSWC, &
       DetermineDayNr, &
       FromGravelMassToGravelVolume, &
       GetCCiActual, &
       GetCCiprev, &
       GetCCiTopEarlySen, &
       GetClimRecord, &
       GetCompartment, &
       GetCompartment_theta, &
       GetCrop, &
       GetCRsalt, &
       GetCRwater, &
       GetDaySubmerged, &
       GetDrain, &
       GetEact, &
       GetECdrain, &
       GetECiAqua, &
       GetECstorage, &
       Geteffectiverain, &
       GetEndSeason, &
       GetEpot, &
       GetETo, &
       GetEToRecord, &
       GetEvapoEntireSoilSurface, &
       GetFullFileNameProgramParameters, &
       GetGenerateDepthMode, &
       GetGenerateTimeMode, &
       GetInfiltrated, &
       GetIniPercTAW, &
       GetIrriAfterSeason, &
       GetIrriBeforeSeason, &
       GetIrriECw, &
       GetIrriFirstDayNr, &
       GetIrrigation, &
       GetIrriInfoLastDay, &
       GetIrriMethod, &
       GetIrriMode, &
       GetManagement, &
       GetManagement_Cuttings, &
       GetMaxPlotNew, &
       GetMaxPlotTr, &
       GetNrCompartments, &
       GetOnset, &
       GetOut1Wabal, &
       GetOut2Crop, &
       GetOut3Prof, &
       GetOut4Salt, &
       GetOut5CompWC, &
       GetOut6CompEC, &
       GetOut7Clim, &
       GetOut8Irri, &
       GetOutDaily, &
       GetOutputAggregate, &
       GetPart1Mult, &
       GetPart2Eval, &
       GetPerennialPeriod, &
       GetPreDay, &
       GetRain, &
       GetRainRecord, &
       GetRootingDepth, &
       GetRootZoneSalt, &
       GetRootZoneWC_Actual, &
       GetRootZoneWC_FC, &
       GetRootZoneWC_Leaf, &
       GetRootZoneWC_SAT, &
       GetRootZoneWC_Sen, &
       GetRootZoneWC_Thresh, &
       GetRootZoneWC_WP, &
       GetRootZoneWC_ZtopAct, &
       GetRootZoneWC_ZtopFC, &
       GetRootZoneWC_ZtopThresh, &
       GetRootZoneWC_ZtopWP, &
       GetRunoff, &
       GetSaltInfiltr, &
       GetSimulation, &
       GetSimulation_MultipleRunConstZrx, &
       GetSimulation_MultipleRunConstZrx, &
       GetSimulation_MultipleRunWithKeepSWC, &
       GetSimulation_MultipleRunWithKeepSWC, &
       GetSimulParam, &
       GetSoil, &
       GetSoilLayer, &
       GetSumWaBal, &
       GetSurf0, &
       GetSurfaceStorage, &
       GetTact, &
       GetTemperatureRecord, &
       GetTmax, &
       GetTmin, &
       GetTotalSaltContent, &
       GetTotalWaterContent, &
       GetTpot, &
       GetZiAqua, &
       InitializeGlobalStrings, &
       IrriMode_Generate, &
       IrriMode_Manual, &
       SetClimFile, &
       SetClimRecord_DataType, &
       SetClimRecord_fromd, &
       SetClimRecord_FromDayNr, &
       SetClimRecord_fromm, &
       SetClimRecord_FromString, &
       SetClimRecord_fromy, &
       SetClimRecord_FromY, &
       SetClimRecord_NrObs, &
       SetClimRecord_tod, &
       SetClimRecord_ToDayNr, &
       SetClimRecord_tom, &
       SetClimRecord_ToString, &
       SetClimRecord_toy, &
       SetCO2Description, &
       SetEToRecord_DataType, &
       SetEToRecord_fromd, &
       SetEToRecord_FromDayNr, &
       SetEToRecord_fromm, &
       SetEToRecord_FromString, &
       SetEToRecord_fromy, &
       SetEToRecord_FromY, &
       SetEToRecord_NrObs, &
       SetEToRecord_tod, &
       SetEToRecord_ToDayNr, &
       SetEToRecord_tom, &
       SetEToRecord_ToString, &
       SetEToRecord_toy, &
       SetFullFileNameProgramParameters, &
       SetMultipleProjectDescription, &
       SetNrCompartments, &
       SetPathNameList, &
       SetPathNameOutp, &
       SetPathNameParam, &
       SetPathNameProg, &
       SetPathNameSimul, &
       SetProfFilefull, &
       SetRainRecord_DataType, &
       SetRainRecord_fromd, &
       SetRainRecord_FromDayNr, &
       SetRainRecord_fromm, &
       SetRainRecord_FromString, &
       SetRainRecord_fromy, &
       SetRainRecord_FromY, &
       SetRainRecord_NrObs, &
       SetRainRecord_tod, &
       SetRainRecord_ToDayNr, &
       SetRainRecord_tom, &
       SetRainRecord_ToString, &
       SetRainRecord_toy, &
       SetSimulation_MultipleRun, &
       SetSimulation_MultipleRunConstZrx, &
       SetSimulation_MultipleRunWithKeepSWC, &
       SetSimulation_NrRuns, &
       SetSimulation_Storage_CropString, &
       SetSoil, &
       SetSoilLayer, &
       SetTemperatureRecord_DataType, &
       SetTemperatureRecord_fromd, &
       SetTemperatureRecord_FromDayNr, &
       SetTemperatureRecord_fromm, &
       SetTemperatureRecord_FromString, &
       SetTemperatureRecord_fromy, &
       SetTemperatureRecord_FromY, &
       SetTemperatureRecord_NrObs, &
       SetTemperatureRecord_tod, &
       SetTemperatureRecord_ToDayNr, &
       SetTemperatureRecord_tom, &
       SetTemperatureRecord_ToString, &
       SetTemperatureRecord_toy, &
       SetTmax, &
       SetTmaxTnxReference12MonthsRun, &
       SetTmin, &
       SetTminTnxReference12MonthsRun, &
       SetTnxReferenceFile, &
       SetTnxReferenceYear, &
       typeproject_typeprm

  use ac74_initialsettings, only: InitializeSettings

  use ac74_kinds, only: intEnum, &
       int8

  use ac74_project_input, only: allocate_project_input, &
       set_project_input

  use ac74_run, only: tDaysZmin, &
       tGDDZmin, &
       GetDayNrCropEnded, &
       fIrri_close, &
       GetalfaHI, &
       GetalfaHIAdj, &
       GetBin, &
       GetBout, &
       GetBprevSum, &
       GetCCiActualWeedInfested, &
       GetCCoTotal, &
       GetCCxCropWeedsNoSFstress, &
       GetCCxTotal, &
       GetCCxWitheredTpotNoS, &
       GetCDCTotal, &
       GetCGCref, &
       GetCO2i, &
       GetCoeffb0, &
       GetCoeffb0Salt, &
       GetCoeffb1, &
       GetCoeffb1Salt, &
       GetCoeffb2, &
       GetCoeffb2Salt, &
       GetCutInfoRecord1, &
       GetCutInfoRecord2, &
       GetDayFraction, &
       GetDayLastCut, &
       GetDayNr1Eval, &
       GetDayNrEval, &
       GetDayNri, &
       GetFracBiomassPotSF, &
       GetfWeedNoS, &
       GetGDDayFraction, &
       GetGDDayi, &
       GetGDDCDCTotal, &
       GetGDDCGCref, &
       GetGDDTadj, &
       GetGwTable, &
       GetHItimesAT, &
       GetHItimesAT1, &
       GetHItimesAT2, &
       GetHItimesBEF, &
       GetIrriInfoRecord1, &
       GetIrriInfoRecord2, &
       GetIrriInterval, &
       GetLineNrEval, &
       GetNextSimFromDayNr, &
       GetNoYear, &
       GetNrCut, &
       GetPlotVarCrop, &
       GetPreviousBmob, &
       GetPreviousBsto, &
       GetPreviousDayNr, &
       GetPreviousStressLevel, &
       GetPreviousSum, &
       GetPreviousSumETo, &
       GetPreviousSumGDD, &
       GetScorAT1, &
       GetScorAT2, &
       GetStageCode, &
       GetStartMode, &
       GetStressLeaf, &
       GetStressSenescence, &
       GetStressSFadjNEW, &
       GetStressTot, &
       GetSumETo, &
       GetSumGDD, &
       GetSumGDDcuts, &
       GetSumGDDPrev, &
       GetSumInterval, &
       GetSumKci, &
       GetSumKcTop, &
       GetSumKcTopStress, &
       GetTadj, &
       GetTimeSenescence, &
       GetTransfer, &
       GetWaterTableInProfile, &
       GetWeedRCi, &
       GetYprevSum, &
       GetZeval, &
       GetZiprev, &
       InitializeRunPart1, &
       InitializeSimulationRunPart2, &
       SetTmax, &
       SetTmin

  use ac74_startunit, only: GetRequestDailyResults, &
       GetRequestParticularResults, &
       GetSimulation_NrRuns, &
       GetTimeAggregationResults, &
       LoadProgramParametersProjectPlugIn

  use AC74_lsmMod
  use ac74_tools, only: AC74_plan_run, AC74_set_run_input

  use LIS_constantsMod, only: LIS_CONST_TKFRZ
  use LIS_coreMod,   only: LIS_rc, LIS_surface, LIS_domain
  use LIS_metforcingMod, only: LIS_forc
#if (defined MF_MERRA2)
  use merra2_forcingMod,  only: merra2_struc
#endif
#if (defined MF_NLDAS2)
  use nldas2_forcingMod,  only: nldas2_struc
#endif
#if (defined MF_NLDAS20)
  use nldas20_forcingMod, only: nldas20_struc
#endif
  use LIS_fileIOMod, only: LIS_read_param
  use LIS_logMod, only: LIS_logunit, &
       LIS_endrun

  use module_sf_aclsm_74, only: &
       WP, SAT, FC, INFRATE, SD, CL, SI

  !
  ! !DESCRIPTION:
  !
  !  This routine is the entry point to set up the parameters
  !  required for AC74.
  !  The routines invoked are:
  !  \begin{description}
  !  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
  !    retrieves LIS parameter data from NetCDF file
  !  \item[AC74\_read\_MULTILEVEL\_param](\ref{AC74_read_MULTILEVEL_param}) \newline
  !    retrieves MULTILEVEL spatial parameter from NetCDF file
  !  \item[AC74\_read\_croptype](\ref{AC74_read_croptype}) \newline
  !    retrieves the crop type from LIS_param and table, then assigns
  !    a crop type to each tile.
  !  \end{description}
  !EOP

  implicit none

  integer           :: mtype
  integer           :: t, k, n, l
  integer           :: col, row
  real, allocatable :: placeholder(:,:)

  real              :: Z_surf, cl_tmp, si_tmp, sd_tmp, InfRate_tmp
  integer           :: REW, descr

  integer(intEnum) :: TheProjectType

  logical ::  ProgramParametersAvailable

  integer :: WindowStart
  integer :: m
  logical :: valid_forcing
  real    :: elevdiff
  real, parameter :: lapse = -0.0065 ! lapse rate (K/m), as for the met forcings
  logical :: MultipleRunWithKeepSWC_temp
  real    :: MultipleRunConstZrx_temp

  logical, external :: AC74_param_has_var
  external :: ac74_read_croptype
  external :: ac74_read_multilevel_param
  external :: soil_parm_ac74

  mtype = LIS_rc%lsm_index

  do n=1, LIS_rc%nnest
     ! allocate memory for place holder for #n nest
     allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))

     !------------------------------------!
     ! reading spatial spatial parameters !
     !------------------------------------!
     ! croptype
     ! Read crop type:

     call ac74_read_croptype(n)

     ! soiltype takes value from the LIS built-in parameter soilt
     if(LIS_rc%usetexturemap(n) .ne. 'none') then
        write(LIS_logunit,*) "AC74: retrieve parameter SOILTYPE from LIS"
        do t=1, LIS_rc%npatch(n, mtype)
           AC74_struc(n)%ac74(t)%soiltype= LIS_surface(n, mtype)%tile(t)%soilt
        enddo
     else
        ! read: soiltype
        write(LIS_logunit,*) "AC74: reading parameter SOILTYPE from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AC74_struc(n)%LDT_ncvar_soiltype), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           AC74_struc(n)%ac74(t)%soiltype = placeholder(col, row)
        enddo
     endif

     ! Height of the forcings (T and q, u and v), from the met forcing source.
     ! Only the height of the last forcing source is considered (overlay).
     if (trim(LIS_rc%metforc_blend_alg) .ne. "overlay") then
        write(LIS_logunit,*) "[ERR] AC74 only runs with overlay option for the met forcing blending."
        write(LIS_logunit,*) "[ERR] Program stopping ..."
        call LIS_endrun
     endif
     AC74_struc(n)%forchgt_tq = 2
     AC74_struc(n)%forchgt_uv = 2
     do m = 1, LIS_rc%nmetforc
        valid_forcing = .false.
#if (defined MF_MERRA2)
        if (LIS_rc%metforc(m) == "MERRA2") then
           ! 2 m fields only
           valid_forcing = (merra2_struc(n)%uselml == 0) .and. (merra2_struc(n)%use2mwind == 1)
           AC74_struc(n)%forchgt_tq = 2
           AC74_struc(n)%forchgt_uv = 2
        endif
#endif
#if (defined MF_ERA5)
        if (LIS_rc%metforc(m) == "ERA5") then
           valid_forcing = .true.
           AC74_struc(n)%forchgt_tq = 10
           AC74_struc(n)%forchgt_uv = 10
        endif
#endif
#if (defined MF_NLDAS2)
        if (LIS_rc%metforc(m) == "NLDAS2 grib") then
           ! 2 m (T and q) and 10 m (u and v) fields only
           valid_forcing = (nldas2_struc(n)%model_level_data == 0)
           AC74_struc(n)%forchgt_tq = 2
           AC74_struc(n)%forchgt_uv = 10
        endif
#endif
#if (defined MF_NLDAS20)
        if (LIS_rc%metforc(m) == "NLDAS2 netcdf") then
           ! 2 m (T and q) and 10 m (u and v) fields only
           valid_forcing = (nldas20_struc(n)%model_level_data == 0)
           AC74_struc(n)%forchgt_tq = 2
           AC74_struc(n)%forchgt_uv = 10
        endif
#endif
        if (.not. valid_forcing) then
           write(LIS_logunit,*) "[ERR] AC74 only runs with the following met forcings: "// &
                "ERA5, MERRA2 (2 m fields), NLDAS2 (2/10 m fields)"
           write(LIS_logunit,*) "[ERR] Program stopping ..."
           call LIS_endrun
        endif
        if (LIS_rc%met_ecor(m) /= "lapse-rate") then
           write(LIS_logunit,*) "[ERR] AC74 only runs with lapse-rate correction turned ON."
           write(LIS_logunit,*) "[ERR] Program stopping ..."
           call LIS_endrun
        endif
     enddo

     ! Tmin and Tmax climatology of each met forcing source (LDT), at the
     ! elevation of the forcing: lapse-rate corrected to 2 m above the tile,
     ! as the forcing itself. With overlay, the last valid source is kept.
     do m = 1, LIS_rc%nmetforc
        write(LIS_logunit,*) "[INFO] AC74: reading parameter AC_Tmin_clim_"//trim(LIS_rc%metforc(m))//" from ",&
             trim(LIS_rc%paramfile(n))
        do k = 1, 12
           call AC74_read_MULTILEVEL_param(n, trim(AC74_struc(n)%LDT_ncvar_tmincli_monthly)//"_"// &
                trim(LIS_rc%metforc(m)), k, placeholder)
           do t = 1, LIS_rc%npatch(n, mtype)
              col = LIS_surface(n, mtype)%tile(t)%col
              row = LIS_surface(n, mtype)%tile(t)%row
              if (placeholder(col, row) .ne. LIS_rc%udef) then
                 elevdiff = (LIS_domain(n)%tile(t)%elev + 2) &
                      - (LIS_forc(n,m)%modelelev(LIS_domain(n)%tile(t)%index) + AC74_struc(n)%forchgt_tq)
                 ! Climatology is rounded to 2 decimals in AquaCrop
                 AC74_struc(n)%ac74(t)%tmincli_monthly(k) = &
                      anint((placeholder(col, row) + lapse*elevdiff)*100)/100 - LIS_CONST_TKFRZ
              endif
           enddo
        enddo

        write(LIS_logunit,*) "[INFO] AC74: reading parameter AC_Tmax_clim_"//trim(LIS_rc%metforc(m))//" from ",&
             trim(LIS_rc%paramfile(n))
        do k = 1, 12
           call AC74_read_MULTILEVEL_param(n, trim(AC74_struc(n)%LDT_ncvar_tmaxcli_monthly)//"_"// &
                trim(LIS_rc%metforc(m)), k, placeholder)
           do t = 1, LIS_rc%npatch(n, mtype)
              col = LIS_surface(n, mtype)%tile(t)%col
              row = LIS_surface(n, mtype)%tile(t)%row
              if (placeholder(col, row) .ne. LIS_rc%udef) then
                 elevdiff = (LIS_domain(n)%tile(t)%elev + 2) &
                      - (LIS_forc(n,m)%modelelev(LIS_domain(n)%tile(t)%index) + AC74_struc(n)%forchgt_tq)
                 ! Climatology is rounded to 2 decimals in AquaCrop
                 AC74_struc(n)%ac74(t)%tmaxcli_monthly(k) = &
                      anint((placeholder(col, row) + lapse*elevdiff)*100)/100 - LIS_CONST_TKFRZ
              endif
           enddo
        enddo
     enddo

     ! Crop GDD parameters from the LIS param file, when all are there.
     ! Note: not produced by LDT, they are added externally to the file.
     ! They are only used where the crop file has 'External' for them.
     AC74_struc(n)%GDDfromLDT = &
          AC74_param_has_var(n, AC74_struc(n)%LDT_ncvar_gdd_cco) .and. &
          AC74_param_has_var(n, AC74_struc(n)%LDT_ncvar_gdd_maturity) .and. &
          AC74_param_has_var(n, AC74_struc(n)%LDT_ncvar_gdd_senescence) .and. &
          AC74_param_has_var(n, AC74_struc(n)%LDT_ncvar_gdd_maxr) .and. &
          AC74_param_has_var(n, AC74_struc(n)%LDT_ncvar_cgc) .and. &
          AC74_param_has_var(n, AC74_struc(n)%LDT_ncvar_cdc)
     if (AC74_struc(n)%GDDfromLDT) then
        write(LIS_logunit,*) "[INFO] AC74: reading crop GDD parameters from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AC74_struc(n)%LDT_ncvar_gdd_cco), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
           AC74_struc(n)%ac74(t)%GDDaysToGermination = nint(placeholder( &
                LIS_surface(n, mtype)%tile(t)%col, LIS_surface(n, mtype)%tile(t)%row))
        enddo
        call LIS_read_param(n, trim(AC74_struc(n)%LDT_ncvar_gdd_maturity), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
           AC74_struc(n)%ac74(t)%GDDaysToHarvest = nint(placeholder( &
                LIS_surface(n, mtype)%tile(t)%col, LIS_surface(n, mtype)%tile(t)%row))
        enddo
        call LIS_read_param(n, trim(AC74_struc(n)%LDT_ncvar_gdd_senescence), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
           AC74_struc(n)%ac74(t)%GDDaysToSenescence = nint(placeholder( &
                LIS_surface(n, mtype)%tile(t)%col, LIS_surface(n, mtype)%tile(t)%row))
        enddo
        call LIS_read_param(n, trim(AC74_struc(n)%LDT_ncvar_gdd_maxr), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
           AC74_struc(n)%ac74(t)%GDDaysToMaxRooting = nint(placeholder( &
                LIS_surface(n, mtype)%tile(t)%col, LIS_surface(n, mtype)%tile(t)%row))
        enddo
        call LIS_read_param(n, trim(AC74_struc(n)%LDT_ncvar_cgc), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
           AC74_struc(n)%ac74(t)%GDDCGC = placeholder( &
                LIS_surface(n, mtype)%tile(t)%col, LIS_surface(n, mtype)%tile(t)%row)
        enddo
        call LIS_read_param(n, trim(AC74_struc(n)%LDT_ncvar_cdc), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
           AC74_struc(n)%ac74(t)%GDDCDC = placeholder( &
                LIS_surface(n, mtype)%tile(t)%col, LIS_surface(n, mtype)%tile(t)%row)
        enddo
     endif

     ! Relative biomass (soil fertility, fraction) from the LIS param file.
     ! Note: not produced by LDT, added externally to the file. Without it,
     ! the soil fertility stress of the management file is used.
     AC74_struc(n)%BrelfromLDT = AC74_param_has_var(n, AC74_struc(n)%LDT_ncvar_brel)
     do t = 1, LIS_rc%npatch(n, mtype)
        AC74_struc(n)%ac74(t)%Brel = -1. ! not supplied
     enddo
     if (AC74_struc(n)%BrelfromLDT) then
        write(LIS_logunit,*) "[INFO] AC74: reading relative biomass from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AC74_struc(n)%LDT_ncvar_brel), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
           AC74_struc(n)%ac74(t)%Brel = placeholder( &
                LIS_surface(n, mtype)%tile(t)%col, LIS_surface(n, mtype)%tile(t)%row)
        enddo
     endif
     deallocate(placeholder)
     ! Read soil table
     call SOIL_PARM_AC74(AC74_struc(n)%soil_tbl_name)
     TheProjectType = typeproject_typeprm

     ! Create AquaCrop 'Run' Structure: a single run, filled per tile with
     ! its planned run (AC74_set_run_input) before each initialization
     call DetermineDayNr(LIS_rc%sda, LIS_rc%smo, LIS_rc%syr, AC74_struc(n)%DayNrLISStart)
     call DetermineDayNr(LIS_rc%eda, LIS_rc%emo, LIS_rc%eyr, AC74_struc(n)%DayNrLISEnd)

     ! A coldstart with a planting criterion cannot start within the window
     if ((trim(LIS_rc%startcode) .eq. "coldstart") .and. &
          (AC74_struc(n)%Temp_crit .or. AC74_struc(n)%Rainfall_crit)) then
        call DetermineDayNr(AC74_struc(n)%Crop_AnnualStartDay, &
             AC74_struc(n)%Crop_AnnualStartMonth, LIS_rc%syr, WindowStart)
        if ((AC74_struc(n)%DayNrLISStart > WindowStart) .and. &
             (AC74_struc(n)%DayNrLISStart < WindowStart + AC74_struc(n)%crit_window)) then
           write(LIS_logunit, *) '[ERR] AC74 coldstart: LIS starts within the planting window'
           write(LIS_logunit, *) '[ERR] start LIS on or before the starting day of crop period'
           call LIS_endrun()
        endif
     endif

     call allocate_project_input(1)
     l = 1
     call set_project_input(l, 'Simulation_YearSeason', 1_int8)
     ! Note: '(External)' input sources can vary spatially (e.g. soil, crop, meteo)
     call set_project_input(l, 'Description', ' LIS ')
     call set_project_input(l, 'Climate_Info', '(External)')
     call set_project_input(l, 'Climate_Filename', '(External)')
     call set_project_input(l, 'Climate_Directory', '(None)')
     call set_project_input(l, 'VersionNr', 7.4)
     call set_project_input(l, 'Temperature_Info', '(None)')
     call set_project_input(l, 'Temperature_Filename', '(External)')
     call set_project_input(l, 'Temperature_Directory', '(None)')
     call set_project_input(l, 'ETo_Info', '(None)')
     call set_project_input(l, 'ETo_Filename', '(External)')
     call set_project_input(l, 'ETo_Directory', '(None)')
     call set_project_input(l, 'Rain_Info', ' LIS ')
     call set_project_input(l, 'Rain_Filename', '(External)')
     call set_project_input(l, 'Rain_Directory', '(None)')
     call set_project_input(l, 'CO2_Info', ' LIS ')
     call set_project_input(l, 'CO2_Filename', trim(AC74_struc(n)%CO2_Filename))
     call set_project_input(l, 'CO2_Directory', trim(AC74_struc(n)%PathNameSimul))
     call set_project_input(l, 'Calendar_Info', '(None)')
     call set_project_input(l, 'Calendar_Filename', '(None)')
     call set_project_input(l, 'Calendar_Directory', '(None)')
     call set_project_input(l, 'Crop_Info', '(None)')
     call set_project_input(l, 'Crop_Filename', '(External)')
     call set_project_input(l, 'Crop_Directory', trim(AC74_struc(n)%PathCropFiles))
     call set_project_input(l, 'Irrigation_Info', ' LIS ')
     call set_project_input(l, 'Irrigation_Filename', trim(AC74_struc(n)%Irrigation_Filename))
     call set_project_input(l, 'Irrigation_Directory', trim(AC74_struc(n)%PathNameSimul))
     call set_project_input(l, 'Management_Info', ' LIS ')
     call set_project_input(l, 'Management_Filename',  trim(AC74_struc(n)%Management_Filename))
     call set_project_input(l, 'Management_Directory', trim(AC74_struc(n)%PathNameSimul))
     call set_project_input(l, 'GroundWater_Info', '(None)')
     call set_project_input(l, 'GroundWater_Filename', '(None)')
     call set_project_input(l, 'GroundWater_Directory', '(None)')
     call set_project_input(l, 'Soil_Info', '(External)')
     call set_project_input(l, 'Soil_Filename', '(External)')
     call set_project_input(l, 'Soil_Directory', '(External)')
     call set_project_input(l, 'SWCIni_Info', '(None)')
     call set_project_input(l, 'SWCIni_Directory', '(None)')
     call set_project_input(l, 'OffSeason_Info', '(None)')
     call set_project_input(l, 'OffSeason_Filename', '(None)')
     call set_project_input(l, 'OffSeason_Directory', '(None)')
     call set_project_input(l, 'Observations_Info', '(None)')
     call set_project_input(l, 'Observations_Filename', '(None)')
     call set_project_input(l, 'Observations_Directory', '(None)')

     do t = 1, LIS_rc%npatch(n, mtype)

        col = LIS_surface(n, mtype)%tile(t)%col
        row = LIS_surface(n, mtype)%tile(t)%row

        ! Initialize global strings
        call InitializeGlobalStrings()

        call SetPathNameSimul(trim(AC74_struc(n)%PathNameSimul)) ! needed for InitializeSettings
        ! Initializes all settings to default, later overwritten in main
        call InitializeSettings(use_default_soil_file=.false.,use_default_crop_file=.false.)
        ! It uses the defaults anyway
        ! First run of the tile, from the start of LIS
        ! (crop parameters are read when calling InitializeRunPart1)
        AC74_struc(n)%ac74(t)%irun = 1
        call AC74_plan_run(n, t, AC74_struc(n)%DayNrLISStart)
        call AC74_set_run_input(n, t)
        call CheckForKeepSWC(MultipleRunWithKeepSWC_temp, &
             MultipleRunConstZrx_temp)
        call SetSimulation_MultipleRunWithKeepSWC(MultipleRunWithKeepSWC_temp)
        call SetSimulation_MultipleRunConstZrx(MultipleRunConstZrx_temp)
!!! Start the Program AC74
        call SetPathNameOutp('')
        call SetPathNameList('')
        call SetPathNameParam('')
        call SetPathNameProg('')


        call GetTimeAggregationResults()
        call GetRequestDailyResults()
        call GetRequestParticularResults()

        call SetClimFile('(External)')
        call SetCO2Description('')
        call SetProfFilefull('External')

        call SetClimRecord_DataType(0_int8)
        call SetClimRecord_fromd(LIS_rc%sda)
        call SetClimRecord_fromdaynr(AC74_struc(n)%DayNrLISStart)
        call SetClimRecord_fromm(LIS_rc%smo)
        call SetClimRecord_fromstring('any date')
        call SetClimRecord_fromy(LIS_rc%syr)
        call SetClimRecord_NrObs(999)
        call SetClimRecord_tod(LIS_rc%eda)
        call SetClimRecord_todaynr(AC74_struc(n)%DayNrLISEnd)
        call SetClimRecord_tom(LIS_rc%emo)
        call SetClimRecord_tostring('any date')
        call SetClimRecord_toy(LIS_rc%eyr)

        ! SetTemperatureRecord
        call SetTemperatureRecord_DataType(0_int8)
        call SetTemperatureRecord_fromd(LIS_rc%sda)
        call SetTemperatureRecord_fromdaynr(AC74_struc(n)%DayNrLISStart)
        call SetTemperatureRecord_fromm(LIS_rc%smo)
        call SetTemperatureRecord_fromstring('any date')
        call SetTemperatureRecord_fromy(LIS_rc%syr)
        call SetTemperatureRecord_NrObs(999)
        call SetTemperatureRecord_tod(LIS_rc%eda)
        call SetTemperatureRecord_todaynr(AC74_struc(n)%DayNrLISEnd)
        call SetTemperatureRecord_tom(LIS_rc%emo)
        call SetTemperatureRecord_tostring('any date')
        call SetTemperatureRecord_toy(LIS_rc%eyr)

        ! SetEToRecord
        call SetEToRecord_DataType(0_int8)
        call SetEToRecord_fromd(LIS_rc%sda)
        call SetEToRecord_fromdaynr(AC74_struc(n)%DayNrLISStart)
        call SetEToRecord_fromm(LIS_rc%smo)
        call SetEToRecord_fromstring('any date')
        call SetEToRecord_fromy(LIS_rc%syr)
        call SetEToRecord_NrObs(999)
        call SetEToRecord_tod(LIS_rc%eda)
        call SetEToRecord_todaynr(AC74_struc(n)%DayNrLISEnd)
        call SetEToRecord_tom(LIS_rc%emo)
        call SetEToRecord_tostring('any date')
        call SetEToRecord_toy(LIS_rc%eyr)

        ! SetRainRecord
        call SetRainRecord_DataType(0_int8)
        call SetRainRecord_fromd(LIS_rc%sda)
        call SetRainRecord_fromdaynr(AC74_struc(n)%DayNrLISStart)
        call SetRainRecord_fromm(LIS_rc%smo)
        call SetRainRecord_fromstring('any date')
        call SetRainRecord_fromy(LIS_rc%syr)
        call SetRainRecord_NrObs(999)
        call SetRainRecord_tod(LIS_rc%eda)
        call SetRainRecord_todaynr(AC74_struc(n)%DayNrLISEnd)
        call SetRainRecord_tom(LIS_rc%emo)
        call SetRainRecord_tostring('any date')
        call SetRainRecord_toy(LIS_rc%eyr)

        call SetSimulation_Storage_CropString('None')

        call SetMultipleProjectDescription('undefined')
        ProgramParametersAvailable = .false.
        call SetFullFileNameProgramParameters("(None)")
        call LoadProgramParametersProjectPlugIn(&
             GetFullFileNameProgramParameters(), &
             ProgramParametersAvailable)
        call SetSimulation_MultipleRun(.true.)
        call SetSimulation_NrRuns(1)
        MultipleRunWithKeepSWC_temp = GetSimulation_MultipleRunWithKeepSWC()
        MultipleRunConstZrx_temp = GetSimulation_MultipleRunConstZrx()

        ! Set AC74 soil parameters based on soiltype
        ! Note: up to 5 soil layers with the same texture
        AC74_struc(n)%ac74(t)%SoilLayer(1)%wp = WP(AC74_struc(n)%ac74(t)%soiltype) * 100
        AC74_struc(n)%ac74(t)%SoilLayer(1)%sat = SAT(AC74_struc(n)%ac74(t)%soiltype) * 100
        AC74_struc(n)%ac74(t)%SoilLayer(1)%fc = FC(AC74_struc(n)%ac74(t)%soiltype) * 100
        AC74_struc(n)%ac74(t)%SoilLayer(1)%InfRate = INFRATE(AC74_struc(n)%ac74(t)%soiltype)
        sd_tmp = sd(AC74_struc(n)%ac74(t)%soiltype)*100
        cl_tmp = cl(AC74_struc(n)%ac74(t)%soiltype)*100
        si_tmp = si(AC74_struc(n)%ac74(t)%soiltype)*100

        ! define default CN
        if (AC74_struc(n)%ac74(t)%SoilLayer(1)%InfRate>864) then
           AC74_struc(n)%ac74(t)%Soil%CNvalue = 46
        elseif (AC74_struc(n)%ac74(t)%SoilLayer(1)%InfRate>=347) then
           AC74_struc(n)%ac74(t)%Soil%CNvalue = 61
        elseif (AC74_struc(n)%ac74(t)%SoilLayer(1)%InfRate>=36) then
           AC74_struc(n)%ac74(t)%Soil%CNvalue = 72
        else
           AC74_struc(n)%ac74(t)%Soil%CNvalue = 77
        endif

        ! define default REW (only top layer will be used)
        Z_surf = 0.04
        REW=nint(10.0*(AC74_struc(n)%ac74(t)%SoilLayer(1)%fc-(AC74_struc(n)%ac74(t)%SoilLayer(1)%wp/2.0))*Z_surf)
        if (REW < 0) REW = 0
        if (REW > 15) REW = 15
        AC74_struc(n)%ac74(t)%Soil%REW = REW

        !  associate soil class with USDA soil type for soil description
        if ((1.5 * cl_tmp + si_tmp) < 15.0) then
           descr = 0
        elseif (((1.5 * cl_tmp + si_tmp) >= 15) .and. ((2 * cl_tmp + si_tmp) <= 30)) then
           descr = 1
        elseif ((cl_tmp >= 7 .and. cl_tmp < 20 .and. sd_tmp >= 52) .and. (2 * cl_tmp + si_tmp >= 30)) then
           descr = 2
        elseif ((cl_tmp < 7) .and. (si_tmp < 50) .and. (2 * cl_tmp + si_tmp >= 30)) then
           descr = 2
        elseif ((cl_tmp >= 7) .and. (cl_tmp < 27) .and. (si_tmp >= 28) .and. (si_tmp < 50) .and. (sd_tmp < 52)) then
           descr = 3
        elseif ((si_tmp >= 50) .and. (cl_tmp >= 12) .and. (cl_tmp < 27)) then
           descr = 4
        elseif ((si_tmp >= 50) .and. (si_tmp < 80) .and. (cl_tmp < 12)) then
           descr = 4
        elseif ((si_tmp >= 80) .and. (cl_tmp < 12)) then
           descr = 5
        elseif ((cl_tmp >= 20) .and. (cl_tmp < 35) .and. (si_tmp < 28) .and. (sd_tmp > 45)) then
           descr = 6
        elseif ((cl_tmp >= 27) .and. (cl_tmp < 40) .and. (sd_tmp >= 20) .and. (sd_tmp < 45)) then
           descr = 7
        elseif ((cl_tmp >= 27) .and. (cl_tmp < 40) .and. (sd_tmp < 20)) then
           descr = 8
        elseif ((cl_tmp >= 35) .and. (cl_tmp < 55) .and. (sd_tmp >= 45) .and. (sd_tmp < 65)) then
           descr = 9
        elseif ((cl_tmp >= 40) .and. (si_tmp >= 40)) then
           descr = 10
        elseif ((cl_tmp >= 40) .and. (sd_tmp < 45) .and. (si_tmp < 40)) then
           descr = 11
        else
           write(LIS_logunit, *) 'no soil texture found'
           write(LIS_logunit, *) 'program stopping ...'
           call LIS_endrun
        end if
        InfRate_tmp = AC74_struc(n)%ac74(t)%SoilLayer(1)%InfRate
        if ((descr == 0) .or. (descr == 1) .or. (descr == 2)) then
           AC74_struc(n)%ac74(t)%SoilLayer(1)%CRa = -0.3112-10**(-5)*InfRate_tmp
           AC74_struc(n)%ac74(t)%SoilLayer(1)%CRb = -1.4936+0.2416*log(InfRate_tmp)
        elseif ((descr == 3) .or. (descr == 4) .or. (descr == 5)) then
           AC74_struc(n)%ac74(t)%SoilLayer(1)%CRa = -0.4986-9*10**(-5)*InfRate_tmp
           AC74_struc(n)%ac74(t)%SoilLayer(1)%CRb = -2.1320+0.4778*log(InfRate_tmp)
        elseif ((descr == 6) .or. (descr == 7) .or. (descr == 9)) then
           AC74_struc(n)%ac74(t)%SoilLayer(1)%CRa = -0.5677-4*10**(-5)*InfRate_tmp
           AC74_struc(n)%ac74(t)%SoilLayer(1)%CRb = -3.7189+0.5922*log(InfRate_tmp)
        elseif ((descr == 8) .or. (descr == 10) .or. (descr == 11)) then
           AC74_struc(n)%ac74(t)%SoilLayer(1)%CRa = -0.6366+8*10**(-4)*InfRate_tmp
           AC74_struc(n)%ac74(t)%SoilLayer(1)%CRb = -1.9165+0.7163*log(InfRate_tmp)
        endif
        ! Set GravelMass and Penetrability
        AC74_struc(n)%ac74(t)%SoilLayer(1)%Penetrability = 100.0
        AC74_struc(n)%ac74(t)%SoilLayer(1)%GravelMass = 10.0
        ! Set GravelVol
        AC74_struc(n)%ac74(t)%SoilLayer(1)%GravelVol = &
             FromGravelMassToGravelVolume(AC74_struc(n)%ac74(t)%SoilLayer(1)%sat,&
             AC74_struc(n)%ac74(t)%SoilLayer(1)%GravelMass)

        AC74_struc(n)%ac74(t)%SoilLayer(1)%Description = 'soil type from LIS'


        AC74_struc(n)%ac74(t)%SoilLayer(1)%Thickness = AC74_struc(n)%Thickness(1)

        ! loop if more than 2 layers
        if (AC74_struc(n)%NrSoilLayers.gt.1) then
           do l = 2, AC74_struc(n)%NrSoilLayers
              AC74_struc(n)%ac74(t)%SoilLayer(l) = AC74_struc(n)%ac74(t)%SoilLayer(1)
              AC74_struc(n)%ac74(t)%SoilLayer(l)%Thickness = AC74_struc(n)%Thickness(l)
              AC74_struc(n)%ac74(t)%SoilLayer(l)%GravelMass = 0.0
              AC74_struc(n)%ac74(t)%SoilLayer(l)%GravelVol = 0.0
           enddo
        endif

        if (AC74_struc(n)%NrSoilLayers.lt.5) then !repeat the last layer
           do l = AC74_struc(n)%NrSoilLayers+1, 5
              AC74_struc(n)%ac74(t)%SoilLayer(l) = AC74_struc(n)%ac74(t)%SoilLayer(AC74_struc(n)%NrSoilLayers)
           enddo
        endif

        AC74_struc(n)%ac74(t)%Soil%NrSoilLayers = AC74_struc(n)%NrSoilLayers
        AC74_struc(n)%ac74(t)%NrCompartments = AC74_struc(n)%max_No_compartments

        ! Set soil global variables
        call SetSoilLayer(AC74_struc(n)%ac74(t)%SoilLayer)
        call SetSoil(AC74_struc(n)%ac74(t)%Soil)
        call SetNrCompartments(AC74_struc(n)%ac74(t)%NrCompartments)

        call SetPathNameProg('')
        !
        call SetMultipleProjectDescription('undefined')
        ProgramParametersAvailable = .false.
        ! Running with default program parameters
        call SetFullFileNameProgramParameters("(None)")
        call LoadProgramParametersProjectPlugIn(&
             GetFullFileNameProgramParameters(), &
             ProgramParametersAvailable)
        call SetSimulation_MultipleRun(.true.)
        call SetSimulation_NrRuns(1)

        call SetClimRecord_DataType(0_int8)
        call SetClimRecord_fromd(0)
        call SetClimRecord_fromdaynr(AC74_struc(n)%DayNrLISStart)
        call SetClimRecord_fromm(0)
        call SetClimRecord_fromstring('any date')
        call SetClimRecord_fromy(LIS_rc%syr)
        call SetClimRecord_NrObs(999)
        call SetClimRecord_tod(0)
        call SetClimRecord_todaynr(AC74_struc(n)%DayNrLISEnd)
        call SetClimRecord_tom(0)
        call SetClimRecord_tostring('any date')
        call SetClimRecord_toy(0)
        call SetClimFile('(External)')

        AC74_struc(n)%ac74(t)%WPi = 0.

        ! No forcing yet: climatology of the start month for day 1
        ! (a coldstart initializes run 1 again in main with the forcing)
        call SetTmin(AC74_struc(n)%ac74(t)%tmincli_monthly(LIS_rc%smo))
        call SetTmax(AC74_struc(n)%ac74(t)%tmaxcli_monthly(LIS_rc%smo))

        ! Set Tmin and Tmax reference to compute the stress realtions
        call SetTminTnxReference12MonthsRun(AC74_struc(n)%ac74(t)%tmincli_monthly(:))
        call SetTmaxTnxReference12MonthsRun(AC74_struc(n)%ac74(t)%tmaxcli_monthly(:))

        ! Set reference year for CO2 for stress functions
        call SetTnxReferenceYear(AC74_struc(n)%tempcli_refyr)
        call SetTnxReferenceFile('(External)')

        ! InitializeRunPart1 (the single project run)
        call InitializeRunPart1(1_int8, TheProjectType)
        call InitializeSimulationRunPart2()
        AC74_struc(n)%ac74(t)%InitializeRun = 0
        AC74_struc(n)%AC74(t)%crop = GetCrop()
        AC74_struc(n)%ac74(t)%SeasonNr = 0
        AC74_struc(n)%ac74(t)%SeasonYield = 0.
        AC74_struc(n)%ac74(t)%CropActive = 0.
        ! Close irrigation file after Run Initialization
        ! Note: only 2 irrigation records can be passed
        if((GetIrriMode().eq.IrriMode_Generate)&
             .or.(GetIrriMode().eq.IrriMode_Manual)) then
           call fIrri_close()
        endif

        ! Set AC74_struc after Initialization
        AC74_struc(n)%AC74(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
        AC74_struc(n)%AC74(t)%RootZoneWC_FC = GetRootZoneWC_FC()
        AC74_struc(n)%AC74(t)%RootZoneWC_WP = GetRootZoneWC_WP()
        AC74_struc(n)%AC74(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
        AC74_struc(n)%AC74(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
        AC74_struc(n)%AC74(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
        AC74_struc(n)%AC74(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
        AC74_struc(n)%AC74(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
        AC74_struc(n)%AC74(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
        AC74_struc(n)%AC74(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
        AC74_struc(n)%AC74(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
        AC74_struc(n)%AC74(t)%Compartment = GetCompartment()
        AC74_struc(n)%AC74(t)%TotalSaltContent = GetTotalSaltContent()
        AC74_struc(n)%AC74(t)%TotalWaterContent = GetTotalWaterContent()
        AC74_struc(n)%AC74(t)%effectiverain = Geteffectiverain()
        AC74_struc(n)%AC74(t)%SumWaBal = GetSumWaBal()
        AC74_struc(n)%AC74(t)%RootZoneSalt = GetRootZoneSalt()
        AC74_struc(n)%AC74(t)%Simulation = GetSimulation()
        AC74_struc(n)%AC74(t)%IrriInterval = GetIrriInterval()
        AC74_struc(n)%AC74(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
        AC74_struc(n)%AC74(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
        AC74_struc(n)%AC74(t)%Irrigation = GetIrrigation()
        AC74_struc(n)%AC74(t)%IrriBeforeSeason = GetIrriBeforeSeason()
        AC74_struc(n)%AC74(t)%IrriAfterSeason = GetIrriAfterSeason()
        AC74_struc(n)%AC74(t)%SoilLayer = GetSoilLayer()
        AC74_struc(n)%AC74(t)%daynri = GetDayNri()
        AC74_struc(n)%AC74(t)%NrCompartments = GetNrCompartments()
        do l=1, AC74_struc(n)%AC74(t)%NrCompartments
           AC74_struc(n)%AC74(t)%smc(l) = GetCompartment_theta(l)
        enddo
        AC74_struc(n)%AC74(t)%IrriECw = GetIrriECw()
        AC74_struc(n)%AC74(t)%Management = GetManagement()
        AC74_struc(n)%AC74(t)%PerennialPeriod = GetPerennialPeriod()
        AC74_struc(n)%AC74(t)%simulparam = GetSimulParam()
        AC74_struc(n)%AC74(t)%Cuttings = GetManagement_Cuttings()
        AC74_struc(n)%AC74(t)%onset = GetOnset()
        AC74_struc(n)%AC74(t)%endseason = GetEndSeason()
        AC74_struc(n)%AC74(t)%Soil = GetSoil()
        AC74_struc(n)%AC74(t)%TemperatureRecord = GetTemperatureRecord()
        AC74_struc(n)%AC74(t)%ClimRecord = GetClimRecord()
        AC74_struc(n)%AC74(t)%RainRecord = GetRainRecord()
        AC74_struc(n)%AC74(t)%EToRecord = GetEToRecord()

        AC74_struc(n)%AC74(t)%GenerateTimeMode = GetGenerateTimeMode()
        AC74_struc(n)%AC74(t)%GenerateDepthMode = GetGenerateDepthMode()
        AC74_struc(n)%AC74(t)%IrriMode = GetIrriMode()
        AC74_struc(n)%AC74(t)%IrriMethod = GetIrriMethod()
        AC74_struc(n)%AC74(t)%DaySubmerged = GetDaySubmerged()
        AC74_struc(n)%AC74(t)%MaxPlotNew = GetMaxPlotNew()
        AC74_struc(n)%AC74(t)%IrriFirstDayNr = GetIrriFirstDayNr()
        AC74_struc(n)%AC74(t)%ZiAqua = GetZiAqua()
        AC74_struc(n)%AC74(t)%IniPercTAW = GetIniPercTAW()
        AC74_struc(n)%AC74(t)%MaxPlotTr = GetMaxPlotTr()
        AC74_struc(n)%AC74(t)%OutputAggregate = GetOutputAggregate()

        AC74_struc(n)%AC74(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
        AC74_struc(n)%AC74(t)%PreDay = GetPreDay()
        AC74_struc(n)%AC74(t)%OutDaily = GetOutDaily()
        AC74_struc(n)%AC74(t)%Out1Wabal = GetOut1Wabal()
        AC74_struc(n)%AC74(t)%Out2Crop = GetOut2Crop()
        AC74_struc(n)%AC74(t)%Out3Prof = GetOut3Prof()
        AC74_struc(n)%AC74(t)%Out4Salt = GetOut4Salt()
        AC74_struc(n)%AC74(t)%Out5CompWC = GetOut5CompWC()
        AC74_struc(n)%AC74(t)%Out6CompEC = GetOut6CompEC()
        AC74_struc(n)%AC74(t)%Out7Clim = GetOut7Clim()
        AC74_struc(n)%AC74(t)%Out8Irri = GetOut8Irri()
        AC74_struc(n)%AC74(t)%Part1Mult = GetPart1Mult()
        AC74_struc(n)%AC74(t)%Part2Eval = GetPart2Eval()

        !
        AC74_struc(n)%AC74(t)%CCiActual = GetCCiActual()
        AC74_struc(n)%AC74(t)%CCiprev = GetCCiprev()
        AC74_struc(n)%AC74(t)%CCiTopEarlySen = GetCCiTopEarlySen()
        AC74_struc(n)%AC74(t)%CRsalt = GetCRsalt ()
        AC74_struc(n)%AC74(t)%CRwater = GetCRwater()
        AC74_struc(n)%AC74(t)%ECdrain = GetECdrain()
        AC74_struc(n)%AC74(t)%ECiAqua = GetECiAqua()
        AC74_struc(n)%AC74(t)%ECstorage = GetECstorage()
        AC74_struc(n)%AC74(t)%Eact = GetEact()
        AC74_struc(n)%AC74(t)%Epot = GetEpot()
        AC74_struc(n)%AC74(t)%ETo = GetETo()
        AC74_struc(n)%AC74(t)%Drain = GetDrain()
        AC74_struc(n)%AC74(t)%Infiltrated = GetInfiltrated()
        AC74_struc(n)%AC74(t)%prcp = GetRain()
        AC74_struc(n)%AC74(t)%RootingDepth = GetRootingDepth()
        AC74_struc(n)%AC74(t)%Runoff = GetRunoff()
        AC74_struc(n)%AC74(t)%SaltInfiltr = GetSaltInfiltr()
        AC74_struc(n)%AC74(t)%Surf0 = GetSurf0()
        AC74_struc(n)%AC74(t)%SurfaceStorage = GetSurfaceStorage()
        AC74_struc(n)%AC74(t)%Tact = GetTact()
        AC74_struc(n)%AC74(t)%Tpot = GetTpot()
        AC74_struc(n)%AC74(t)%TactWeedInfested = 0. !not ini in AC GetTactWeedInfested()
        AC74_struc(n)%AC74(t)%Tmax = GetTmax()
        AC74_struc(n)%AC74(t)%Tmin = GetTmin()


        AC74_struc(n)%AC74(t)%GwTable = GetGwTable()
        AC74_struc(n)%AC74(t)%PlotVarCrop = GetPlotVarCrop()
        AC74_struc(n)%AC74(t)%StressTot = GetStressTot()
        AC74_struc(n)%AC74(t)%CutInfoRecord1 = GetCutInfoRecord1()
        AC74_struc(n)%AC74(t)%CutInfoRecord2 = GetCutInfoRecord2()
        AC74_struc(n)%AC74(t)%Transfer = GetTransfer()
        AC74_struc(n)%AC74(t)%PreviousSum = GetPreviousSum()
        AC74_struc(n)%AC74(t)%Tadj = GetTadj()
        AC74_struc(n)%AC74(t)%GDDTadj = GetGDDTadj()
        AC74_struc(n)%AC74(t)%DayLastCut = GetDayLastCut()
        AC74_struc(n)%AC74(t)%NrCut = GetNrCut()
        AC74_struc(n)%AC74(t)%SumInterval = GetSumInterval()
        AC74_struc(n)%AC74(t)%PreviousStressLevel = GetPreviousStressLevel()
        AC74_struc(n)%AC74(t)%StressSFadjNEW = GetStressSFadjNEW()
        AC74_struc(n)%AC74(t)%Bin = GetBin()
        AC74_struc(n)%AC74(t)%Bout = GetBout()
        AC74_struc(n)%AC74(t)%GDDayi = GetGDDayi()
        AC74_struc(n)%AC74(t)%CO2i = GetCO2i()
        AC74_struc(n)%AC74(t)%FracBiomassPotSF = GetFracBiomassPotSF()
        AC74_struc(n)%AC74(t)%SumETo = GetSumETo()
        AC74_struc(n)%AC74(t)%SumGDD = GetSumGDD()
        AC74_struc(n)%AC74(t)%Ziprev = GetZiprev()
        AC74_struc(n)%AC74(t)%SumGDDPrev = GetSumGDDPrev()
        AC74_struc(n)%AC74(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
        AC74_struc(n)%AC74(t)%Coeffb0 = GetCoeffb0()
        AC74_struc(n)%AC74(t)%Coeffb1 = GetCoeffb1()
        AC74_struc(n)%AC74(t)%Coeffb2 = GetCoeffb2()
        AC74_struc(n)%AC74(t)%Coeffb0Salt = GetCoeffb0Salt()
        AC74_struc(n)%AC74(t)%Coeffb1Salt = GetCoeffb1Salt()
        AC74_struc(n)%AC74(t)%Coeffb2Salt = GetCoeffb2Salt()
        AC74_struc(n)%AC74(t)%StressLeaf = GetStressLeaf()
        AC74_struc(n)%AC74(t)%StressSenescence = GetStressSenescence()
        AC74_struc(n)%AC74(t)%DayFraction = GetDayFraction()
        AC74_struc(n)%AC74(t)%GDDayFraction = GetGDDayFraction()
        AC74_struc(n)%AC74(t)%CGCref = GetCGCref()
        AC74_struc(n)%AC74(t)%GDDCGCref = GetGDDCGCref()
        AC74_struc(n)%AC74(t)%TimeSenescence = GetTimeSenescence()
        AC74_struc(n)%AC74(t)%SumKcTop = GetSumKcTop()
        AC74_struc(n)%AC74(t)%SumKcTopStress = GetSumKcTopStress()
        AC74_struc(n)%AC74(t)%SumKci = GetSumKci()
        AC74_struc(n)%AC74(t)%CCoTotal = GetCCoTotal()
        AC74_struc(n)%AC74(t)%CCxTotal = GetCCxTotal()
        AC74_struc(n)%AC74(t)%CDCTotal = GetCDCTotal()
        AC74_struc(n)%AC74(t)%GDDCDCTotal = GetGDDCDCTotal()
        AC74_struc(n)%AC74(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
        AC74_struc(n)%AC74(t)%WeedRCi = GetWeedRCi()
        AC74_struc(n)%AC74(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
        AC74_struc(n)%AC74(t)%fWeedNoS = GetfWeedNoS()
        AC74_struc(n)%AC74(t)%Zeval = GetZeval()
        AC74_struc(n)%AC74(t)%BprevSum = GetBprevSum()
        AC74_struc(n)%AC74(t)%YprevSum = GetYprevSum()
        AC74_struc(n)%AC74(t)%SumGDDcuts = GetSumGDDcuts()
        AC74_struc(n)%AC74(t)%HItimesBEF = GetHItimesBEF()
        AC74_struc(n)%AC74(t)%ScorAT1 = GetScorAT1()
        AC74_struc(n)%AC74(t)%ScorAT2 = GetScorAT2()
        AC74_struc(n)%AC74(t)%HItimesAT1 = GetHItimesAT1()
        AC74_struc(n)%AC74(t)%HItimesAT2 = GetHItimesAT2()
        AC74_struc(n)%AC74(t)%HItimesAT = GetHItimesAT()
        AC74_struc(n)%AC74(t)%alfaHI = GetalfaHI()
        AC74_struc(n)%AC74(t)%alfaHIAdj = GetalfaHIAdj()
        AC74_struc(n)%AC74(t)%NextSimFromDayNr = GetNextSimFromDayNr()
        AC74_struc(n)%AC74(t)%DayNr1Eval = GetDayNr1Eval()
        AC74_struc(n)%AC74(t)%DayNrEval = GetDayNrEval()
        AC74_struc(n)%AC74(t)%LineNrEval = GetLineNrEval()
        AC74_struc(n)%AC74(t)%PreviousSumETo = GetPreviousSumETo()
        AC74_struc(n)%AC74(t)%PreviousSumGDD = GetPreviousSumGDD()
        AC74_struc(n)%AC74(t)%PreviousBmob = GetPreviousBmob()
        AC74_struc(n)%AC74(t)%PreviousBsto = GetPreviousBsto()
        AC74_struc(n)%AC74(t)%StageCode = GetStageCode()
        AC74_struc(n)%AC74(t)%PreviousDayNr = GetPreviousDayNr()
        AC74_struc(n)%AC74(t)%NoYear = GetNoYear()
        AC74_struc(n)%AC74(t)%WaterTableInProfile = GetWaterTableInProfile()
        AC74_struc(n)%AC74(t)%StartMode = GetStartMode()
        AC74_struc(n)%AC74(t)%NrRuns = GetSimulation_NrRuns()
        AC74_struc(n)%AC74(t)%TheProjectType = TheProjectType
        AC74_struc(n)%AC74(t)%DayNrCropEnded = GetDayNrCropEnded()
        AC74_struc(n)%AC74(t)%IrriInfoLastDay = GetIrriInfoLastDay()
        AC74_struc(n)%AC74(t)%tDaysZmin = tDaysZmin
        AC74_struc(n)%AC74(t)%tGDDZmin = tGDDZmin

        ! Check for irrigation (irrigation file management)
        if(AC74_struc(n)%ac74(t)%IrriMode.eq.IrriMode_Manual)then
           if(AC74_struc(n)%ac74(t)%IrriInfoRecord1%NoMoreInfo)then
              AC74_struc(n)%ac74(t)%irri_lnr = 9
           else
              AC74_struc(n)%ac74(t)%irri_lnr = 10
           endif
        elseif(AC74_struc(n)%ac74(t)%IrriMode.eq.IrriMode_Generate)then
           if(AC74_struc(n)%ac74(t)%IrriInfoRecord1%NoMoreInfo)then
              AC74_struc(n)%ac74(t)%irri_lnr = 11 + AC74_struc(n)%IrriLastDayLine
           else
              AC74_struc(n)%ac74(t)%irri_lnr = 12 + AC74_struc(n)%IrriLastDayLine
           endif
        else ! no irrigation, set to 0
           AC74_struc(n)%ac74(t)%irri_lnr = 0
        endif

        ! On a coldstart, initialize run 1 again in main with the forcing of
        ! day 1 (a restart reads the run of the tile from the restart file)
        if (trim(LIS_rc%startcode) .eq. "coldstart") then
           AC74_struc(n)%ac74(t)%InitializeRun = 1
        endif

     enddo ! do t = 1, LIS_rc%npatch(n, mtype)
  enddo
end subroutine AC74_setup

!BOP
!
! !ROUTINE: AC74_read_MULTILEVEL_param
!  \label{AC74_read_MULTILEVEL_param}
!
! !REVISION HISTORY:
!  04 NOV Louise Busschaert; Initial implementation
!
! !INTERFACE:
subroutine AC74_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)

  ! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_localPet,   &
       LIS_ews_halo_ind, LIS_ewe_halo_ind, &
       LIS_nss_halo_ind, LIS_nse_halo_ind
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  ! !ARGUMENTS:
  integer, intent(in)          :: n
  character(len=*), intent(in) :: ncvar_name
  integer, intent(in)          :: level
  real, intent(out)            :: placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n))
  ! !DESCRIPTION:
  !  This subroutine reads MULTILEVEL parameters from the LIS
  !  NetCDF parameter data file
  !
  !  The arguments are:
  !  \begin{description}
  !   \item[n]
  !    index of n
  !   \item[level]
  !    level index (month, quarter, soil layer, snow layer) of the data to be read
  !   \item[array]
  !    array containing returned values
  !   \end{description}
  !
  !EOP

  integer       :: ios, nid, param_ID, nc_ID, nr_ID, dimids(3)
  integer       :: nc, nr, nlevel
  real, pointer :: level_data(:, :, :)
  logical       :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then
     write(LIS_logunit, *) 'Reading '//trim(ncvar_name)//' map for level ', level

     ! open NetCDF parameter file
     ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
     call LIS_verify(ios, 'Error in nf90_open in AC74_read_MULTILEVEL_param')

     ! inquire the ID of east-west dimension
     ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
     call LIS_verify(ios, 'Error in nf90_inq_dimid in AC74_read_MULTILEVEL_param')

     ! inquire the ID of north-south dimension
     ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
     call LIS_verify(ios, 'Error in nf90_inq_dimid in AC74_read_MULTILEVEL_param')

     ! inquire the length of east-west dimension
     ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
     call LIS_verify(ios, 'Error in nf90_inquire_dimension in AC74_read_MULTILEVEL_param')

     ! inquire the length of north-south dimension
     ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
     call LIS_verify(ios, 'Error in nf90_inquire_dimension in AC74_read_MULTILEVEL_param')

     ! inquire the ID of parameter.
     ios = nf90_inq_varid(nid, Trim(ncvar_name), param_ID)
     call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

     ! inquire the IDs of all dimensions. The third dimension is the level dimension
     ios = nf90_inquire_variable(nid, param_ID, dimids = dimids)
     call LIS_verify(ios, trim(ncvar_name)//' failed to inquire dimensions')

     ! inquire the length of the level dimension
     ios = nf90_inquire_dimension(nid, dimids(3), len=nlevel)
     call LIS_verify(ios, trim(ncvar_name)//' failed to inquire the length of the 3rd dimension')

     ! allocate memory
     allocate(level_data (LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

     ! inquire the variable ID of parameter
     ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
     call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

     ! read parameter
     ios = nf90_get_var(nid, param_ID, level_data)
     call LIS_verify(ios, 'Error in nf90_get_var in AC74_read_MULTILEVEL_param')

     ! close netcdf file
     ios = nf90_close(nid)
     call LIS_verify(ios, 'Error in nf90_close in AC74_read_MULTILEVEL_param')

     ! grab parameter at specific level
     placeholder(:, :) = &
          level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
          LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

     ! free memory
     deallocate(level_data)

  else
     write(LIS_logunit, *) 'MULTILEVEL parameter data file: ', LIS_rc%paramfile(n), ' does not exist'
     write(LIS_logunit, *) 'program stopping ...'
     call LIS_endrun
  endif

#endif

end subroutine AC74_read_MULTILEVEL_param

!BOP
!
! !ROUTINE: AC74_read_croptype
!  \label{AC74_read_croptype}
!
! !REVISION HISTORY:
!  04 NOV 2024; Louise Busschaert, initial implementation
!
! !INTERFACE:
subroutine ac74_read_croptype(n)

  ! !USES
  use AC74_lsmMod
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_coreMod
  use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  ! !ARGUMENTS
  integer, intent(in)     :: n

  ! !DESCRIPTION:
  !  This subroutine checks the crop classifcation and assigns a
  !  a crop (corresponding to a crop file .CRO) to each tile.
  !
  !  The arguments are:
  !  \begin{description}
  !   \item[n]
  !    index of n
  !   \end{description}
  !
  !EOP
  character(len=LIS_CONST_PATH_LEN) :: crop_path
  integer                 :: ftn
  integer                 :: t,col,row,IINDEX,CT,n_crops
  integer                 :: nid, ios, croptypeId
  integer                 :: rc
  logical                 :: file_exists
  character(len=100),allocatable :: croptypes(:)
  real, allocatable       :: l_croptype(:,:)
  real, allocatable       :: glb_croptype(:,:)
  character*128 :: mess
  ! __________________________________________________________________________

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  call ESMF_ConfigFindLabel(LIS_config,"AquaCrop.7.4 crop library directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config,crop_path,rc=rc)
  call LIS_verify(rc,'AquaCrop.7.4 crop library directory: not specified')

  !- Set path to crop files
  AC74_struc(n)%PathCropFiles = trim(crop_path)//"AC_Crop.Files/"

  !- Read in crop table
  ftn = LIS_getNextUnitNumber()
  open(ftn, FILE=trim(crop_path)//"AC_Crop.Inventory",FORM='FORMATTED',STATUS='OLD',IOSTAT=rc)
  call LIS_verify(rc,'AC74_setup.F: failure opening AC_Crop.Inventory')

  read(ftn,*) n_crops

  allocate(character(len=100) :: croptypes(n_crops))

  write( mess , *) 'AC_Crop.Inventory contains ', n_crops, ' types'
  if (LIS_masterproc) then
     write(LIS_logunit, *) trim(mess)
  endif
  read(ftn,*)
  do CT=1,n_crops
     read(ftn,*) IINDEX, croptypes(CT)
  enddo
  call LIS_releaseUnitNumber(ftn)

  allocate(l_croptype(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  !- Read in crop type map file (specified in LIS parameter input file)
  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then
     ios = nf90_open(path=LIS_rc%paramfile(n),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in the lis input netcdf file')

     write(LIS_logunit,*) "[INFO] Reading in the crop type field ... "

     allocate(glb_croptype(LIS_rc%gnc(n),LIS_rc%gnr(n)))

     ios = nf90_inq_varid(nid,'AC_CROPT',croptypeId)
     call LIS_verify(ios,'nf90_inq_varid failed for AC_CROPT')

     ios = nf90_get_var(nid, croptypeId, glb_croptype)
     call LIS_verify(ios,'nf90_get_var failed for AC_CROPT')

     ios = nf90_close(nid)
     call LIS_verify(ios,'nf90_close failed in AC74_setup')

     l_croptype(:,:) = glb_croptype(&
          LIS_ews_halo_ind(n,LIS_localPet+1): &
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1))
     deallocate(glb_croptype)

     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
        row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row

        AC74_struc(n)%ac74(t)%cropt = trim(croptypes(nint(l_croptype(col,row))))
     enddo
     deallocate(croptypes)

  else
     write(LIS_logunit,*) "[ERR] The croptype inventory: ",&
          LIS_rc%paramfile(n)," does not exist."
     write(LIS_logunit,*) "[ERR] Program stopping ..."
     call LIS_endrun
  endif
  deallocate(l_croptype)
#endif

end subroutine ac74_read_croptype


!BOP
!
! !ROUTINE: AC74_param_has_var
!  \label{AC74_param_has_var}
!
! !REVISION HISTORY:
!  19 OCT 2026 Louise Busschaert; Initial implementation (from AC72 GMD branch)
!
! !INTERFACE:
logical function AC74_param_has_var(n, pname)

  ! !USES:
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  !
  ! !DESCRIPTION:
  !  True when the LIS param file of nest n has the variable pname.
  !EOP

  implicit none

  integer, intent(in)          :: n
  character(len=*), intent(in) :: pname

  integer :: nid, varid, ios
  logical :: file_exists

  AC74_param_has_var = .false.

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=trim(LIS_rc%paramfile(n)), exist=file_exists)
  if (.not. file_exists) return

  ios = nf90_open(trim(LIS_rc%paramfile(n)), NF90_NOWRITE, nid)
  if (ios /= NF90_NOERR) return

  ios = nf90_inq_varid(nid, trim(pname), varid)
  if (ios == NF90_NOERR) AC74_param_has_var = .true.

  ios = nf90_close(nid)
  call LIS_verify(ios, 'Error in nf90_close in AC74_param_has_var')
#endif
end function AC74_param_has_var
