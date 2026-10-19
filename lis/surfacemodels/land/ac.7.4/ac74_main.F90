!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: AC74_main
! \label{AC74_main}
!
! !REVISION HISTORY:
!   04 NOV 2024, Louise Busschaert; initial implementation
!
! !INTERFACE:
subroutine AC74_main(n)
  ! !USES:
!!! MB_AC70
  use ac74_global, only: DegreesDay, &
       GetCCiActual, &
       GetCCiprev, &
       GetCCiTopEarlySen, &
       GetClimRecord, &
       GetCompartment, &
       GetCompartment_theta, &
       GetCrop, &
       GetCrop_Day1, &
       GetCrop_DayN, &
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
       GetSimulation_FromDayNr, &
       GetSimulation_SumGDD, &
       GetSimulation_SumGDD, &
       GetSimulation_SumGDDfromDay1, &
       GetSimulation_SumGDDfromDay1, &
       GetSimulation_ToDayNr, &
       GetSimulation_ToDayNr, &
       GetSimulParam, &
       GetSimulParam_GDDMethod, &
       GetSoil, &
       GetSoilLayer, &
       GetSumWaBal, &
       GetSurf0, &
       GetSurfaceStorage, &
       GetTact, &
       GetTactWeedInfested, &
       GetTemperatureRecord, &
       GetTmax, &
       GetTmin, &
       GetTotalSaltContent, &
       GetTotalWaterContent, &
       GetTpot, &
       GetZiAqua, &
       IrriMode_Generate, &
       IrriMode_Manual, &
       IrriMode_NoIrri, &
       ModeCycle_GDDays, &
       SetCCiActual, &
       SetCCiprev, &
       SetCCiTopEarlySen, &
       SetClimFile, &
       SetClimRecord_DataType, &
       SetClimRecord_fromd, &
       SetClimRecord_fromdaynr, &
       SetClimRecord_fromm, &
       SetClimRecord_fromstring, &
       SetClimRecord_fromy, &
       SetClimRecord_NrObs, &
       SetClimRecord_tod, &
       SetClimRecord_todaynr, &
       SetClimRecord_tom, &
       SetClimRecord_tostring, &
       SetClimRecord_toy, &
       SetCompartment, &
       SetCompartment_theta, &
       SetCrop, &
       SetCRsalt, &
       SetCRwater, &
       SetDaySubmerged, &
       SetDrain, &
       SetEact, &
       SetECdrain, &
       SetECiAqua, &
       SetECstorage, &
       Seteffectiverain, &
       SetEndSeason, &
       SetEpot, &
       SetETo, &
       SetEvapoEntireSoilSurface, &
       SetGenerateDepthMode, &
       SetGenerateTimeMode, &
       SetInfiltrated, &
       SetIniPercTAW, &
       SetIrriAfterSeason, &
       SetIrriBeforeSeason, &
       SetIrriECw, &
       SetIrriFirstDayNr, &
       SetIrriInfoLastDay, &
       SetIrrigation, &
       SetIrriMethod, &
       SetIrriMode, &
       SetManagement, &
       SetManagement_Cuttings, &
       SetMaxPlotNew, &
       SetMaxPlotTr, &
       SetNrCompartments, &
       SetOnset, &
       SetOut1Wabal, &
       SetOut2Crop, &
       SetOut3Prof, &
       SetOut4Salt, &
       SetOut5CompWC, &
       SetOut6CompEC, &
       SetOut7Clim, &
       SetOut8Irri, &
       SetOutDaily, &
       SetOutputAggregate, &
       SetPart1Mult, &
       SetPart2Eval, &
       SetPerennialPeriod, &
       SetPreDay, &
       SetRain, &
       SetRootingDepth, &
       SetRootZoneSalt, &
       SetRootZoneWC_Actual, &
       SetRootZoneWC_FC, &
       SetRootZoneWC_Leaf, &
       SetRootZoneWC_SAT, &
       SetRootZoneWC_Sen, &
       SetRootZoneWC_Thresh, &
       SetRootZoneWC_WP, &
       SetRootZoneWC_ZtopAct, &
       SetRootZoneWC_ZtopFC, &
       SetRootZoneWC_ZtopThresh, &
       SetRootZoneWC_ZtopWP, &
       SetRunoff, &
       SetSaltInfiltr, &
       SetSimulation, &
       SetSimulation_SumGDD, &
       SetSimulation_SumGDDfromDay1, &
       SetSimulParam, &
       SetSoil, &
       SetSoilLayer, &
       SetSumWaBal, &
       SetSurf0, &
       SetSurfaceStorage, &
       SetTact, &
       SetTactWeedInfested, &
       SetTmax, &
       SetTmaxTnxReference12MonthsRun, &
       SetTmin, &
       SetTminTnxReference12MonthsRun, &
       SetTnxReferenceFile, &
       SetTnxReferenceYear, &
       SetTotalSaltContent, &
       SetTotalWaterContent, &
       SetTpot, &
       SetZiAqua, &
       undef_int

  use ac74_kinds, only: int32, &
       int8

  use ac74_run, only: tDaysZmin, &
       SeasonEndDayNr, &
       tGDDZmin, &
       GetDayNrCropEnded, &
       SetDayNrCropEnded, &
       AdvanceOneTimeStep, &
       fIrri_close, &
       fIrri_open, &
       firri_read, &
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
       GetCrop_Tbase, &
       GetCrop_Tupper, &
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
       GetIrriInfoRecord1, &
       GetIrriInfoRecord1_NoMoreInfo, &
       GetIrriInfoRecord2, &
       GetIrriInfoRecord2, &
       GetIrriInterval, &
       GetLineNrEval, &
       GetNextSimFromDayNr, &
       GetNoMoreCrop, &
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
       SetalfaHI, &
       SetalfaHIAdj, &
       SetBin, &
       SetBout, &
       SetBprevSum, &
       SetCCiActualWeedInfested, &
       SetCCoTotal, &
       SetCCxCropWeedsNoSFstress, &
       SetCCxTotal, &
       SetCCxWitheredTpotNoS, &
       SetCDCTotal, &
       SetCGCref, &
       SetCO2i, &
       SetCoeffb0, &
       SetCoeffb0Salt, &
       SetCoeffb1, &
       SetCoeffb1Salt, &
       SetCoeffb2, &
       SetCoeffb2Salt, &
       SetCutInfoRecord1, &
       SetCutInfoRecord2, &
       SetDayFraction, &
       SetDayLastCut, &
       SetDayNr1Eval, &
       SetDayNrEval, &
       SetDayNri, &
       SetFracBiomassPotSF, &
       SetfWeedNoS, &
       SetGDDayFraction, &
       SetGDDayi, &
       SetGDDCDCTotal, &
       SetGDDCGCref, &
       SetGDDTadj, &
       SetGwTable, &
       SetHItimesAT, &
       SetHItimesAT1, &
       SetHItimesAT2, &
       SetHItimesBEF, &
       SetIrriInfoRecord1, &
       SetIrriInfoRecord2, &
       SetIrriInterval, &
       SetLineNrEval, &
       SetNextSimFromDayNr, &
       SetNoMoreCrop, &
       SetNoYear, &
       SetNrCut, &
       SetPlotVarCrop, &
       SetPreviousBmob, &
       SetPreviousBsto, &
       SetPreviousDayNr, &
       SetPreviousStressLevel, &
       SetPreviousSum, &
       SetPreviousSumETo, &
       SetPreviousSumGDD, &
       SetScorAT1, &
       SetScorAT2, &
       SetStageCode, &
       SetStartMode, &
       SetStressLeaf, &
       SetStressSenescence, &
       SetStressSFadjNEW, &
       SetStressTot, &
       SetSumETo, &
       SetSumGDD, &
       SetSumGDDcuts, &
       SetSumGDDPrev, &
       SetSumInterval, &
       SetSumKci, &
       SetSumKcTop, &
       SetSumKcTopStress, &
       SetTadj, &
       SetTimeSenescence, &
       SetTransfer, &
       SetWaterTableInProfile, &
       SetWeedRCi, &
       SetYprevSum, &
       SetZeval, &
       SetZiprev


  use ac74_tempprocessing, only: CreateTnxReference365Days

  use AC74_lsmMod
  use ac74_tools
  use LIS_constantsMod, only: LIS_CONST_TKFRZ
  use LIS_coreMod, only: LIS_rc, LIS_surface, LIS_domain
  use LIS_histDataMod
  use LIS_logMod, only     : LIS_logunit, LIS_endrun
  use LIS_timeMgrMod, only : LIS_isAlarmRinging

  implicit none

  ! !ARGUMENTS:
  integer, intent(in)  :: n

  integer              :: t
  integer              :: i
  real                 :: dt
  real                 :: lat, lon
  real                 :: tmp_elev
  integer              :: row, col
  logical              :: alarmCheck
  integer              :: l
  integer              :: irr_record_flag, DNr ! for irri file management
  character(250)       :: TempStr

  real                 :: tmp_pres, tmp_precip, tmp_tmax, tmp_tmin   ! Weather Forcing
  real                 :: tmp_tdew, tmp_swrad, tmp_wind, tmp_eto     ! Weather Forcing

  ! For AdvanceOneTimeStep
  real                 :: tmp_wpi
  !
  ! !DESCRIPTION:
  !  This is the entry point for calling the AC74 physics.
  !  This routine calls the {\AdvanceOneTimeStep} routine that runs AquaCrop

  !  The arguments are:
  !  \begin{description}
  !  \item[n]
  !    index of the nest
  !  \end{description}
  !EOP

  ! define variables for AC74
  ! check AC74 alarm. If alarm is ring, run model.
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "AC74 model alarm")
  if (alarmCheck) Then
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        dt = LIS_rc%ts
        row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
        col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
        lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
        lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
        tmp_elev = LIS_domain(n)%tile(t)%elev

        !!------ This Block Is Where We Obtain Weather Forcing ------------------------------!!
        ! retrieve forcing data from AC74_struc(n)%ac74(t) and assign to local variables

        ! PRES: Daily average surface pressure (kPa)
        tmp_pres      = (AC74_struc(n)%ac74(t)%psurf / AC74_struc(n)%forc_count) / 1000

        ! PRECIP: Total daily precipitation (rain+snow) (mm)
        tmp_precip    = (AC74_struc(n)%ac74(t)%prcp / AC74_struc(n)%forc_count) * 3600. * 24. !Convert from kg/ms2 to mm/d

        ! TMAX: maximum daily air temperature (degC)
        tmp_tmax      = AC74_struc(n)%ac74(t)%tmax - LIS_CONST_TKFRZ !Convert from K to C

        ! TMIN: minimum daily air temperature (degC)
        tmp_tmin      = AC74_struc(n)%ac74(t)%tmin - LIS_CONST_TKFRZ !Convert from K to C

        ! TDEW: average daily dewpoint temperature (degC)
        tmp_tdew      = (AC74_struc(n)%ac74(t)%tdew / AC74_struc(n)%forc_count) - LIS_CONST_TKFRZ !Convert from K to C

        ! SW_RAD: daily total incoming solar radiation (MJ/(m2d))
        tmp_swrad     = (AC74_struc(n)%ac74(t)%swdown / AC74_struc(n)%forc_count) * 0.0864 !Convert from W/m2 to MJ/(m2d)

        ! Wind: daily average wind speed (m/s)
        tmp_wind      = AC74_struc(n)%ac74(t)%wndspd / AC74_struc(n)%forc_count

        ! check validity of PRES
        if(tmp_pres .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable PRES in AC74"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif

        ! check validity of PRECIP
        if(tmp_precip .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable PRECIP in AC74"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif

        ! check validity of TMAX
        if(tmp_tmax .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable TMAX in AC74"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif

        ! check validity of TMIN
        if(tmp_tmin .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable TMIN in AC74"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif

        ! check validity of TDEW
        if(tmp_tdew .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable TDEW in AC74"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif

        ! check validity of SW_RAD
        if(tmp_swrad .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable SW_RAD in AC74"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif

        ! check validity of WIND
        if(tmp_wind .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable WIND in AC74"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif

        ! Call ETo subroutine
        call ac74_ETo_calc(tmp_pres, tmp_tmax, tmp_tmin, tmp_tdew, &
             tmp_wind, tmp_swrad, &
             tmp_elev, lat, tmp_eto)
        AC74_struc(n)%ac74(t)%eto = tmp_eto

        ! setting all global variables
        call SetRootZoneWC_Actual(AC74_struc(n)%ac74(t)%RootZoneWC_Actual)
        call SetRootZoneWC_FC(AC74_struc(n)%ac74(t)%RootZoneWC_FC)
        call SetRootZoneWC_WP(AC74_struc(n)%ac74(t)%RootZoneWC_WP)
        call SetRootZoneWC_SAT(AC74_struc(n)%ac74(t)%RootZoneWC_SAT)
        call SetRootZoneWC_Leaf(AC74_struc(n)%ac74(t)%RootZoneWC_Leaf)
        call SetRootZoneWC_Thresh(AC74_struc(n)%ac74(t)%RootZoneWC_Thresh)
        call SetRootZoneWC_Sen(AC74_struc(n)%ac74(t)%RootZoneWC_Sen)
        call SetRootZoneWC_ZtopAct(AC74_struc(n)%ac74(t)%RootZoneWC_ZtopAct)
        call SetRootZoneWC_ZtopFC(AC74_struc(n)%ac74(t)%RootZoneWC_ZtopFC)
        call SetRootZoneWC_ZtopWP(AC74_struc(n)%ac74(t)%RootZoneWC_ZtopWP)
        call SetRootZoneWC_ZtopThresh(AC74_struc(n)%ac74(t)%RootZoneWC_ZtopThresh)
        call SetCompartment(AC74_struc(n)%ac74(t)%Compartment)
        call SetTotalSaltContent(AC74_struc(n)%ac74(t)%TotalSaltContent)
        call SetTotalWaterContent(AC74_struc(n)%ac74(t)%TotalWaterContent)
        call Seteffectiverain(AC74_struc(n)%ac74(t)%effectiverain)
        call SetSumWaBal(AC74_struc(n)%ac74(t)%SumWaBal)
        call SetRootZoneSalt(AC74_struc(n)%ac74(t)%RootZoneSalt)
        call SetSimulation(AC74_struc(n)%ac74(t)%Simulation)
        call SetIrriInterval(AC74_struc(n)%ac74(t)%IrriInterval)
        call SetIrriInfoRecord1(AC74_struc(n)%ac74(t)%IrriInfoRecord1)
        call SetIrriInfoRecord2(AC74_struc(n)%ac74(t)%IrriInfoRecord2)
        call SetIrrigation(AC74_struc(n)%ac74(t)%Irrigation)
        do l=1, AC74_struc(n)%ac74(t)%NrCompartments
           call SetCompartment_theta(l,AC74_struc(n)%ac74(t)%smc(l))
        enddo
        call SetIrriECw(AC74_struc(n)%ac74(t)%IrriECw)
        call SetManagement(AC74_struc(n)%ac74(t)%Management)
        call SetPerennialPeriod(AC74_struc(n)%ac74(t)%PerennialPeriod)
        call SetSimulParam(AC74_struc(n)%ac74(t)%simulparam)
        call SetManagement_Cuttings(AC74_struc(n)%ac74(t)%Cuttings)
        call SetOnset(AC74_struc(n)%ac74(t)%onset)
        call SetEndSeason(AC74_struc(n)%ac74(t)%endseason)
        call SetCrop(AC74_struc(n)%ac74(t)%crop)
        call SetSoil(AC74_struc(n)%ac74(t)%Soil)
        call SetIrriBeforeSeason(AC74_struc(n)%ac74(t)%IrriBeforeSeason)
        call SetIrriAfterSeason(AC74_struc(n)%ac74(t)%IrriAfterSeason)
        call SetSoilLayer(AC74_struc(n)%ac74(t)%soillayer)
        call SetDayNri(AC74_struc(n)%ac74(t)%daynri)

        call SetGenerateTimeMode(AC74_struc(n)%ac74(t)%GenerateTimeMode)
        call SetGenerateDepthMode(AC74_struc(n)%ac74(t)%GenerateDepthMode)
        call SetIrriMode(AC74_struc(n)%ac74(t)%IrriMode)
        call SetDaySubmerged(AC74_struc(n)%ac74(t)%DaySubmerged)
        call SetMaxPlotNew(AC74_struc(n)%ac74(t)%MaxPlotNew)
        call SetNrCompartments(AC74_struc(n)%ac74(t)%NrCompartments)
        call SetIrriFirstDayNr(AC74_struc(n)%ac74(t)%IrriFirstDayNr)
        call SetZiAqua(AC74_struc(n)%ac74(t)%ZiAqua)
        call SetIniPercTAW(AC74_struc(n)%ac74(t)%IniPercTAW)
        call SetMaxPlotTr(AC74_struc(n)%ac74(t)%MaxPlotTr)
        call SetOutputAggregate(AC74_struc(n)%ac74(t)%OutputAggregate)

        call SetEvapoEntireSoilSurface(AC74_struc(n)%ac74(t)%EvapoEntireSoilSurface)
        call SetPreDay(AC74_struc(n)%ac74(t)%PreDay)
        call SetOutDaily(AC74_struc(n)%ac74(t)%OutDaily)
        call SetOut1Wabal(AC74_struc(n)%ac74(t)%Out1Wabal)
        call SetOut2Crop(AC74_struc(n)%ac74(t)%Out2Crop)
        call SetOut3Prof(AC74_struc(n)%ac74(t)%Out3Prof)
        call SetOut4Salt(AC74_struc(n)%ac74(t)%Out4Salt)
        call SetOut5CompWC(AC74_struc(n)%ac74(t)%Out5CompWC)
        call SetOut6CompEC(AC74_struc(n)%ac74(t)%Out6CompEC)
        call SetOut7Clim(AC74_struc(n)%ac74(t)%Out7Clim)
        call SetOut8Irri(AC74_struc(n)%ac74(t)%Out8Irri)
        call SetPart1Mult(AC74_struc(n)%ac74(t)%Part1Mult)
        call SetPart2Eval(AC74_struc(n)%ac74(t)%Part2Eval)

        !
        call SetCCiActual(AC74_struc(n)%ac74(t)%CCiActual)
        call SetCCiprev(AC74_struc(n)%ac74(t)%CCiprev)

        call SetCCiTopEarlySen(AC74_struc(n)%ac74(t)%CCiTopEarlySen)
        call SetCRsalt(AC74_struc(n)%ac74(t)%CRsalt)
        call SetCRwater(AC74_struc(n)%ac74(t)%CRwater)
        call SetECdrain(AC74_struc(n)%ac74(t)%ECdrain)
        call SetEciAqua(AC74_struc(n)%ac74(t)%ECiAqua)
        call SetECstorage(AC74_struc(n)%ac74(t)%ECstorage)
        call SetEact(AC74_struc(n)%ac74(t)%Eact)
        call SetEpot(AC74_struc(n)%ac74(t)%Epot)
        call SetDrain(AC74_struc(n)%ac74(t)%Drain)
        call SetInfiltrated(AC74_struc(n)%ac74(t)%Infiltrated)
        call SetRootingDepth(AC74_struc(n)%ac74(t)%RootingDepth)
        call SetRunoff(AC74_struc(n)%ac74(t)%Runoff)
        call SetSaltInfiltr(AC74_struc(n)%ac74(t)%SaltInfiltr)
        call SetSurf0(AC74_struc(n)%ac74(t)%Surf0)
        call SetSurfaceStorage(AC74_struc(n)%ac74(t)%SurfaceStorage)
        call SetTact(AC74_struc(n)%ac74(t)%Tact)
        call SetTpot(AC74_struc(n)%ac74(t)%Tpot)
        call SetTactWeedInfested(AC74_struc(n)%ac74(t)%TactWeedInfested)

        call SetGwTable(AC74_struc(n)%ac74(t)%GwTable)
        call SetPlotVarCrop(AC74_struc(n)%ac74(t)%PlotVarCrop)
        call SetStressTot(AC74_struc(n)%ac74(t)%StressTot)
        call SetCutInfoRecord1(AC74_struc(n)%ac74(t)%CutInfoRecord1)
        call SetCutInfoRecord2(AC74_struc(n)%ac74(t)%CutInfoRecord2)
        call SetTransfer(AC74_struc(n)%ac74(t)%Transfer)
        call SetPreviousSum(AC74_struc(n)%ac74(t)%PreviousSum)
        call SetTadj(AC74_struc(n)%ac74(t)%Tadj)
        call SetGDDTadj(AC74_struc(n)%ac74(t)%GDDTadj)
        call SetDayLastCut(AC74_struc(n)%ac74(t)%DayLastCut)
        call SetNrCut(AC74_struc(n)%ac74(t)%NrCut)
        call SetSumInterval(AC74_struc(n)%ac74(t)%SumInterval)
        call SetPreviousStressLevel(AC74_struc(n)%ac74(t)%PreviousStressLevel)
        call SetStressSFadjNEW(AC74_struc(n)%ac74(t)%StressSFadjNEW)
        call SetBin(AC74_struc(n)%ac74(t)%Bin)
        call SetBout(AC74_struc(n)%ac74(t)%Bout)
        call SetCO2i(AC74_struc(n)%ac74(t)%CO2i)
        call SetFracBiomassPotSF(AC74_struc(n)%ac74(t)%FracBiomassPotSF)
        call SetSumETo(AC74_struc(n)%ac74(t)%SumETo)
        call SetSumGDD(AC74_struc(n)%ac74(t)%SumGDD)
        call SetZiprev(AC74_struc(n)%ac74(t)%Ziprev)
        call SetSumGDDPrev(AC74_struc(n)%ac74(t)%SumGDDPrev)
        call SetCCxWitheredTpotNoS(AC74_struc(n)%ac74(t)%CCxWitheredTpotNoS)
        call SetCoeffb0(AC74_struc(n)%ac74(t)%Coeffb0)
        call SetCoeffb1(AC74_struc(n)%ac74(t)%Coeffb1)
        call SetCoeffb2(AC74_struc(n)%ac74(t)%Coeffb2)
        call SetCoeffb0Salt(AC74_struc(n)%ac74(t)%Coeffb0Salt)
        call SetCoeffb1Salt(AC74_struc(n)%ac74(t)%Coeffb1Salt)
        call SetCoeffb2Salt(AC74_struc(n)%ac74(t)%Coeffb2Salt)
        call SetStressLeaf(AC74_struc(n)%ac74(t)%StressLeaf)
        call SetStressSenescence(AC74_struc(n)%ac74(t)%StressSenescence)
        call SetDayFraction(AC74_struc(n)%ac74(t)%DayFraction)
        call SetGDDayFraction(AC74_struc(n)%ac74(t)%GDDayFraction)
        call SetCGCref(AC74_struc(n)%ac74(t)%CGCref)
        call SetGDDCGCref(AC74_struc(n)%ac74(t)%GDDCGCref)
        call SetTimeSenescence(AC74_struc(n)%ac74(t)%TimeSenescence)
        call SetSumKcTop(AC74_struc(n)%ac74(t)%SumKcTop)
        call SetSumKcTopStress(AC74_struc(n)%ac74(t)%SumKcTopStress)
        call SetSumKci(AC74_struc(n)%ac74(t)%SumKci)
        call SetCCoTotal(AC74_struc(n)%ac74(t)%CCoTotal)
        call SetCCxTotal(AC74_struc(n)%ac74(t)%CCxTotal)
        call SetCDCTotal(AC74_struc(n)%ac74(t)%CDCTotal)
        call SetGDDCDCTotal(AC74_struc(n)%ac74(t)%GDDCDCTotal)
        call SetCCxCropWeedsNoSFstress(AC74_struc(n)%ac74(t)%CCxCropWeedsNoSFstress)
        call SetWeedRCi(AC74_struc(n)%ac74(t)%WeedRCi)
        call SetCCiActualWeedInfested(AC74_struc(n)%ac74(t)%CCiActualWeedInfested)
        call SetfWeedNoS(AC74_struc(n)%ac74(t)%fWeedNoS)
        call SetZeval(AC74_struc(n)%ac74(t)%Zeval)
        call SetBprevSum(AC74_struc(n)%ac74(t)%BprevSum)
        call SetYprevSum(AC74_struc(n)%ac74(t)%YprevSum)
        call SetSumGDDcuts(AC74_struc(n)%ac74(t)%SumGDDcuts)
        call SetHItimesBEF(AC74_struc(n)%ac74(t)%HItimesBEF)
        call SetScorAT1(AC74_struc(n)%ac74(t)%ScorAT1)
        call SetScorAT2(AC74_struc(n)%ac74(t)%ScorAT2)
        call SetHItimesAT1(AC74_struc(n)%ac74(t)%HItimesAT1)
        call SetHItimesAT2(AC74_struc(n)%ac74(t)%HItimesAT2)
        call SetHItimesAT(AC74_struc(n)%ac74(t)%HItimesAT)
        call SetalfaHI(AC74_struc(n)%ac74(t)%alfaHI)
        call SetalfaHIAdj(AC74_struc(n)%ac74(t)%alfaHIAdj)
        call SetNextSimFromDayNr(AC74_struc(n)%ac74(t)%NextSimFromDayNr)
        call SetDayNr1Eval(AC74_struc(n)%ac74(t)%DayNr1Eval)
        call SetDayNrEval(AC74_struc(n)%ac74(t)%DayNrEval)
        call SetLineNrEval(int(AC74_struc(n)%ac74(t)%LineNrEval,kind=int32))
        call SetPreviousSumETo(AC74_struc(n)%ac74(t)%PreviousSumETo)
        call SetPreviousSumGDD(AC74_struc(n)%ac74(t)%PreviousSumGDD)
        call SetPreviousBmob(AC74_struc(n)%ac74(t)%PreviousBmob)
        call SetPreviousBsto(AC74_struc(n)%ac74(t)%PreviousBsto)
        call SetStageCode(AC74_struc(n)%ac74(t)%StageCode)
        call SetPreviousDayNr(AC74_struc(n)%ac74(t)%PreviousDayNr)
        call SetNoYear(AC74_struc(n)%ac74(t)%NoYear)
        call SetWaterTableInProfile(AC74_struc(n)%ac74(t)%WaterTableInProfile)
        call SetStartMode(AC74_struc(n)%ac74(t)%StartMode)
        call SetGDDayi(AC74_struc(n)%ac74(t)%GDDayi)
        call SetNoMoreCrop(AC74_struc(n)%AC74(t)%NoMoreCrop)
        call SetDayNrCropEnded(AC74_struc(n)%ac74(t)%DayNrCropEnded)
        call SetIrriInfoLastDay(AC74_struc(n)%ac74(t)%IrriInfoLastDay)
        tDaysZmin = AC74_struc(n)%ac74(t)%tDaysZmin
        tGDDZmin = AC74_struc(n)%ac74(t)%tGDDZmin

        ! Fixed var
        call SetOut3Prof(.true.) ! needed for correct rootzone sm
        call SetOutDaily(.true.)

        ! Tmin and Tmax reference (climatology) of the tile
        call SetTnxReferenceFile('(External)')
        call SetTminTnxReference12MonthsRun(AC74_struc(n)%ac74(t)%tmincli_monthly(:))
        call SetTmaxTnxReference12MonthsRun(AC74_struc(n)%ac74(t)%tmaxcli_monthly(:))
        ! Daily reference, needed to stop irrigating the last days of a GDD crop
        if ((AC74_struc(n)%ac74(t)%IrriMode.eq.IrriMode_Generate) &
             .and. (AC74_struc(n)%ac74(t)%IrriInfoLastDay.ne.undef_int) &
             .and. (AC74_struc(n)%ac74(t)%crop%ModeCycle.eq.ModeCycle_GDDays)) then
           call CreateTnxReference365Days()
        endif

        ! Set climate variables
        ! Round them to 4 digits after the comma as done in the AC standalone
        tmp_precip = anint(tmp_precip*10000)/10000
        call SetRain(tmp_precip)
        tmp_tmin = anint(tmp_tmin*10000)/10000
        call SetTmin(tmp_tmin)
        tmp_tmax = anint(tmp_tmax*10000)/10000
        call SetTmax(tmp_tmax)
        tmp_eto = anint(tmp_eto*10000)/10000
        call SetETo(tmp_eto)

        ! Planting criteria with the forcing of today: when met, a new run
        ! starts today (initialized below)
        call AC74_plant_search(n, t, GetDayNri(), tmp_tmin, tmp_precip)

        ! SumGDD calculation needed only for second day when not done within InitializeSimulationRunPart2
        if (GetDayNri()>GetSimulation_FromDayNr()) then
           ! Sum of GDD at end of first day
           call SetGDDayi(DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), GetTmin(), &
                GetTmax(), GetSimulParam_GDDMethod()))
           if (GetDayNri() >= GetCrop_Day1()) then
              call SetSimulation_SumGDD(GetSimulation_SumGDD() + GetGDDayi())
              call SetSimulation_SumGDDfromDay1(GetSimulation_SumGDDfromDay1() + &
                   GetGDDayi())
           end if
        end if

        ! Start irrigation block
        call SetIrriMode(AC74_struc(n)%ac74(t)%IrriMode)
        irr_record_flag = 0
        if(AC74_struc(n)%ac74(t)%IrriMode.ne.IrriMode_NoIrri) then
           ! Irrigation ON for tile
           call SetIrriAfterSeason(AC74_struc(n)%ac74(t)%IrriAfterSeason)
           call SetIrriBeforeSeason(AC74_struc(n)%ac74(t)%IrriBeforeSeason)
           call SetIrriECw(AC74_struc(n)%ac74(t)%IrriECw)
           call SetIrriInfoRecord1(AC74_struc(n)%ac74(t)%IrriInfoRecord1)
           call SetIrriInfoRecord2(AC74_struc(n)%ac74(t)%IrriInfoRecord2)
           call SetIrriMethod(AC74_struc(n)%ac74(t)%IrriMethod)
           call SetGenerateDepthMode(AC74_struc(n)%ac74(t)%GenerateDepthMode)
           call SetGenerateTimeMode(AC74_struc(n)%ac74(t)%GenerateTimeMode)

           ! For IrriMode_Generate and IrriMode_Manual
           ! Check if new record needs to be read
           if((AC74_struc(n)%ac74(t)%IrriMode.eq.IrriMode_Generate) &
                .and. (GetDayNri() >= GetCrop_Day1()) &
                .and. (GetDayNri() <= GetCrop_DayN())) then
              if((AC74_struc(n)%ac74(t)%daynri-AC74_struc(n)%ac74(t)%Crop%Day1+1)&
                   .gt.AC74_struc(n)%ac74(t)%IrriInfoRecord1%ToDay) then
                 irr_record_flag = 1
              endif
           elseif(AC74_struc(n)%ac74(t)%IrriMode.eq.IrriMode_Manual) then
              ! Check start date of schedule
              if(AC74_struc(n)%ac74(t)%IrriFirstDayNr.eq.undef_int)then
                 DNr = AC74_struc(n)%ac74(t)%daynri &
                      - AC74_struc(n)%ac74(t)%Crop%Day1 + 1
              else
                 DNr = AC74_struc(n)%ac74(t)%daynri &
                      - AC74_struc(n)%ac74(t)%IrriFirstDayNr + 1
              endif
              if(AC74_struc(n)%ac74(t)%IrriInfoRecord1%TimeInfo.eq.DNr)then
                 irr_record_flag = 1
              endif
           endif

           ! re-open irrigation file and read the previous records
           if(irr_record_flag.eq.1)then
              call fIrri_open(trim(AC74_struc(n)%PathNameSimul)&
                   //trim(AC74_struc(n)%Irrigation_Filename), 'r')
              do i=1,AC74_struc(n)%ac74(t)%irri_lnr
                 TempStr = fIrri_read()
              enddo
              AC74_struc(n)%ac74(t)%irri_lnr = AC74_struc(n)%ac74(t)%irri_lnr + 1
           endif
        endif
        ! End irrigation block

!!! initialize the planned run of the tile (start of LIS, end of a season, or planting)
        if (AC74_struc(n)%ac74(t)%InitializeRun.eq.1) then
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

           ! Planned run of the tile in the project run
           ! (crop parameters are read when calling InitializeRunPart1)
           call AC74_set_run_input(n, t)

           ! Set reference year for CO2 for stress functions
           ! (Tmin and Tmax reference are set above)
           call SetTnxReferenceYear(AC74_struc(n)%tempcli_refyr)

           ! InitializeRunPart (the single project run)
           call InitializeRunPart1(1_int8, AC74_struc(n)%ac74(t)%TheProjectType)
           call InitializeSimulationRunPart2()
           AC74_struc(n)%ac74(t)%HarvestNow = .false. ! Initialize to false

           ! Overwrite the SMC to avoid problems when restarting
           do l=1, AC74_struc(n)%ac74(t)%NrCompartments
              call SetCompartment_theta(l,AC74_struc(n)%ac74(t)%smc(l))
           enddo

           ! Irrigaton file management after InitializeRun
           if(GetIrriMode().ne.IrriMode_NoIrri) then
              call fIrri_close()
              ! Check for irrigation (irrigation file management)
              if(GetIrriMode().eq.IrriMode_Manual)then
                 if(GetIrriInfoRecord1_NoMoreInfo())then
                    AC74_struc(n)%ac74(t)%irri_lnr = 9
                 else
                    AC74_struc(n)%ac74(t)%irri_lnr = 10
                 endif
              elseif(GetIrriMode().eq.IrriMode_Generate)then
                 if(AC74_struc(n)%ac74(t)%IrriInfoRecord1%NoMoreInfo)then
                    AC74_struc(n)%ac74(t)%irri_lnr = 11 + AC74_struc(n)%IrriLastDayLine
                 else
                    AC74_struc(n)%ac74(t)%irri_lnr = 12 + AC74_struc(n)%IrriLastDayLine
                 endif
              else ! no irrigation, set to 0
                 AC74_struc(n)%ac74(t)%irri_lnr = 0
              endif
           endif
           ! End irrigation block
           AC74_struc(n)%ac74(t)%InitializeRun = 0 ! Initialization done
        end if

        ! Run AC
        tmp_wpi = AC74_struc(n)%ac74(t)%WPi
        call AdvanceOneTimeStep(tmp_wpi, AC74_struc(n)%ac74(t)%HarvestNow)
        AC74_struc(n)%ac74(t)%WPi = tmp_wpi

        ! Close irri file if opened
        if(irr_record_flag.eq.1)then
           call fIrri_close()
        endif

        ! Get all the ac74 variables and store in AC74_struc
        do l=1, AC74_struc(n)%ac74(t)%NrCompartments
           AC74_struc(n)%ac74(t)%smc(l) = GetCompartment_theta(l)
        enddo

        AC74_struc(n)%ac74(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
        AC74_struc(n)%ac74(t)%RootZoneWC_FC = GetRootZoneWC_FC()
        AC74_struc(n)%ac74(t)%RootZoneWC_WP = GetRootZoneWC_WP()
        AC74_struc(n)%ac74(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
        AC74_struc(n)%ac74(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
        AC74_struc(n)%ac74(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
        AC74_struc(n)%ac74(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
        AC74_struc(n)%ac74(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
        AC74_struc(n)%ac74(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
        AC74_struc(n)%ac74(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
        AC74_struc(n)%ac74(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
        AC74_struc(n)%ac74(t)%Compartment = GetCompartment()
        AC74_struc(n)%ac74(t)%TotalSaltContent = GetTotalSaltContent()
        AC74_struc(n)%ac74(t)%TotalWaterContent = GetTotalWaterContent()
        AC74_struc(n)%ac74(t)%effectiverain = Geteffectiverain()
        AC74_struc(n)%ac74(t)%SumWaBal = GetSumWaBal()
        AC74_struc(n)%ac74(t)%RootZoneSalt = GetRootZoneSalt()
        AC74_struc(n)%ac74(t)%Simulation = GetSimulation()
        AC74_struc(n)%ac74(t)%IrriInterval = GetIrriInterval()
        AC74_struc(n)%ac74(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
        AC74_struc(n)%ac74(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
        AC74_struc(n)%ac74(t)%Irrigation = GetIrrigation()
        AC74_struc(n)%ac74(t)%IrriBeforeSeason = GetIrriBeforeSeason()
        AC74_struc(n)%ac74(t)%IrriAfterSeason = GetIrriAfterSeason()
        AC74_struc(n)%ac74(t)%SoilLayer = GetSoilLayer()
        AC74_struc(n)%ac74(t)%daynri = GetDayNri()
        do l=1, AC74_struc(n)%ac74(t)%NrCompartments
           AC74_struc(n)%ac74(t)%smc(l) = GetCompartment_theta(l)
        enddo
        AC74_struc(n)%ac74(t)%IrriECw = GetIrriECw()
        AC74_struc(n)%ac74(t)%Management = GetManagement()
        AC74_struc(n)%ac74(t)%PerennialPeriod = GetPerennialPeriod()
        AC74_struc(n)%ac74(t)%simulparam = GetSimulParam()
        AC74_struc(n)%ac74(t)%Cuttings = GetManagement_Cuttings()
        AC74_struc(n)%ac74(t)%onset = GetOnset()
        AC74_struc(n)%ac74(t)%endseason = GetEndSeason()
        AC74_struc(n)%ac74(t)%crop = GetCrop()
        AC74_struc(n)%ac74(t)%Soil = GetSoil()
        AC74_struc(n)%ac74(t)%TemperatureRecord = GetTemperatureRecord()
        AC74_struc(n)%ac74(t)%ClimRecord = GetClimRecord()
        AC74_struc(n)%ac74(t)%RainRecord = GetRainRecord()
        AC74_struc(n)%ac74(t)%EToRecord = GetEToRecord()

        AC74_struc(n)%ac74(t)%GenerateTimeMode = GetGenerateTimeMode()
        AC74_struc(n)%ac74(t)%GenerateDepthMode = GetGenerateDepthMode()
        AC74_struc(n)%ac74(t)%IrriMode = GetIrriMode()
        AC74_struc(n)%ac74(t)%IrriMethod = GetIrriMethod()
        AC74_struc(n)%ac74(t)%DaySubmerged = GetDaySubmerged()
        AC74_struc(n)%ac74(t)%MaxPlotNew = GetMaxPlotNew()
        AC74_struc(n)%ac74(t)%NrCompartments = GetNrCompartments()
        AC74_struc(n)%ac74(t)%IrriFirstDayNr = GetIrriFirstDayNr()
        AC74_struc(n)%ac74(t)%ZiAqua = GetZiAqua()
        AC74_struc(n)%ac74(t)%IniPercTAW = GetIniPercTAW()
        AC74_struc(n)%ac74(t)%MaxPlotTr = GetMaxPlotTr()
        AC74_struc(n)%ac74(t)%OutputAggregate = GetOutputAggregate()

        AC74_struc(n)%ac74(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
        AC74_struc(n)%ac74(t)%PreDay = GetPreDay()
        AC74_struc(n)%ac74(t)%OutDaily = GetOutDaily()
        AC74_struc(n)%ac74(t)%Out1Wabal = GetOut1Wabal()
        AC74_struc(n)%ac74(t)%Out2Crop = GetOut2Crop()
        AC74_struc(n)%ac74(t)%Out3Prof = GetOut3Prof()
        AC74_struc(n)%ac74(t)%Out4Salt = GetOut4Salt()
        AC74_struc(n)%ac74(t)%Out5CompWC = GetOut5CompWC()
        AC74_struc(n)%ac74(t)%Out6CompEC = GetOut6CompEC()
        AC74_struc(n)%ac74(t)%Out7Clim = GetOut7Clim()
        AC74_struc(n)%ac74(t)%Out8Irri = GetOut8Irri()
        AC74_struc(n)%ac74(t)%Part1Mult = GetPart1Mult()
        AC74_struc(n)%ac74(t)%Part2Eval = GetPart2Eval()

        !
        AC74_struc(n)%ac74(t)%CCiActual = GetCCiActual()
        AC74_struc(n)%ac74(t)%CCiprev = GetCCiprev()
        AC74_struc(n)%ac74(t)%CCiTopEarlySen = GetCCiTopEarlySen()
        AC74_struc(n)%ac74(t)%CRsalt = GetCRsalt ()
        AC74_struc(n)%ac74(t)%CRwater = GetCRwater()
        AC74_struc(n)%ac74(t)%ECdrain = GetECdrain()
        AC74_struc(n)%ac74(t)%ECiAqua = GetECiAqua()
        AC74_struc(n)%ac74(t)%ECstorage = GetECstorage()
        AC74_struc(n)%ac74(t)%Eact = GetEact()
        AC74_struc(n)%ac74(t)%Epot = GetEpot()
        AC74_struc(n)%ac74(t)%ETo = GetETo()
        AC74_struc(n)%ac74(t)%Drain = GetDrain()
        AC74_struc(n)%ac74(t)%Infiltrated = GetInfiltrated()
        AC74_struc(n)%ac74(t)%prcp = GetRain()
        AC74_struc(n)%ac74(t)%RootingDepth = GetRootingDepth()
        AC74_struc(n)%ac74(t)%Runoff = GetRunoff()
        AC74_struc(n)%ac74(t)%SaltInfiltr = GetSaltInfiltr()
        AC74_struc(n)%ac74(t)%Surf0 = GetSurf0()
        AC74_struc(n)%ac74(t)%SurfaceStorage = GetSurfaceStorage()
        AC74_struc(n)%ac74(t)%Tact = GetTact()
        AC74_struc(n)%ac74(t)%Tpot = GetTpot()
        AC74_struc(n)%ac74(t)%TactWeedInfested = GetTactWeedInfested()
        AC74_struc(n)%ac74(t)%tmax = GetTmax()
        AC74_struc(n)%ac74(t)%tmin =GetTmin()


        AC74_struc(n)%ac74(t)%GwTable = GetGwTable()
        AC74_struc(n)%ac74(t)%PlotVarCrop = GetPlotVarCrop()
        AC74_struc(n)%ac74(t)%StressTot = GetStressTot()
        AC74_struc(n)%ac74(t)%CutInfoRecord1 = GetCutInfoRecord1()
        AC74_struc(n)%ac74(t)%CutInfoRecord2 = GetCutInfoRecord2()
        AC74_struc(n)%ac74(t)%Transfer = GetTransfer()
        AC74_struc(n)%ac74(t)%PreviousSum = GetPreviousSum()
        AC74_struc(n)%ac74(t)%Tadj = GetTadj()
        AC74_struc(n)%ac74(t)%GDDTadj = GetGDDTadj()
        AC74_struc(n)%ac74(t)%DayLastCut = GetDayLastCut()
        AC74_struc(n)%ac74(t)%NrCut = GetNrCut()
        AC74_struc(n)%ac74(t)%SumInterval = GetSumInterval()
        AC74_struc(n)%ac74(t)%PreviousStressLevel = GetPreviousStressLevel()
        AC74_struc(n)%ac74(t)%StressSFadjNEW = GetStressSFadjNEW()
        AC74_struc(n)%ac74(t)%Bin = GetBin()
        AC74_struc(n)%ac74(t)%Bout = GetBout()
        AC74_struc(n)%ac74(t)%GDDayi = GetGDDayi()
        AC74_struc(n)%ac74(t)%CO2i = GetCO2i()
        AC74_struc(n)%ac74(t)%FracBiomassPotSF = GetFracBiomassPotSF()
        AC74_struc(n)%ac74(t)%SumETo = GetSumETo()
        AC74_struc(n)%ac74(t)%SumGDD = GetSumGDD()
        AC74_struc(n)%ac74(t)%Ziprev = GetZiprev()
        AC74_struc(n)%ac74(t)%SumGDDPrev = GetSumGDDPrev()
        AC74_struc(n)%ac74(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
        AC74_struc(n)%ac74(t)%Coeffb0 = GetCoeffb0()
        AC74_struc(n)%ac74(t)%Coeffb1 = GetCoeffb1()
        AC74_struc(n)%ac74(t)%Coeffb2 = GetCoeffb2()
        AC74_struc(n)%ac74(t)%Coeffb0Salt = GetCoeffb0Salt()
        AC74_struc(n)%ac74(t)%Coeffb1Salt = GetCoeffb1Salt()
        AC74_struc(n)%ac74(t)%Coeffb2Salt = GetCoeffb2Salt()
        AC74_struc(n)%ac74(t)%StressLeaf = GetStressLeaf()
        AC74_struc(n)%ac74(t)%StressSenescence = GetStressSenescence()
        AC74_struc(n)%ac74(t)%DayFraction = GetDayFraction()
        AC74_struc(n)%ac74(t)%GDDayFraction = GetGDDayFraction()
        AC74_struc(n)%ac74(t)%CGCref = GetCGCref()
        AC74_struc(n)%ac74(t)%GDDCGCref = GetGDDCGCref()
        AC74_struc(n)%ac74(t)%TimeSenescence = GetTimeSenescence()
        AC74_struc(n)%ac74(t)%SumKcTop = GetSumKcTop()
        AC74_struc(n)%ac74(t)%SumKcTopStress = GetSumKcTopStress()
        AC74_struc(n)%ac74(t)%SumKci = GetSumKci()
        AC74_struc(n)%ac74(t)%CCoTotal = GetCCoTotal()
        AC74_struc(n)%ac74(t)%CCxTotal = GetCCxTotal()
        AC74_struc(n)%ac74(t)%CDCTotal = GetCDCTotal()
        AC74_struc(n)%ac74(t)%GDDCDCTotal = GetGDDCDCTotal()
        AC74_struc(n)%ac74(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
        AC74_struc(n)%ac74(t)%WeedRCi = GetWeedRCi()
        AC74_struc(n)%ac74(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
        AC74_struc(n)%ac74(t)%fWeedNoS = GetfWeedNoS()
        AC74_struc(n)%ac74(t)%Zeval = GetZeval()
        AC74_struc(n)%ac74(t)%BprevSum = GetBprevSum()
        AC74_struc(n)%ac74(t)%YprevSum = GetYprevSum()
        AC74_struc(n)%ac74(t)%SumGDDcuts = GetSumGDDcuts()
        AC74_struc(n)%ac74(t)%HItimesBEF = GetHItimesBEF()
        AC74_struc(n)%ac74(t)%ScorAT1 = GetScorAT1()
        AC74_struc(n)%ac74(t)%ScorAT2 = GetScorAT2()
        AC74_struc(n)%ac74(t)%HItimesAT1 = GetHItimesAT1()
        AC74_struc(n)%ac74(t)%HItimesAT2 = GetHItimesAT2()
        AC74_struc(n)%ac74(t)%HItimesAT = GetHItimesAT()
        AC74_struc(n)%ac74(t)%alfaHI = GetalfaHI()
        AC74_struc(n)%ac74(t)%alfaHIAdj = GetalfaHIAdj()
        AC74_struc(n)%ac74(t)%NextSimFromDayNr = GetNextSimFromDayNr ()
        AC74_struc(n)%ac74(t)%DayNr1Eval = GetDayNr1Eval()
        AC74_struc(n)%ac74(t)%DayNrEval = GetDayNrEval()
        AC74_struc(n)%ac74(t)%LineNrEval = GetLineNrEval()
        AC74_struc(n)%ac74(t)%PreviousSumETo = GetPreviousSumETo()
        AC74_struc(n)%ac74(t)%PreviousSumGDD = GetPreviousSumGDD()
        AC74_struc(n)%ac74(t)%PreviousBmob = GetPreviousBmob()
        AC74_struc(n)%ac74(t)%PreviousBsto = GetPreviousBsto()
        AC74_struc(n)%ac74(t)%StageCode = GetStageCode()
        AC74_struc(n)%ac74(t)%PreviousDayNr = GetPreviousDayNr()
        AC74_struc(n)%ac74(t)%NoYear = GetNoYear()
        AC74_struc(n)%ac74(t)%WaterTableInProfile = GetWaterTableInProfile()
        AC74_struc(n)%ac74(t)%StartMode = GetStartMode()
        AC74_struc(n)%AC74(t)%NoMoreCrop = GetNoMoreCrop()
        AC74_struc(n)%ac74(t)%DayNrCropEnded = GetDayNrCropEnded()
        AC74_struc(n)%ac74(t)%IrriInfoLastDay = GetIrriInfoLastDay()
        AC74_struc(n)%ac74(t)%tDaysZmin = tDaysZmin
        AC74_struc(n)%ac74(t)%tGDDZmin = tGDDZmin

        ! End of the season: the crop matured or died, or the end of the
        ! cropping period is reached
        ! (DayNri - 1 because DayNri is already for next day)
        ! CropActive: sown or planted (1), or fallow (0), from the planting day of the
        ! run (AquaCrop's Crop_Day1 moves to the germination when delayed)
        AC74_struc(n)%ac74(t)%CropActive = 0.
        if ((GetDayNri()-1) >= AC74_struc(n)%ac74(t)%RunCropDay1) then
           AC74_struc(n)%ac74(t)%CropActive = 1.
           if (GetNoMoreCrop() .or. ((GetDayNri()-1) >= SeasonEndDayNr()) &
                .or. ((GetDayNri()-1) >= GetSimulation_ToDayNr())) then
              ! keep the yield of the season until the next one ends
              AC74_struc(n)%ac74(t)%SeasonYield = AC74_struc(n)%ac74(t)%SumWaBal%YieldPart
              AC74_struc(n)%ac74(t)%SeasonNr = AC74_struc(n)%ac74(t)%SeasonNr + 1
              ! next run from tomorrow: fallow until the next planting
              AC74_struc(n)%ac74(t)%irun = AC74_struc(n)%ac74(t)%irun + 1
              call AC74_plan_run(n, t, GetDayNri())
              AC74_struc(n)%ac74(t)%InitializeRun = 1
           endif
        endif

        ! Diagnostic output variables
        ![ 1] output variable: smc (unit=m^3 m-3 ). ***  volumetric soil moisture
        do i=1, AC74_struc(n)%ac74(t)%NrCompartments
           call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = AC74_struc(n)%ac74(t)%smc(i), &
                vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
        end do
        ![ 2] output variable: biomass (unit=t/ha).  *** cummulative biomass
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Biomass, value = AC74_struc(n)%ac74(t)%SumWaBal%Biomass, &
             vlevel=1, unit="t ha-1", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 3] output variable: biomass (unit=mm).  *** actual rootzone water content
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_RootZoneWC_Actual, value = AC74_struc(n)%ac74(t)%RootZoneWC_Actual, &
             vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 4] output variable: RootZoneWC_WP (unit=mm).  *** rootzone water content at wilting point
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_RootZoneWC_WP, value = AC74_struc(n)%ac74(t)%RootZoneWC_WP, &
             vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 5] output variable: RootZoneWC_FC (unit=mm).  *** rootzone water content at field capacity
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_RootZoneWC_FC, value = AC74_struc(n)%ac74(t)%RootZoneWC_FC, &
             vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 6] output variable: Tact (unit=mm).  *** actual transpiration
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Tact, value = AC74_struc(n)%ac74(t)%Tact, &
             vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 7] output variable: Eact (unit=mm).  *** actual evaporation
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Eact, value = AC74_struc(n)%ac74(t)%Eact, &
             vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 8] output variable: AC74ETo (unit=mm).  *** reference evapotranspiration (Penman-Monteith)
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_ETo, value = AC74_struc(n)%ac74(t)%eto, &
             vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 9] output variable: RootingDepth (unit=m).  *** rooting depth
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_RootingDepth, value = AC74_struc(n)%ac74(t)%RootingDepth, &
             vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 10] output variable: CCiActual (unit=-).  *** canopy cover
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_CCiActual, value = AC74_struc(n)%ac74(t)%CCiActual, &
             vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 11] output variable: AC74Tmin (unit=deg C).  *** daily minimum temperature
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Tmin, value = tmp_tmin, &
             vlevel=1, unit="degC", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 12] output variable: AC74Tmax (unit=deg C).  *** daily maximum temperature
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Tmax, value = tmp_tmax, &
             vlevel=1, unit="degC", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 13] output variable: AC74Rain (unit=mm).  *** precipitation rate
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Rain, value = tmp_precip, &
             vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 13] output variable: yield (unit=t ha-1).  *** yield
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Yield, value = AC74_struc(n)%ac74(t)%SumWaBal%YieldPart, &
             vlevel=1, unit="t ha-1", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 14] output variable: irrigation (unit=mm).  *** irrigation
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Irrigation, value = AC74_struc(n)%ac74(t)%Irrigation, &
             vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 15] output variable: StExp (unit=%).  *** expansion stress
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_StExp, value = GetStressLeaf(), &
             vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 16] output variable: StSen (unit=%).  *** senescence stress
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_StSen, value = GetStressSenescence(), &
             vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 17] output variable: SeasonYield (unit=t ha-1).  *** yield of the last completed season
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_SeasonYield, value = AC74_struc(n)%ac74(t)%SeasonYield, &
             vlevel=1, unit="t/ha", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 18] output variable: SeasonNr (unit=-).  *** number of completed seasons
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_SeasonNr, value = real(AC74_struc(n)%ac74(t)%SeasonNr), &
             vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
        ![ 19] output variable: CropActive (unit=-).  *** crop sown or planted (1), or fallow (0)
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_CropActive, value = AC74_struc(n)%ac74(t)%CropActive, &
             vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

        !  Reset forcings
        AC74_struc(n)%ac74(t)%tair = 0.0
        AC74_struc(n)%ac74(t)%tmax = 0.0
        AC74_struc(n)%ac74(t)%tmin = 0.0
        AC74_struc(n)%ac74(t)%tdew = 0.0
        AC74_struc(n)%ac74(t)%wndspd = 0.0
        AC74_struc(n)%ac74(t)%psurf = 0.0
        AC74_struc(n)%ac74(t)%prcp = 0.0
        AC74_struc(n)%ac74(t)%eto = 0.0
        AC74_struc(n)%ac74(t)%swdown = 0.0
     enddo ! end of tile (t) loop
     ! reset forcing counter to be zero
     AC74_struc(n)%forc_count = 0
  endif ! end of alarmCheck loop

end subroutine AC74_main
