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
! !ROUTINE: AC74_readcrd
! \label{AC74_readcrd}
!
! !REVISION HISTORY:
!   04 NOV 2024, Louise Busschaert; initial implementation for AC74
!
! !INTERFACE:
#include "LIS_misc.h"
subroutine AC74_readcrd()

  ! !USES:
  use AC74_lsmMod, only    : AC74_struc
  use ESMF
  use LIS_coreMod, only    : LIS_rc , LIS_config
  use LIS_logMod, only     : LIS_logunit , LIS_verify, LIS_endrun
  use LIS_timeMgrMod, only : LIS_parseTimeString
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  !
  ! !DESCRIPTION:
  !
  !  This routine reads the options specific to AC74 model from
  !  the LIS configuration file.
  !
  !EOP

  implicit none

  integer      :: rc
  integer      :: n, i
  character*10 :: time
  integer :: ios, funit
  real    :: IrriVersionNr
  integer, allocatable :: nids(:)
  character*32 :: soil_scheme_name,str

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  allocate(nids(LIS_rc%nnest))

  write(LIS_logunit, *) "[INFO] Start reading LIS configuration file for AquaCrop.7.4 model"
  ! open NetCDF parameter file for reading global attributes
  do n=1,LIS_rc%nnest
     ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE,ncid=nids(n))
     call LIS_verify(ios,'Error in nf90_open in '//trim(LIS_rc%paramfile(n))//' in AC74_readcrd')
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 model timestep:", rc = rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
     call LIS_verify(rc, "AquaCrop.7.4 model timestep: not defined")
     call LIS_parseTimeString(time, AC74_struc(n)%ts)
     if (AC74_struc(n)%ts.ne.86400) then
        write(LIS_logunit, *) "Fatal error: AquaCrop.7.4 only runs with a daily time step ..."
        write(LIS_logunit, *) "Please select AquaCrop.7.4 model timestep: 1da"
        call LIS_endrun()
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 restart output interval:", rc = rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
     call LIS_verify(rc,"AquaCrop.7.4 restart output interval: not defined")
     call LIS_parseTimeString(time, AC74_struc(n)%rstInterval)
  enddo

  ! First day of year of cropping period
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 starting day of crop period:", rc = rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Crop_AnnualStartDay, rc=rc)
     call LIS_verify(rc, "AquaCrop.7.4 starting day of crop period: not defined")
  enddo


  ! First month of year of cropping period
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 starting month of crop period:", rc = rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Crop_AnnualStartMonth, rc=rc)
     call LIS_verify(rc, "AquaCrop.7.4 starting month of crop period: not defined")
  enddo

  ! Maximum length of the crop cycle (end of the cropping period)
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 maximum length of crop cycle:", rc = rc)
  do n=1, LIS_rc%nnest
     if (rc == 0) then
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%MaxCropLength, rc=rc)
     else
        write(LIS_logunit, *) '[INFO] AC74 maximum length of crop cycle: not defined, 365 days'
        AC74_struc(n)%MaxCropLength = 365
     endif
     if ((AC74_struc(n)%MaxCropLength < 1) .or. (AC74_struc(n)%MaxCropLength > 365)) then
        write(LIS_logunit, *) '[ERR] AquaCrop.7.4 maximum length of crop cycle: between 1 and 365 days'
        call LIS_endrun()
     endif
  enddo

  ! Planting criteria (optional): within a window from the planting day of
  ! the main season
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 temperature criterion for cropping start:", rc = rc)
  do n=1, LIS_rc%nnest
     AC74_struc(n)%Temp_crit = .false.
     if (rc == 0) call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Temp_crit, rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 rainfall criterion for cropping start:", rc = rc)
  do n=1, LIS_rc%nnest
     AC74_struc(n)%Rainfall_crit = .false.
     if (rc == 0) call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Rainfall_crit, rc=rc)
  enddo

  do n=1, LIS_rc%nnest
     if (AC74_struc(n)%Temp_crit) then
        write(LIS_logunit, *) '[INFO] AC74 temperature criterion for cropping start'
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 temperature criterion Tmin:", rc = rc)
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Temp_crit_tmin, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 temperature criterion Tmin: not defined")
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 temperature criterion consecutive days:", rc = rc)
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Temp_crit_days, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 temperature criterion consecutive days: not defined")
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 temperature criterion occurrence:", rc = rc)
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Temp_crit_occurrence, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 temperature criterion occurrence: not defined")
        if ((AC74_struc(n)%Temp_crit_days < 1) .or. (AC74_struc(n)%Temp_crit_occurrence < 1)) then
           write(LIS_logunit, *) '[ERR] AquaCrop.7.4 temperature criterion: days and occurrence at least 1'
           call LIS_endrun()
        endif
     endif
     if (AC74_struc(n)%Rainfall_crit) then
        write(LIS_logunit, *) '[INFO] AC74 rainfall criterion for cropping start'
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 rainfall criterion amount:", rc = rc)
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Rainfall_crit_amount, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 rainfall criterion amount: not defined")
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 rainfall criterion days:", rc = rc)
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Rainfall_crit_days, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 rainfall criterion days: not defined")
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 rainfall criterion occurrence:", rc = rc)
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%Rainfall_crit_occurrence, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 rainfall criterion occurrence: not defined")
        if ((AC74_struc(n)%Rainfall_crit_days < 1) .or. (AC74_struc(n)%Rainfall_crit_occurrence < 1)) then
           write(LIS_logunit, *) '[ERR] AquaCrop.7.4 rainfall criterion: days and occurrence at least 1'
           call LIS_endrun()
        endif
     endif
     AC74_struc(n)%crit_window = 1
     if (AC74_struc(n)%Temp_crit .or. AC74_struc(n)%Rainfall_crit) then
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 criterion search window length:", rc = rc)
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%crit_window, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 criterion search window length: not defined")
        if ((AC74_struc(n)%crit_window < 1) .or. &
             (AC74_struc(n)%crit_window + AC74_struc(n)%MaxCropLength > 366)) then
           write(LIS_logunit, *) '[ERR] AquaCrop.7.4 criterion search window length: at least 1, '// &
                'and with the maximum length of crop cycle within one year'
           call LIS_endrun()
        endif
     endif
  enddo

  ! PathNameSimul
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 input path:", rc = rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          AC74_struc(n)%PathNameSimul, rc=rc)
     call LIS_verify(rc, "AquaCrop.7.4 input path: not defined")
  enddo

  ! CO2_Filename
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 CO2_Filename:", rc = rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          AC74_struc(n)%CO2_Filename, rc=rc)
     call LIS_verify(rc, "AquaCrop.7.4 CO2_Filename: not defined")
  enddo

  ! Management_Filename
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 Management_Filename:", rc = rc)
  do n=1, LIS_rc%nnest
     if (rc == 0) then
        call ESMF_ConfigGetAttribute(LIS_config, &
             AC74_struc(n)%Management_Filename, rc=rc)
        ! change lis none to AquaCrop (None)
        if ((AC74_struc(n)%Management_Filename .eq. 'none') .or. &
             (AC74_struc(n)%Management_Filename .eq. 'None')) then
           AC74_struc(n)%Management_Filename = '(None)'
        endif
     else
        write(LIS_logunit, *)'[INFO] AC74 Management_Filename: not defined, default management'
        AC74_struc(n)%Management_Filename = '(None)'
     endif
  enddo

  ! Irrigation_Filename
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 Irrigation_Filename:", rc = rc)
  do n=1, LIS_rc%nnest
     if (rc == 0) then
        call ESMF_ConfigGetAttribute(LIS_config, &
             AC74_struc(n)%Irrigation_Filename, rc=rc)
        ! change lis none to AquaCrop (None)
        if ((AC74_struc(n)%Irrigation_Filename .eq. 'none') .or. &
             (AC74_struc(n)%Irrigation_Filename .eq. 'None')) then
           AC74_struc(n)%Irrigation_Filename = '(None)'
        endif
     else
        write(LIS_logunit, *)'[INFO] AC74 Irrigation_Filename: not defined, no irrigation'
        AC74_struc(n)%Irrigation_Filename = '(None)'
     endif
     ! From version 7.3, a generate schedule has one more header line
     ! (last irrigation day)
     AC74_struc(n)%IrriLastDayLine = 0
     if (AC74_struc(n)%Irrigation_Filename .ne. '(None)') then
        open(newunit=funit, file=trim(AC74_struc(n)%PathNameSimul) &
             //trim(AC74_struc(n)%Irrigation_Filename), status='old', &
             action='read', iostat=ios)
        call LIS_verify(ios, 'AC74 Irrigation_Filename: cannot open '// &
             trim(AC74_struc(n)%Irrigation_Filename))
        read(funit, *) ! description
        read(funit, *) IrriVersionNr
        close(funit)
        if (nint(IrriVersionNr*10) >= 73) AC74_struc(n)%IrriLastDayLine = 1
     endif
  enddo

  ! AquaCrop model soil parameter table
  call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 soil parameter table:", rc = rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%soil_tbl_name, rc=rc)
     call LIS_verify(rc, "AquaCrop.7.4 soil parameter table: not defined")
  enddo


  ! soil classification scheme
  do n=1, LIS_rc%nnest
     ios = nf90_get_att(nids(n), NF90_GLOBAL, 'SOILTEXT_SCHEME', soil_scheme_name)
     call LIS_verify(ios, 'Error in nf90_get_att: SOILTEXT_SCHEME')
     if (trim(soil_scheme_name) .eq. "STATSGO") then
        AC74_struc(n)%soil_scheme_name = "STAS"
     else
        write(LIS_logunit, *) "Fatal error: currently, only STATSGO soil scheme is supported by AquaCrop!"
        call LIS_endrun()
     endif
  enddo

  ! Read Nrsoil layers from LDT
  do n=1, LIS_rc%nnest
     ios = nf90_get_att(nids(n), NF90_GLOBAL, 'SOIL_LAYERS', AC74_struc(n)%NrSoilLayers)
     call LIS_verify(ios,'Error in nf90_get_att in ac74_readcrd')
     write(LIS_logunit,*)"[INFO] Number of soil layers from LDT: ",AC74_struc(n)%NrSoilLayers
     allocate(AC74_struc(n)%Thickness(AC74_struc(n)%NrSoilLayers))
     allocate(AC74_struc(n)%init_smc(AC74_struc(n)%NrSoilLayers))

     ! Get thickness of soil layers from LDT
     do i=1, AC74_struc(n)%NrSoilLayers
        write (str, '(i0)') i
        ios = nf90_get_att(nids(n), NF90_GLOBAL, "THICKNESS_LAYER_"//trim(str), &
             AC74_struc(n)%Thickness(i))
        call LIS_verify(ios,'Error in nf90_get_att in ac74_readcrd')
        write(LIS_logunit,*)"[INFO] AC74 THICKNESS_LAYER_"//trim(str), &
             AC74_struc(n)%Thickness(i)
     enddo
  enddo

  ! Reference year for climatology (stress functions)
  do n=1, LIS_rc%nnest
     ios = nf90_get_att(nids(n), NF90_GLOBAL, 'AC_CLIM_REF_YEAR', AC74_struc(n)%tempcli_refyr)
     call LIS_verify(ios,'Error in nf90_get_att in ac74_readcrd')
     write(LIS_logunit,*)"[INFO] Reference year for climatology from LDT: ",AC74_struc(n)%tempcli_refyr
  enddo

  do n=1, LIS_rc%nnest
     AC74_struc(n)%dt = AC74_struc(n)%ts
  enddo

  ! Note: the height of the forcings (T and q, u and v) is given by the met
  ! forcing source (AC74_setup)

  do n=1,LIS_rc%nnest
     ios = nf90_close(nids(n))
     call LIS_verify(ios,'Error in nf90_close in '//trim(LIS_rc%paramfile(n))//' in AC74_readcrd')
  enddo

  ! The following lines hard code the LDT NetCDF variable names.
  do n=1, LIS_rc%nnest
     AC74_struc(n)%LDT_ncvar_soiltype = 'AC72_SOILTYPE'
     AC74_struc(n)%LDT_ncvar_tmincli_monthly = 'AC_Tmin_clim'
     AC74_struc(n)%LDT_ncvar_tmaxcli_monthly = 'AC_Tmax_clim'
     ! not produced by LDT, added externally to the LIS param file
     AC74_struc(n)%LDT_ncvar_gdd_cco = 'GDD_CCo'
     AC74_struc(n)%LDT_ncvar_gdd_maturity = 'GDD_maturity'
     AC74_struc(n)%LDT_ncvar_gdd_senescence = 'GDD_senescence'
     AC74_struc(n)%LDT_ncvar_gdd_maxr = 'GDD_maxR'
     AC74_struc(n)%LDT_ncvar_cgc = 'CGC'
     AC74_struc(n)%LDT_ncvar_cdc = 'CDC'
     AC74_struc(n)%LDT_ncvar_brel = 'AC_Brel'
  enddo

  ! set default restart format to netcdf
  do n=1,LIS_rc%nnest
     AC74_struc(n)%rformat = "netcdf"
  enddo
  ! restart run, read restart file
  if (trim(LIS_rc%startcode) == "restart") then
     Call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 restart file:", rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%rfile, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 restart file: not defined")
     enddo

     Call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 restart file format:", rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%rformat, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 restart file format: not defined")
     enddo
     ! cold start run, read initial state variables
  else
     ! volumetric liquid soil moisture
     call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.4 initial liquid soil moistures:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC74_struc(n)%init_smc, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.4 initial liquid soil moistures: not defined")
     enddo
  endif

  deallocate(nids)

  write(LIS_logunit, *) "[INFO] Finish reading LIS configuration file for AquaCrop.7.4 model"

#endif
end subroutine AC74_readcrd
