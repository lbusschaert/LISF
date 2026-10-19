!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module AC74_lsmMod
  !BOP
  !
  ! !MODULE: AC74_lsmMod
  !
  ! !DESCRIPTION:
  !
  ! This module provides the definition of derived data type used to
  ! control the operation of AC74 model. It also provides the entry method
  ! for the initialization of AC74-specific variables. The derived
  ! data type {\tt AC74\_struc} includes the variables that specify
  ! the runtime options and other control variables as described below:
  !
  ! \begin{description}
  ! \item[rfile]
  !  name of the AC74 restart file
  ! \item[rformat]
  !  format of restart file (binary or netcdf) for AC74
  ! \item[LDT\_ncvar\_soiltype]
  !   LDT NetCDF variable name for soil type index
  ! \item[LDT\_ncvar\_gdd\_*, LDT\_ncvar\_cgc, LDT\_ncvar\_cdc]
  !   NetCDF variable names of the crop GDD parameters in the LIS param file.
  !   Not produced by LDT: added externally. Used where the crop file has
  !   'External' for these parameters
  ! \item[LDT\_ncvar\_brel]
  !   NetCDF variable name of the relative biomass (soil fertility, fraction)
  !   in the LIS param file. Not produced by LDT: added externally
  ! \item[GDDfromLDT, BrelfromLDT]
  !   the GDD parameters and the relative biomass are in the LIS param file
  ! \item[ts]
  !   AC74 model time step in second
  ! \item[count]
  !   variable to keep track of the number of timesteps before an output
  ! \item[rstInterval]
  !   restart writing interval
  ! \item[outInterval]
  !   output writing interval
  ! \item[ac74]
  !  AC74 model specific variables
  ! \item[forc\_count]
  !   counter of forcing data
  ! \item[soil\_tbl\_name]
  !   AquaCrop model soil parameter table
  ! \item[soil\_scheme\_name]
  !   soil classification scheme
  ! \item[dt]
  !   time step in seconds
  ! \item[PathNameSimul]
  !   path of simulation folder
  ! \item[PathCropFiles]
  !   path of crop files folder
  ! \item[CO2_Filename]
  !   CO2 filename (should be in simul folder)
  ! \item[Management_Filename]
  !   Management filename (should be in simul folder)
  ! \item[Irrigation_Filename]
  !   Irrigation filename (should be in simul folder)
  ! \item[IrriLastDayLine]
  !   1 if the irrigation file (v7.3+) has the last irrigation day line, else 0
  ! \item[Crop_AnnualStartDay]
  !   annual planting day of the main season
  ! \item[Crop_AnnualStartMonth]
  !   annual planting month of the main season
  ! \item[MaxCropLength]
  !   maximum length of the crop cycle (days), end of the cropping period
  ! \item[DayNrLISStart, DayNrLISEnd]
  !   AquaCrop day numbers of the start and end of the LIS run
  ! \item[Temp\_crit]
  !   planting by a temperature criterion
  ! \item[Temp\_crit\_tmin, Temp\_crit\_days, Temp\_crit\_occurrence]
  !   Tmin threshold (degC), number of consecutive days and number of
  !   occurrences of the temperature criterion
  ! \item[Rainfall\_crit]
  !   planting by a rainfall criterion
  ! \item[Rainfall\_crit\_amount, Rainfall\_crit\_days, Rainfall\_crit\_occurrence]
  !   rainfall (mm) to exceed within a number of days, and number of
  !   occurrences of the rainfall criterion
  ! \item[crit\_window]
  !   length (days) of the planting window, from the planting day of the
  !   main season; the crop is planted on its last day if not met before
  ! \item[NrSoilLayers]
  !   number of soil layers !!different from compartments!!
  ! \item[max_No_compartments]
  !   maximum number of compartments (=12)
  ! \item[Thickness]
  !   thickness of soil layers
  ! \item[forchgt\_tq]
  !   height (m) of the forcings T and q, from the met forcing source
  ! \item[forchgt\_uv]
  !   height (m) of the forcings u and v, from the met forcing source
  ! \item[tempcli_refyr]
  !   reference year for climatology record for stress functions

  !
  ! !REVISION HISTORY:
  !  04 NOV 2024, Louise Busschaert; initial implementation for AC74
  !
  ! !USES:
  use AC74_module, only: AC74dec
  use ac74_utils, only: set_fatal_handler, set_warning_handler
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
  !-------------------------------------------------------------------------
  ! PUBLIC MEMBER FUNCTIONS
  !-------------------------------------------------------------------------
  public :: AC74_ini
  !-------------------------------------------------------------------------
  ! PUBLIC TYPES
  !-------------------------------------------------------------------------
  public :: AC74_struc
  !EOP
  type, public :: AC74_type_dec
     character(len=LIS_CONST_PATH_LEN) :: rfile
     character*256      :: rformat
     !-------------------------------------------------------------------------
     ! Parameter file names
     !-------------------------------------------------------------------------
     character*128      :: LDT_ncvar_soiltype
     character*128      :: LDT_ncvar_tmincli_monthly
     character*128      :: LDT_ncvar_tmaxcli_monthly
     ! Crop GDD parameters and relative biomass of soil fertility, read from
     ! the LIS param file when there (not produced by LDT, added externally)
     character*128      :: LDT_ncvar_gdd_cco
     character*128      :: LDT_ncvar_gdd_maturity
     character*128      :: LDT_ncvar_gdd_senescence
     character*128      :: LDT_ncvar_gdd_maxr
     character*128      :: LDT_ncvar_cgc
     character*128      :: LDT_ncvar_cdc
     character*128      :: LDT_ncvar_brel
     logical            :: GDDfromLDT
     logical            :: BrelfromLDT
     !-------------------------------------------------------------------------
     ! ts, Count, rstInterval, outInterval
     !-------------------------------------------------------------------------
     real               :: ts
     integer            :: count
     real               :: rstInterval
     integer            :: outInterval
     integer            :: forc_count
     !-------------------------------------------------------------------------
     ! Initial Model State for cold start
     !-------------------------------------------------------------------------
     real, pointer      :: init_smc(:)

     !-------------------------------------------------------------------------
     ! Constant Parameter
     !-------------------------------------------------------------------------
     character(len=LIS_CONST_PATH_LEN) :: soil_tbl_name
     character(len=256) :: soil_scheme_name
     real               :: dt

     character(len=LIS_CONST_PATH_LEN) :: PathNameOutp
     character(len=LIS_CONST_PATH_LEN) :: PathNameSimul
     character(len=LIS_CONST_PATH_LEN) :: PathNameList
     character(len=LIS_CONST_PATH_LEN) :: PathNameParam
     character(len=LIS_CONST_PATH_LEN) :: PathCropFiles
     character(len=LIS_CONST_PATH_LEN) :: CO2_Filename
     character(len=LIS_CONST_PATH_LEN) :: Management_Filename
     character(len=LIS_CONST_PATH_LEN) :: Irrigation_Filename
     integer            :: IrriLastDayLine
     integer            :: Crop_AnnualStartDay
     integer            :: Crop_AnnualStartMonth
     integer            :: MaxCropLength
     integer            :: DayNrLISStart
     integer            :: DayNrLISEnd
     logical            :: Temp_crit
     real               :: Temp_crit_tmin
     integer            :: Temp_crit_days
     integer            :: Temp_crit_occurrence
     logical            :: Rainfall_crit
     real               :: Rainfall_crit_amount
     integer            :: Rainfall_crit_days
     integer            :: Rainfall_crit_occurrence
     integer            :: crit_window
     integer            :: NrSoilLayers
     integer            :: max_No_compartments
     real, pointer      :: Thickness(:)
     integer            :: forchgt_tq
     integer            :: forchgt_uv
     integer            :: tempcli_refyr
     type(AC74dec), pointer :: ac74(:)
  end type AC74_type_dec

  type(AC74_type_dec), pointer :: AC74_struc(:)

contains

  !BOP
  !
  ! !ROUTINE: AC74_ini
  ! \label{AC74_ini}
  !
  ! !INTERFACE:
  subroutine AC74_ini()
    ! !USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
    use LIS_timeMgrMod, only: LIS_update_timestep, &
         LIS_registerAlarm

    ! !DESCRIPTION:
    !
    !  This routine creates the datatypes and allocates memory for AC74-specific
    !  variables. It also invokes the routine to read the runtime specific options
    !  for AC74 from the configuration file.
    !
    !  The routines invoked are:
    !  \begin{description}
    !   \item[AC74\_readcrd](\ref{AC74_readcrd}) \newline
    !    reads the runtime options for AC74 model
    !  \end{description}
    !EOP
    implicit none

    integer  :: n, t

    external :: AC74_readcrd

    ! errors and warnings of AquaCrop go to the LIS log
    call set_fatal_handler(AC74_fatal)
    call set_warning_handler(AC74_warn)

    ! allocate memory for nest
    allocate(AC74_struc(LIS_rc%nnest))

    ! read configuation information from lis.config file
    call AC74_readcrd()

    do n=1, LIS_rc%nnest
       ! allocate memory for all tiles in current nest
       allocate(AC74_struc(n)%ac74(LIS_rc%npatch(n, LIS_rc%lsm_index)))
       !------------------------------------------------------------------------
       ! allocate memory for vector variables passed to model interfaces
       ! TODO: check the following allocation statements carefully!
       !------------------------------------------------------------------------
       AC74_struc(n)%max_No_compartments = 12 ! hard coded
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          allocate(AC74_struc(n)%ac74(t)%smc(AC74_struc(n)%max_No_compartments))
          allocate(AC74_struc(n)%ac74(t)%tmincli_monthly(12))
          allocate(AC74_struc(n)%ac74(t)%tmaxcli_monthly(12))
          ! rainfall of the last days, for the planting criterion
          if (AC74_struc(n)%Rainfall_crit) then
             allocate(AC74_struc(n)%ac74(t)%RainCritBuf(AC74_struc(n)%Rainfall_crit_days))
          endif
       enddo

       ! initialize forcing variables to zeros
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          AC74_struc(n)%ac74(t)%tair = 0.0
          AC74_struc(n)%ac74(t)%tmax = 0.0
          AC74_struc(n)%ac74(t)%tmin = 0.0
          AC74_struc(n)%ac74(t)%tdew = 0.0
          AC74_struc(n)%ac74(t)%wndspd = 0.0
          AC74_struc(n)%ac74(t)%psurf = 0.0
          AC74_struc(n)%ac74(t)%prcp = 0.0
          AC74_struc(n)%ac74(t)%eto = 0.0

          !LB: Initialize HarvestNow (new from AC7.1)
          AC74_struc(n)%ac74(t)%HarvestNow = .false.

       enddo ! end of tile (t) loop
       !------------------------------------------------------------------------
       ! Model timestep Alarm
       !------------------------------------------------------------------------
       AC74_struc(n)%forc_count = 0

       call LIS_update_timestep(LIS_rc, n, AC74_struc(n)%ts)

       call LIS_registerAlarm("AC74 model alarm",&
            AC74_struc(n)%ts, &
            AC74_struc(n)%ts)

       call LIS_registerAlarm("AC74 restart alarm", &
            AC74_struc(n)%ts,&
            AC74_struc(n)%rstInterval)
       ! Set number of soil moisture layers in surface model (corresponding to
       ! the compartments for AquaCrop)
       LIS_sfmodel_struc(n)%nsm_layers = AC74_struc(n)%max_No_compartments
       allocate(LIS_sfmodel_struc(n)%lyrthk(AC74_struc(n)%max_No_compartments))
       ! Note, the default compartment size (0.10) is stored in the
       ! SURFACEMODEL output file, but these are spatially variable and included in the
       ! input file produced by LDT
       LIS_sfmodel_struc(n)%lyrthk(:) = 0.10
       LIS_sfmodel_struc(n)%ts = AC74_struc(n)%ts
    enddo
  end subroutine AC74_ini

  !BOP
  !
  ! !ROUTINE: AC74_fatal
  ! \label{AC74_fatal}
  !
  ! !INTERFACE:
  subroutine AC74_fatal(message)
    ! !USES:
    use LIS_logMod, only : LIS_logunit, LIS_endrun
    !
    ! !DESCRIPTION:
    !  Writes an error of AquaCrop to the LIS log and stops LIS.
    !EOP
    implicit none
    character(len=*), intent(in) :: message

    write(LIS_logunit, *) '[ERR] AC74: ', trim(message)
    call LIS_endrun()
  end subroutine AC74_fatal

  !BOP
  !
  ! !ROUTINE: AC74_warn
  ! \label{AC74_warn}
  !
  ! !INTERFACE:
  subroutine AC74_warn(message)
    ! !USES:
    use LIS_logMod, only : LIS_logunit
    !
    ! !DESCRIPTION:
    !  Writes a warning of AquaCrop to the LIS log.
    !EOP
    implicit none
    character(len=*), intent(in) :: message

    write(LIS_logunit, *) '[WARN] AC74: ', trim(message)
  end subroutine AC74_warn
end module AC74_lsmMod
