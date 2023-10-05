!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

!BOP
! 
! !ROUTINE: ac70_getirrigationstates
! \label{ac70_getirrigationstates}
! 
! !INTERFACE:
subroutine ac70_getirrigationstates(n,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac70_lsmMod
						

! !DESCRIPTION:

! Gets the irrigation from the AquaCrop structure and transfers
! it to the irrigation structure to be written out as output.
! All irrigaiton calculations are performed in AquaCrop,
! the obtained values in mm/day are converted to kg m-2 s-1
! to remain consistent with other LSM outputs.
! 
! No scaling is applied (could be later by considering the map of
! irrigated fractions from Salmon et al. (2013).
! Also no rooting deoth is considered as it is defined by the AquaCrop
! crop file (.CRO).
! 
!
! REVISION HISTORY:
!
! Aug 2008: Hiroko Kato; Initial code
! Nov 2012: Sujay Kumar, Incorporated into LIS
! Oct 2023: Louise Busschaert; Incorporated into AC70,
!                              added dynamic irrigation season option
!
!EOP
  implicit none

  integer              :: n,rc
  integer              :: t
  type(ESMF_State)     :: irrigState
  type(ESMF_Field)     :: irrigRateField
  
  real,  pointer       :: irrigRate(:)
  real				   :: gthresh !LB: dynamic irr

  call ESMF_StateGet(irrigState, "Irrigation rate",irrigRateField,rc=rc)
  call LIS_verify(rc,'ESMF_StateGet failed for Irrigation rate')    
  call ESMF_FieldGet(irrigRateField, localDE=0,farrayPtr=irrigRate,rc=rc)
  call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation rate')


!----------------------------------------------------------------------
! Extract irrigation from LSM and convert to the right units
!----------------------------------------------------------------------

    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      ! Extract irrigation from AC70 structure
        irrigRate(t) = AC70_struc(n)%ac70(t)%Irrigation/86400
        
        !---------------------------------------------------------------------
	! For next time step: dynamic irrigation (if option turned on)
	!----------------------------------------------------------------------
	if (LIS_rc%irrigation_dveg .eq. 1) then
		! Calculate threshold
		gthresh = AC70_struc(n)%ac70(t)%Crop%CCini &
			+ (LIS_rc%irrigation_GVFparam1 + LIS_rc%irrigation_GVFparam2*&
			(AC70_struc(n)%ac70(t)%Crop%CCx - AC70_struc(n)%ac70(t)%Crop%CCini)) &
			* (AC70_struc(n)%ac70(t)%Crop%CCx - AC70_struc(n)%ac70(t)%Crop%CCini)
		if (AC70_struc(n)%ac70(t)%CCiActual .ge. gthresh) then
		    ! Irrigation allowed
		    AC70_struc(n)%ac70(t)%IrriInfoRecord1%TimeInfo = &
					    int(LIS_rc%irrigation_thresh)
		else ! very large threshold to block irrigation
		    AC70_struc(n)%ac70(t)%IrriInfoRecord1%TimeInfo = 400 
		endif
	endif
    end do

  end subroutine ac70_getirrigationstates
