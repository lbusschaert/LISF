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
! !ROUTINE: AC74_finalize
! \label{AC74_finalize}
!
! !REVISION HISTORY:
!  04 NOV 2024, Louise Busschaert; initial implementation for AC74
!
! !INTERFACE:
subroutine AC74_finalize(n)

  ! !USES:
  use AC74_lsmMod
  use LIS_coreMod, only : LIS_rc

  !
  ! !DESCRIPTION:
  !
  !  This routine cleans up the allocated memory structures in AC74
  !
  !EOP

  implicit none

  integer, intent(inout) :: n

  do n=1, LIS_rc%nnest
     !deallocate all vars
     ! free memory for ac74, the data at tile level
     deallocate(AC74_struc(n)%ac74)

     ! free memory for constant parameter
     deallocate(AC74_struc(n)%Thickness)

     ! free memory for initial state variable
     deallocate(AC74_struc(n)%init_smc)
  end do ! nest loop

  deallocate(AC74_struc)

end subroutine AC74_finalize

