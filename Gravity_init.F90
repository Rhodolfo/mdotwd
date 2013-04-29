!!****if* source/physics/Gravity/GravityMain/Roche/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!  
!! SYNOPSIS
!!
!!  Gravity_init()
!!
!! DESCRIPTION
!!
!!  This routine initializes the gravitational physics unit for Roche.
!!
!! ARGUMENTS
!!
!!  
!!
!!***

subroutine Gravity_init()

  use Gravity_data
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none
  real :: totmass

#include "constants.h"

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshMe)
  ! Reals
  call PhysicalConstants_get("newton"    , grv_newton)
  call RuntimeParameters_get("mass1"     , grv_m1)
  call RuntimeParameters_get("mass2"     , grv_m2)
  call RuntimeParameters_get("omega"     , grv_ome)
  call RuntimeParameters_get("separation", grv_sepa)
  ! Flags
  call RuntimeParameters_get("useGravity"    , useGravity)
  call RuntimeParameters_get("useCentrifugal", useCentrifugal)
  call RuntimeParameters_get("useCoriolis"   , useCoriolis)

!==============================================================================

!==============================================================================

  totmass = grv_m1 + grv_m2

  grv_x1  = + grv_m2*grv_sepa/totmass
  grv_y1  = 0.
  grv_z1  = 0.

  grv_x2  = - grv_m1*grv_sepa/totmass
  grv_y2  = 0.
  grv_z2  = 0.

  return
end subroutine Gravity_init
