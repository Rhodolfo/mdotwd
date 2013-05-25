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
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none
  real :: Msun 

#include "constants.h"

  ! Everybody should know these
  call Driver_getMype    (MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshMe)
  ! Parameters
  call PhysicalConstants_get("newton", grv_newton)
  call RuntimeParameters_get("mass1" , grv_m1)
  call RuntimeParameters_get("mass2" , grv_m2)
  call RuntimeParameters_get("x1"    , grv_x1)
  call RuntimeParameters_get("x2"    , grv_x2) 
  ! Flags
  call RuntimeParameters_get("useGravity"    , useGravity)
  call RuntimeParameters_get("useCentrifugal", useCentrifugal)
  call RuntimeParameters_get("useCoriolis"   , useCoriolis)
  call RuntimeParameters_get("useMass1"      , useMass1)
  call RuntimeParameters_get("useMass2"      , useMass2)
 
  if (.not.useGravity) return

! ================
! = Derived data =
! ================

  Msun    = 1.9889225d33
  grv_com = (grv_x1*grv_m1 + grv_x2*grv_m2)/(grv_m1 + grv_m2)
  grv_m1  = Msun*grv_m1
  grv_m2  = Msun*grv_m2
  grv_ome = sqrt(grv_newton*(grv_m1+grv_m2)/abs(grv_x1-grv_x2))

  grv_y1  = 0.
  grv_z1  = 0.
  grv_y2  = 0.
  grv_z2  = 0.

  return
end subroutine Gravity_init
