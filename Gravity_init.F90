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
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none
  real :: Msun,m1,m2

#include "constants.h"

  ! Everybody should know these
  call Driver_getMype    (MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshMe)
  ! Parameters
  call PhysicalConstants_get("newton", grv_newton)
  call RuntimeParameters_get("mass1" , m1)
  call RuntimeParameters_get("mass2" , m2)
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
  grv_m1  = Msun*m1
  grv_m2  = Msun*m2

! ================================
! = Extract from Simulation_data = 
! ================================
  
  grv_com = sim_centmass
  grv_ome = sim_omega
  grv_x1  = sim_acc_center
  grv_x2  = sim_don_center

  return
end subroutine Gravity_init
