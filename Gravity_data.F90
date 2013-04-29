!!****if* source/physics/Gravity/GravityMain/Roche/Gravity_data
!!
!! NAME
!!
!!  Gravity_data
!!
!! SYNOPSIS
!!
!!  use Gravity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for the Roche potential.
!!
!! PARAMETERS
!!   grv_m1,grv_m2,grv_ome,grv_sepa,grv_com
!!
!!***

module Gravity_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: grv_m1,grv_x1,grv_y1,grv_z1
  real, save :: grv_m2,grv_x2,grv_y2,grv_z2
  real, save :: grv_ome,grv_sepa

  !! *** Physical Constants *** !!

  real, save :: grv_newton

  !! *** I totally don't know what this is ***

  integer, save :: grv_meshMe, grv_meshNumProcs

  !! *** Flags *** !!

  logical, save :: useGravity
  logical, save :: useCentrifugal
  logical, save :: useCoriolis

end module Gravity_data
