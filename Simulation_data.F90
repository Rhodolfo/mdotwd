!!****if* source/Simulation/SimulationMain/RhoSphere/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the spherical polytrope
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_rhoCenter  Central density
!!  sim_preCenter  Central pressure
!!  sim_temCenter  Central temperature
!!  gamma      Polytrope exponent

module Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"

  !! *** Runtime Parameters *** !!

  ! Binary system parameters here
  real, save :: sim_omega
  real, save :: sim_separ
  real, save :: sim_eggle
  real, save :: sim_ratio
  real, save :: sim_center_of_mass
  real, save :: sim_L1
  ! Accretor parameters go here
  real, save :: sim_acc_mass
  real, save :: sim_acc_radius
  real, save :: sim_acc_center
  real, save :: sim_acc_n
  real, save :: sim_acc_c
  real, save :: sim_acc_rhoc
  ! Donor parameters here
  logical, save :: sim_use_donor
  logical, save :: sim_use_rotation
  real, save :: sim_don_mass
  real, save :: sim_don_radius
  real, save :: sim_don_center
  real, save :: sim_don_n
  real, save :: sim_don_c
  real, save :: sim_don_rhoc
  ! Minimum stuff
  real, save :: sim_smallrho
  real, save :: sim_smallx
  ! Domain boundaries
  real, save :: sim_xmax
  real, save :: sim_xmin
  real, save :: sim_trelax
  real, save :: sim_relaxrate

  !! *** For the Lane-Emdem solver *** !!
  integer, parameter                     :: SIM_NPROFILE=500
  real, dimension(NSPECIES), save        :: sim_xn
  real, dimension(SIM_NPROFILE), save    :: sim_acc_rProf   ,sim_don_rProf
  real, dimension(SIM_NPROFILE), save    :: sim_acc_rhoProf ,sim_don_rhoProf
  real, dimension(SIM_NPROFILE), save    :: sim_acc_pProf   ,sim_don_pProf 
  real, dimension(SIM_NPROFILE), save    :: sim_acc_vProf   ,sim_don_vProf
  real, dimension(SIM_NPROFILE), save    :: sim_acc_cProf   ,sim_don_cProf
  integer, parameter        :: np = 100000
  logical, save :: sim_gCell
  integer, save :: sim_meshMe

end module Simulation_data
