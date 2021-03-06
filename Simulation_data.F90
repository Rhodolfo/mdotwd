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

  ! Flags
  logical, save :: useBdryDon
  logical, save :: useBdryAcc
  logical, save :: useNStrAcc 
  logical, save :: useMassShl
  logical, save :: RedoDon
  logical, save :: fillDon
  logical, save :: fillAcc
  ! Binary system parameters here
  real, save :: sim_omega
  real, save :: sim_separ
  real, save :: sim_eggle
  real, save :: sim_ratio
  real, save :: sim_centmass
  real, save :: sim_L1
  ! Accretor parameters go here
  real, save :: sim_acc_mass
  real, save :: sim_acc_radius
  real, save :: sim_acc_bdry
  real, save :: sim_acc_center
  real, save :: sim_acc_n
  real, save :: sim_acc_c
  real, save :: sim_acc_rhoc
  real, save :: sim_acc_inrho
  real, save :: sim_acc_inpres
  ! Donor parameters here
  real, save :: sim_don_mass
  real, save :: sim_don_radius
  real, save :: sim_don_center
  real, save :: sim_don_bdry
  real, save :: sim_don_n
  real, save :: sim_don_c
  real, save :: sim_don_rhoc
  real, save :: sim_don_inrho
  real, save :: sim_don_inpres
  ! Minimum stuff
  real, save :: sim_smallrho
  real, save :: sim_smallx
  ! Domain boundaries
  real, save :: sim_xmax
  real, save :: sim_xmin
  ! Relaxation vars
  real   , save :: sim_trelax
  real   , save :: sim_relaxrate
  logical, save :: sim_relax
  !! *** For the Lane-Emdem solver *** !!
  integer, parameter                  :: SIM_NPROFILE=500
  real, dimension(NSPECIES)    , save :: sim_xn
  real, dimension(SIM_NPROFILE), save :: sim_acc_rProf   ,sim_don_rProf   ! Radial
  real, dimension(SIM_NPROFILE), save :: sim_acc_rhoProf ,sim_don_rhoProf ! Density
  real, dimension(SIM_NPROFILE), save :: sim_acc_pProf   ,sim_don_pProf   ! Pressure
  real, dimension(SIM_NPROFILE), save :: sim_acc_vProf   ,sim_don_vProf   ! Velocity
  real, dimension(SIM_NPROFILE), save :: sim_acc_cProf   ,sim_don_cProf   ! Sound speed
  real, dimension(SIM_NPROFILE), save :: sim_acc_mProf   ,sim_don_mProf   ! Mass 
  real, dimension(SIM_NPROFILE), save :: sim_acc_bProf   ,sim_don_bProf   ! Boundary criterion
  integer, parameter        :: np = 100000
  logical, save :: sim_gCell
  integer, save :: sim_meshMe

end module Simulation_data
