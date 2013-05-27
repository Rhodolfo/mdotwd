!!****if* source/Driver/DriverMain/Driver_sourceTerms
!!
!! NAME
!!
!!  Driver_sourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_sourceTerms(integer(IN)::blockCount,
!!                     integer(IN)::blockList(blockCount),
!!                     real(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Driver for source terms. Instead of calling all these routines 
!!  from Driver_evolveFlash we call Driver_sourceTerms which then
!!  makes the calls to Cool, Burn, Heat and Stir.  If a unit is not
!!  included in the simulation, the routine will be a stub and return
!!  without doing anything.
!! 
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the stirring operator
!!  dt           : the current timestep
!!
!!***



subroutine Driver_sourceTerms(blockCount, blockList, dt, pass)

! Default Driver_sourceTerms declarations
  use Polytrope_interface, ONLY : Polytrope
  use Driver_data, ONLY: dr_simTime
  use Flame_interface, ONLY : Flame_step
  use Stir_interface, ONLY : Stir
  use Heat_interface, ONLY : Heat
  use Heatexchange_interface, ONLY : Heatexchange
  use Burn_interface, ONLY : Burn
  use Cool_interface, ONLY : Cool
  use Ionize_interface, ONLY : Ionize
  use EnergyDeposition_interface, ONLY : EnergyDeposition
  use Deleptonize_interface, ONLY : Deleptonize
! Non-default declarations
  use Simulation_data, ONLY: sim_tRelax,sim_relax,sim_relaxrate,sim_xn,SIM_NPROFILE,&
                             sim_don_radius,sim_don_center,RedoDon,sim_don_inrho,&
                             sim_don_rProf,sim_don_pProf,sim_don_rhoProf,sim_don_inpres
  use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Driver_interface, ONLY : Driver_getMype
  use Eos_interface, ONLY : Eos



  implicit none
! Non-default includes
#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "Multispecies.h"
! Default Driver_sourceTerms declarations
  real, intent(IN)    :: dt
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN):: blockList
  integer, OPTIONAL, intent(IN):: pass
! Non-default declarations
  integer  ::  i, j, k, l, put, lb, ierr, istat, myPE, n, jLo, jHi
  real     ::  xx, yy, zz, xcenter, ycenter, zcenter
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: xn
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  real, dimension(:,:,:,:),pointer :: solnData
  real, dimension(EOS_NUM) :: eosData
  integer,dimension(MDIM) :: axis
  logical :: gcell = .true.,dotherelax,dotherefill,inbdry,indonor
  real :: relax_rate,dist,frac
  real :: accelx,x,y,r
  real :: rho_zone,pres_zone,ptot,eint,etot,gamma,temp_zone







! RELAX AND REFILL
  if (sim_relax.or.RedoDon) then 
    call Driver_getMype(MESH_COMM,myPE)
    dotherelax = (dr_simTime.lt.sim_tRelax.and.sim_relax.eqv..true.)
    if (dotherelax) then 
      relax_rate = 0.98d0 + 0.02d0*(dr_simTime/sim_tRelax)
      if (myPE.eq. 0) print*,"Relaxing, factor is",relax_rate
    end if
  ! Begin the loop
    do lb = 1, blockCount
      blkLimits = 0
      blkLimitsGC = 0
      call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC,CENTER) 
      sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
      allocate(xCoord(sizeX),stat=istat)
      sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
      allocate(yCoord(sizeY),stat=istat)
      sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
      allocate(zCoord(sizeZ),stat=istat)
      call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
      call Grid_getCellCoords(JAXIS, blockList(lb), CENTER, gcell, yCoord, sizeY)
      call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
      call Grid_getBlkPtr(blockList(lb),solnData,CENTER)
      do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
          do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
          ! RELAXING HERE, reduce by constant factor, then reset energy
            if (dotherelax) then 
              solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k)*relax_rate
              solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)*relax_rate
              solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k)*relax_rate
              solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
              0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
            end if
          ! REFILLING HERE, with full EOS
          ! Checking if within boundary domain
            if (solnData(BDRY_VAR,i,j,k).gt.0.) then 
              inbdry = .true.
              else
              inbdry = .false.
            end if
          ! Checking if within donor radius
            dist = sqrt((xCoord(i)-sim_don_center)**2+yCoord(i)**2+zCoord(i)**2)
            if (dist.le.sim_don_radius) then 
              indonor = .true. 
              else
              indonor = .false.
            end if
          ! Only refill if within donor and outside boundary
            dotherefill = indonor.and.RedoDon.and.(dr_simTime.gt.1e-14)
            if (dotherefill) then 
              call driver_find(sim_don_rProf,SIM_NPROFILE,dist,jLo)
              if (jLo .eq. 0) then
                jLo = 1
                jHi = 1
                frac = 0.
              else if (jLo .eq. SIM_NPROFILE) then
                jLo = SIM_NPROFILE
                jHi = SIM_NPROFILE
                frac = 0.
              else
                jHi = jLo + 1
                frac = (dist - sim_don_rProf(jlo)) &
                     / (sim_don_rProf(jhi)-sim_don_rProf(jlo))
              end if
            ! pres_zone = sim_don_pProf(jlo) + &
            !             frac*(sim_don_pProf(jhi)  - sim_don_pProf(jlo))
              rho_zone  = sim_don_rhoProf(jlo) + &
                          frac*(sim_don_rhoProf(jhi)- sim_don_rhoProf(jlo))
              if (inbdry) then 
            !   pres_zone = sim_don_inpres
                rho_zone  = sim_don_inrho
              end if
            ! eosData(EOS_DENS) = rho_zone
            ! eosData(EOS_PRES) = pres_zone
            ! USE DENS, TEMP to set other variables...
            ! call Eos(MODE_DENS_PRES,1,eosData,sim_xn)
            ! temp_zone = eosData(EOS_TEMP)
            ! rho_zone = eosData(EOS_DENS)
            ! ptot = eosData(EOS_PRES)
            ! eint = eosData(EOS_EINT)
            ! gamma = eosData(EOS_GAMC)
            ! calculate kinetic energy and total energy
            ! etot = eint
            ! fill the flash arrays
            ! solnData(TEMP_VAR,i,j,k) = temp_zone
            ! solnData(DENS_VAR,i,j,k) = rho_zone
            ! solnData(PRES_VAR,i,j,k) = ptot
            ! solnData(EINT_VAR,i,j,k) = eint
            ! solnData(ENER_VAR,i,j,k) = etot
            ! solnData(GAMC_VAR,i,j,k) = gamma
            ! solnData(GAME_VAR,i,j,k) = (ptot/(etot+rho_zone)+1.)
            ! solnData(VELX_VAR,i,j,k) = 0.
            ! solnData(VELY_VAR,i,j,k) = 0.
            ! solnData(VELZ_VAR,i,j,k) = 0.
            ! do n = SPECIES_BEGIN,SPECIES_END
            ! solnData(n,i,j,k) = sim_xn(n)
            ! enddo
            end if
          enddo
        enddo
      enddo
      call Grid_releaseBlkPtr(blockList(lb), solnData)
      deallocate(xCoord)
      deallocate(yCoord)
      deallocate(zCoord)
    enddo
  else
  end if




! Default Driver_sourceTerms calls
  call Polytrope(blockCount, blockList, dt)
  call Stir(blockCount, blockList, dt) 
  call Flame_step(blockCount, blockList, dt)
  call Burn(blockCount, blockList, dt) 
  call Heat(blockCount, blockList, dt, dr_simTime) 
  call Heatexchange(blockCount, blockList, dt)
  call Cool(blockCount, blockList, dt, dr_simTime)
  call Ionize(blockCount, blockList, dt, dr_simTime)
  call EnergyDeposition(blockCount, blockList, dt, dr_simTime, pass)
  call Deleptonize(blockCount, blockList, dt, dr_simTime)



  return
end subroutine Driver_sourceTerms





!******************************************************************************

!  Routine:     sim_find()

!  Description: Given a monotonically increasing table x(nn) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine driver_find (x, nn, x0, i)

  implicit none

! Arguments, LBR guessed intent on these
  integer, intent(IN) :: nn
  integer, intent(OUT):: i
  real, intent(IN)    :: x(nn), x0

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then

     i = 0

  elseif (x0 .gt. x(nn)) then

     i = nn

  else

     il = 1
     ir = nn
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   i = il

  endif

  return
end subroutine driver_find
