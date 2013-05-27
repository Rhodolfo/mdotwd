
subroutine Simulation_initBlock(blockId,myPE)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
    Multispecies_getSumFrac
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getBlkPtr,&
    Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"

  integer,INTENT(in) ::  blockId, myPE

!  Local variables

  real, external :: softening_func
  real :: abar, zbar                 ! something to do with sum of mass fractions
  real :: xx, yy, zz, acc_dist, don_dist
  logical, parameter :: useGuardCell = .TRUE.
  real,pointer :: solnData(:,:,:,:)

  real,allocatable,dimension(:) :: xCoordsCell,yCoordsCell,zCoordsCell
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ

  integer :: i, j, k, n
  integer, dimension(MDIM) :: iPosition   !for putting data with Grid_putData

  integer :: icount, jHi, jLo, put
  integer, parameter :: ifail = -1
  real, allocatable  :: rvec(:)                   ! for the random number generator
  integer            :: rvecSize=0             ! number of random numbers needed,
                                                     ! calculated below
  ! variables needed for the eos call
  real :: temp_zone, rho_zone, vel_zone, frac, pres_zone, bdry_zone
  real :: rho_don, vel_don, pres_don, mass_don, fact_don
  real :: rho_acc, vel_acc, pres_acc, mass_acc, fact_acc
  logical :: bdry_don,bdry_acc
  real :: ptot, eint, etot, gamma
  real, dimension(EOS_NUM)  :: eosData
  real :: sum_p, sum_rho, width, soften, dist, soften_xtra
  soften_xtra = 1.





  ! ----------------------------------------------------------------------------------------------

  ! Get the indices of the blocks
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoordsCell(sizeX))
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoordsCell(sizeY))
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoordsCell(sizeZ))

  call Grid_getBlkPtr(blockID,solnData,CENTER)
#if NSPECIES > 0
  solnData(SPECIES_BEGIN,:,:,:) =  1.0-(NSPECIES-1)*sim_smallX
  solnData(SPECIES_BEGIN+1:SPECIES_END,:,:,:) =     sim_smallX
#endif
  call Grid_getCellCoords(KAXIS, blockId, CENTER, useGuardCell, zCoordsCell, sizeZ)
  call Grid_getCellCoords(JAXIS, blockId, CENTER, useGuardCell, yCoordsCell, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, useGuardCell, xCoordsCell, sizeX)
  
  !..get a blocks worth of random numbers between 0.0 and 1.0
    
  icount = 0

  ! Z-axis
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     if (NDIM == 3) then
        iPosition(3) = k
        zz = zCoordsCell(k)
     else 
        zz = 0.0
     endif
    
     ! Y-axis
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        if (NDIM >= 2) then
           iPosition(2) = j
           yy = yCoordsCell(j)
        else
           yy = 0.0
        endif

        !X-axis
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           iPosition(1) = i
           xx = xCoordsCell(i)
           icount = icount + 1

         ! Compute the distance from the center of each component
           acc_dist = sqrt((xx-sim_acc_center)**2 + (yy)**2 + (zz)**2)
           don_dist = sqrt((xx-sim_don_center)**2 + (yy)**2 + (zz)**2)

         ! ============================
         ! = Accretor Density Profile = 
         ! ============================

             call sim_find(sim_acc_rProf,SIM_NPROFILE,acc_dist,jLo)
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
                frac = (acc_dist - sim_acc_rProf(jlo)) &
                     / (sim_acc_rProf(jhi)-sim_acc_rProf(jlo))
             endif
             pres_acc = sim_acc_pProf(jlo) + &
                        frac*(sim_acc_pProf(jhi)  - sim_acc_pProf(jlo))
             rho_acc  = sim_acc_rhoProf(jlo) + &
                        frac*(sim_acc_rhoProf(jhi)- sim_acc_rhoProf(jlo))
             vel_acc  = sim_acc_vProf(jlo) + &
                        frac*(sim_acc_vProf(jhi) - sim_acc_vProf(jlo))                         
             mass_acc = sim_acc_mProf(jlo) + &
                        frac*(sim_acc_mProf(jhi) - sim_acc_mProf(jlo)) 
             if (acc_dist.gt.sim_acc_radius) then 
               soften   = soften_xtra &
                        * softening_func(acc_dist/sim_acc_radius)
               rho_acc  = rho_acc*soften
             ! pres_acc = pres_acc*(soften**(1.+1./sim_acc_n))
             end if

         ! =========================
         ! = Donor Density Profile = 
         ! =========================

             call sim_find(sim_don_rProf,SIM_NPROFILE,don_dist,jLo)
             soften = 1.0
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
                frac = (don_dist - sim_don_rProf(jlo)) &
                     / (sim_don_rProf(jhi)-sim_don_rProf(jlo))
             endif
             pres_don = sim_don_pProf(jlo) + &
                         frac*(sim_don_pProf(jhi)  - sim_don_pProf(jlo))
             rho_don  = sim_don_rhoProf(jlo) + &
                        frac*(sim_don_rhoProf(jhi)- sim_don_rhoProf(jlo))
             vel_don  = sim_don_vProf(jlo) + &
                        frac*(sim_don_vProf(jhi) - sim_don_vProf(jlo)) 
             mass_don = sim_don_mProf(jlo) + &
                        frac*(sim_don_mProf(jhi) - sim_don_mProf(jlo)) 
             if (don_dist.gt.sim_don_radius) then 
               soften   = soften_xtra &
                        * softening_func(don_dist/sim_don_radius)
 
               rho_don  = rho_don*soften
             ! pres_don = pres_don*(soften**(1.+1./sim_don_n))
             end if



         ! ======================
         ! = Fill in the Fields =
         ! ======================

             if (fillDon) then
             else
               pres_don = 0.0
               rho_don  = 0.0
             end if
             if (fillAcc) then
             else
               pres_acc = 0.0
               rho_acc  = 0.0
             end if
             pres_zone = pres_don + pres_acc
             rho_zone  = rho_don  + rho_acc
             vel_zone = 0.

         ! =======================
         ! = Boundary conditions = 
         ! =======================

             bdry_zone = -1. !! -1 for fluid cells
           ! Decide between mass or radius criterion for boundary
             if (useMassShl) then 
               bdry_don  = (mass_don.lt.sim_don_bdry*sim_don_mass).and.useBdryDon
               bdry_acc  = (mass_acc.lt.sim_acc_bdry*sim_acc_mass).and.useBdryAcc  
             else
               bdry_don  = (don_dist.lt.sim_acc_bdry*sim_don_radius).and.useBdryDon
               bdry_acc  = (acc_dist.lt.sim_don_bdry*sim_acc_radius).and.useBdryAcc 
             end if
           ! Fill in values inside boundary
             if (bdry_don.or.bdry_acc) then 
               bdry_zone = +1. ! 1 for solid cells 
               if (acc_dist.le.sim_acc_radius) then 
                 rho_zone  = sim_acc_inrho
                 pres_zone = sim_acc_inpres
               else if (don_dist.le.sim_don_radius) then 
                 rho_zone  = sim_don_inrho
                 pres_zone = sim_don_inpres
               else 
                 stop "YOU BASTARD"
               end if
               vel_zone  = 0.
             end if
           ! if (rho_zone.le.sim_smallRho) then 
             ! soften    = (sim_smallRho/rho_zone)**(1.+1./sim_don_n)
               rho_zone  = 1.001*max(rho_zone,sim_smallRho)
             ! pres_zone = pres_zone*soften
           ! end if



           ! CALL EOS to fill in pres, eint, etc
           eosData(EOS_DENS) = rho_zone
           eosData(EOS_PRES) = pres_zone
           ! USE DENS, TEMP to set other variables...
           call Eos(MODE_DENS_PRES,1,eosData,sim_xn)
           temp_zone = eosData(EOS_TEMP)
           rho_zone = eosData(EOS_DENS)
           ptot = eosData(EOS_PRES)
           eint = eosData(EOS_EINT)
           gamma = eosData(EOS_GAMC)

           ! calculate kinetic energy and total energy
           etot = eint + 0.5*vel_zone**2

           ! store the values
           ! fill the flash arrays
           solnData(TEMP_VAR,i,j,k) = temp_zone
           solnData(DENS_VAR,i,j,k) = rho_zone
           solnData(PRES_VAR,i,j,k) = ptot
           solnData(EINT_VAR,i,j,k) = eint
           solnData(ENER_VAR,i,j,k) = etot
           solnData(GAMC_VAR,i,j,k) = gamma
           solnData(GAME_VAR,i,j,k) = (ptot/(etot+rho_zone)+1.)
           solnData(VELX_VAR,i,j,k) = 0.
           solnData(VELY_VAR,i,j,k) = 0.
           solnData(VELZ_VAR,i,j,k) = 0.
           solnData(BDRY_VAR,i,j,k) = bdry_zone
           do n = SPECIES_BEGIN,SPECIES_END
              solnData(n,i,j,k) = sim_xn(n)
           enddo
                      
           !..end of 3d loops
        enddo  ! end of k loop
     enddo     ! end of j loop
  enddo        ! end of i loop

  ! cleanup
  deallocate(xCoordsCell)
  deallocate(yCoordsCell)
  deallocate(zCoordsCell)
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)


  return
end subroutine Simulation_initBlock


!******************************************************************************

!  Routine:     sim_find()

!  Description: Given a monotonically increasing table x(nn) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine sim_find (x, nn, x0, i)

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
end subroutine sim_find





function softening_func(x) 
  implicit none
  real :: softening_func,x
  softening_func = 1e-3*(x**(-2))
end function softening_func
