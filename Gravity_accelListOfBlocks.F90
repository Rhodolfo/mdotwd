!!****if* source/physics/Gravity/GravityMain/Roche/Gravity_accelListOfBlocks
!!
!! NAME
!!
!!  Gravity_accelListOfBlocks  
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelListOfBlocks(integer(IN) :: blockCount,
!!                         integer(IN)    :: blockList(blockCount),
!!                         integer(IN)    :: component,
!!                         integer(IN)    :: accelIndex,
!!                         integer(IN),optional :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!! ARGUMENTS
!!
!!   blockCount   - The number of blocks in the list
!!   blockList    - The list of blocks on which to calculate acceleration.
!!   component         - The component of the acceleration to compute.
!!                  Permitted values are IAXIS, JAXIS, KAXIS.  In
!!                     theory, ALLDIR will be implemented.  At the moment
!!                     this routine aborts in a messy way if called with ALLDIR.
!!   accelIndex - variable # to store the acceleration
!!   potentialIndex :   Variable # to take as potential if present, 
!!                     not applicable to pointmass
!!
!!***

subroutine Gravity_accelListOfBlocks (blockCount,blockList,component, &
     accelIndex, potentialIndex)

!==============================================================================

  use Gravity_data, ONLY: grv_m1,grv_x1,grv_y1,grv_z1,&
                          grv_m2,grv_x2,grv_y2,grv_z2,&
                          grv_ome,grv_com,grv_newton,& 
                          useGravity,useCentrifugal,useCoriolis
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_getCellCoords, Grid_releaseBlkPtr
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"

  integer,intent(IN)                      :: blockCount
  integer,dimension(MAXBLOCKS), intent(IN)     :: blockList
  integer, INTENT(in) ::  component
  integer, INTENT(in) ::  accelIndex
  integer,intent(IN),optional :: potentialIndex

!==============================================================================



  integer       :: i, j, k, lb
  integer       :: csize
  logical       :: gcell = .true.
#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_IHI_GC) :: xx
  real,dimension(GRID_JHI_GC) :: yy
  real,dimension(GRID_KHI_GC) :: zz
#else
  real,allocatable,dimension(:) ::xx,yy,zz
#endif
  real, pointer :: solnVec(:,:,:,:)
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real    :: xx1,yy1,zz1,xx2,yy2,zz2,vx,vy
  real    :: dr32_1,tmpa_1,tmpb_1
  real    :: dr32_2,tmpb_2,tmpb_2 
  real    :: g_1,g_2,cent,cori

!==============================================================================
  if (.NOT.useGravity) return

  xx1 = grv_x1
  yy1 = grv_y1
  zz1 = grv_z1
  xx2 = grv_x2
  yy2 = grv_y2
  zz2 = grv_z2

  do lb = 1, blockCount
     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     
#ifndef FIXEDBLOCKSIZE
     allocate(xx(blkLimitsGC(HIGH,IAXIS)))
     allocate(yy(blkLimitsGC(HIGH,JAXIS)))
     allocate(zz(blkLimitsGC(HIGH,KAXIS)))
#endif
     call Grid_getBlkPtr(blockList(lb), solnVec)
!------------------------------------------------------------------------------

     csize = blkLimitsGC(HIGH, IAXIS)
     call Grid_getCellCoords(IAXIS,blockList(lb),CENTER, gcell,xx,csize)
     xx = xx
     yy = 0.
     zz = 0.
     if (NDIM >= 2) then
        csize = blkLimitsGC(HIGH, JAXIS)
        call Grid_getCellCoords(JAXIS,blockList(lb),CENTER,gcell,yy,csize)
        yy = yy
     endif
     if (NDIM == 3) then 
        csize = blkLimitsGC(HIGH, KAXIS)
        call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, gcell, zz, csize)
        zz = zz
     endif

!------------------------------------------------------------------------------
  
     if (component == IAXIS) then                    ! x-axis

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           tmpa_1 = (zz(k)-zz1)**2
           tmpa_2 = (zz(k)-zz2)**2
           do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              tmpb_1 = (yy(j)-yy1)**2
              tmpb_2 = (yy(j)-yy2)**2
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              ! Mass 1 
                dr32_1 = sqrt((xx(i)-xx1)**2+tmpa_1+tmpb_1)
                dr32_1 = dr32_1**3
                g_1    = - grv_newton*grv_m1*(xx(i)-xx1)/dr32_1
              ! Mass 2
                dr32_2 = sqrt((xx(i)-xx2)**2+tmpa_2+tmpb_2)
                dr32_2 = dr32_2**3
                g_2    = - grv_newton*grv_m2*(xx(i)-xx2)/dr32_2
              ! Centrifugal 
                if (useCentrifugal) then 
                  cent = (grv_ome**2)*(xx(i)-grv_com)
                  else
                  cent = 0.
                end if
              ! Coriolis 
                if (useCoriolis) then 
                  vx   = solnVec(VELX_VAR,i,j,k)
                  vy   = solnVec(VELY_VAR,i,j,k)
                  cori = + 2.*grv_ome*vy
                  else
                  cori = 0.
                end if
              ! All together
                solnVec(accelIndex,i,j,k) = g_1 + g_2 + cent + cori
              enddo
           enddo
        enddo

!------------------------------------------------------------------------------
  
     Elseif (component == JAXIS) then                ! y-axis

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           tmpa_1 = (zz(k)-zz1)**2
           tmpa_2 = (zz(k)-zz2)**2
           do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              tmpb_1 = (yy(j)-yy1)**2
              tmpb_2 = (yy(j)-yy2)**2
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              ! Mass 1 
                dr32_1 = sqrt((xx(i)-xx1)**2+tmpa_1+tmpb_1)
                dr32_1 = dr32_1**3
                g_1    = - grv_newton*grv_m1*(yy(i)-yy1)/dr32_1
              ! Mass 2
                dr32_2 = sqrt((xx(i)-xx2)**2+tmpa_2+tmpb_2)
                dr32_2 = dr32_2**3
                g_2    = - grv_newton*grv_m2*(yy(i)-yy2)/dr32_2
              ! Centrifugal 
                if (useCentrifugal) then 
                  cent = (grv_ome**2)*yy(i)
                  else
                  cent = 0.
                end if
              ! Coriolis 
                if (useCoriolis) then 
                  vx   = solnVec(VELX_VAR,i,j,k)
                  vy   = solnVec(VELY_VAR,i,j,k)
                  cori = - 2.*grv_ome*vx
                  else
                  cori = 0.
                end if
              ! All together
                solnVec(accelIndex,i,j,k) = g_1 + g_2 + cent + cori
              enddo
           enddo
        enddo

!------------------------------------------------------------------------------
  
     elseif (component == KAXIS) then                ! z-axis

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           tmpa_1 = (zz(k)-zz1)**2
           tmpa_2 = (zz(k)-zz2)**2
           do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              tmpb_1 = (yy(j)-yy1)**2
              tmpb_2 = (yy(j)-yy2)**2
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              ! Mass 1 
                dr32_1 = sqrt((xx(i)-xx1)**2+tmpa_1+tmpb_1)
                dr32_1 = dr32_1**3
                g_1    = - grv_newton*grv_m1*(zz(i)-zz1)/dr32_1
              ! Mass 2
                dr32_2 = sqrt((xx(i)-xx2)**2+tmpa_2+tmpb_2)
                dr32_2 = dr32_2**3
                g_2    = - grv_newton*grv_m2*(zz(i)-zz2)/dr32_2
              ! Centrifugal 
                cent = 0.
              ! Coriolis 
                cori = 0.
              ! All together
                solnVec(accelIndex,i,j,k) = g_1 + g_2 + cent + cori
              enddo
           enddo
        enddo

!------------------------------------------------------------------------------
  
     else                                        ! ALLAXIS
        call Driver_abortFlash &
             ('[Gravity_accelListOfBlocks] ALLAXIS not supported!')
     endif

!------------------------------------------------------------------------------
     
#ifndef FIXEDBLOCKSIZE
    deallocate(xx)
    deallocate(yy)
    deallocate(zz)
#endif
    call Grid_releaseBlkPtr(blockList(lb),solnVec)
    
 enddo

!==============================================================================

 return

end subroutine Gravity_accelListOfBlocks
