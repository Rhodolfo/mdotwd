!!****if* source/physics/Gravity/GravityMain/Roche/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow 
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneRow(integer(2)(IN) :: pos,
!!                      integer(IN)    :: sweepDir,
!!                      integer(IN)    :: blockID,
!!                      integer(IN)    :: numCells,
!!                      real(:)(OUT)   :: grav,
!!                      integer(IN),optional :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration for a row
!!  of cells in a specified direction in a given block.
!!
!! ARGUMENTS
!!
!!  pos      :  Row indices transverse to the sweep direction
!!  sweepDir :    The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y and SWEEP_X. These values are defined
!!              in constants.h
!!  blockID  :  The local identifier of the block to work on
!!  numCells :  Number of cells to update in grav()
!!  grav()   :   Array to receive result
!!  potentialIndex :  optional, not applicable in pointmass gravity
!! 
!!***

subroutine Gravity_accelOneRow (pos,sweepDir,blockID, numCells, grav, &
      potentialIndex)

!=======================================================================

  use Gravity_data, ONLY: grv_m1,grv_x1,grv_y1,grv_z1,&
                          grv_m2,grv_x2,grv_y2,grv_z2,&
                          grv_ome,grv_com,grv_newton,&
                          useGravity,useCentrifugal,useCoriolis,&
                          useMass1,useMass2
  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkPtr, &
                             Grid_releaseBlkPtr,Grid_getBlkIndexLimits,Grid_getCellCoords

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: sweepDir,blockID,numCells
  integer, dimension(2),INTENT(in) ::pos
  real, dimension(numCells),INTENT(inout) :: grav
  integer,intent(IN),optional :: potentialIndex

!==========================================================================

#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_KHI_GC) :: zz
  real,dimension(GRID_JHI_GC) :: yy
  real,dimension(GRID_IHI_GC) :: xx
  real,POINTER,DIMENSION(:,:,:,:) :: solnVec  
#else
  real,allocatable,dimension(:)   :: xx,yy,zz
  real,dimension(numCells)        :: vx,vy
  real,dimension(MDIM)            :: blockSize
  real,POINTER,DIMENSION(:,:,:,:) :: solnVec 
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
#endif
  real, dimension(numcells)       :: vx,vy
  real    :: xx1,yy1,zz1,xx2,yy2,zz2
  real    :: dr32_1,tmpdr32_1,dr32_2,tmpdr32_2
  real    :: g_1,g_2,cent,cori
  integer :: sizeX,sizeY,sizez
  integer :: ii,j,k
  logical :: gcell = .true.

!==============================================================================

  if (.NOT.useGravity) return

  xx1 = grv_x1
  yy1 = grv_y1
  zz1 = grv_z1
  xx2 = grv_x2
  yy2 = grv_y2
  zz2 = grv_z2

! Getting the values for the fields in FLASH
  j=pos(1)
  k=pos(2)
#ifndef FIXEDBLOCKSIZE
  call Grid_getBlkPhysicalSize(blockID,blockSize)  
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)
  sizeY=blkLimitsGC(HIGH,JAXIS)
  sizeZ=blkLimitsGC(HIGH,KAXIS)
  allocate(xx(sizeX))
  allocate(yy(sizeY))
  allocate(zz(sizeZ))
#else
  call Grid_getBlkPtr(blockID,solnVec) 
  sizeX=GRID_IHI_GC
  sizeY=GRID_JHI_GC
  sizeZ=GRID_KHI_GC
#endif
  call Grid_getBlkPtr(blockID,solnVec)
  zz = 0.
  yy = 0.
  xx = 0.
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zz, sizeZ)
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yy, sizeY)
  endif
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xx, sizeX)
  




  if (sweepDir .eq. SWEEP_X) then                       ! x-component

     tmpdr32_1 = (yy(j)-yy1)**2  + (zz(k)-zz1)**2
     tmpdr32_2 = (yy(j)-yy2)**2  + (zz(k)-zz2)**2
     vx(:) = solnVec(VELX_VAR,:,pos(1),pos(2))
     vy(:) = solnVec(VELY_VAR,:,pos(1),pos(2)) 
     do ii = 1, numCells
      ! Mass 1
        if (useMass1) then 
          dr32_1 = sqrt((xx(ii)-xx1)**2 + tmpdr32_1)
          dr32_1 = (dr32_1)**3
          g_1    = - grv_newton*grv_m1*(xx(ii)-xx1)/dr32_1
          else
          g_1    = 0. 
        end if
      ! Mass 2
        if (useMass2) then
          dr32_2 = sqrt((xx(ii)-xx2)**2 + tmpdr32_2)
          dr32_2 = (dr32_2)**3
          g_2    = - grv_newton*grv_m2*(xx(ii)-xx2)/dr32_2
          else
          g_2    = 0.
        end if
      ! Coriolis
        if (useCentrifugal) then 
          cent = (grv_ome**2)*(xx(ii)-grv_com)
          else
          cent = 0.
        end if
      ! Coriolis
        if (useCoriolis) then 
          cori = + 2.*grv_ome*vy(ii)
          else
          cori = 0.
        end if
      ! All together
        grav(ii) = g_1 + g_2 + cent + cori
     enddo





  else if (sweepDir .eq. SWEEP_Y) then          ! y-component

     tmpdr32_1 = (xx(j)-xx1)**2  + (zz(k)-zz1)**2
     tmpdr32_2 = (xx(j)-xx2)**2  + (zz(k)-zz2)**2
     vx(:) = solnVec(VELX_VAR,pos(1),:,pos(2))
     vy(:) = solnVec(VELY_VAR,pos(1),:,pos(2)) 
     do ii = 1, numCells
      ! Mass 2
        if (useMass1) then
          dr32_1 = sqrt((yy(ii)-yy1)**2 + tmpdr32_1)
          dr32_1 = (dr32_1)**3
          g_1    = - grv_newton*grv_m1*(yy(ii)-yy1)/dr32_1
          else
          g_1    = 0.
        end if
      ! Mass 2
        if(useMass2) then 
          dr32_2 = sqrt((yy(ii)-yy2)**2 + tmpdr32_2)
          dr32_2 = (dr32_2)**3
          g_2    = - grv_newton*grv_m2*(yy(ii)-yy2)/dr32_2
          else
          g_2    = 0.
        end if
      ! Centrifugal
        if (useCentrifugal) then 
          cent = (grv_ome**2)*yy(ii)
          else
          cent = 0.
        end if
      ! Coriolis
        if (useCoriolis) then 
          cori = - 2.*grv_ome*vx(ii)
          else
          cori = 0.
        end if
      ! All together
        grav(ii) = g_1 + g_2 + cent + cori
     enddo





  else if (sweepDir .eq. SWEEP_Z) then          ! z-component

     tmpdr32_1 = (xx(j)-xx1)**2  + (yy(k)-yy1)**2
     tmpdr32_2 = (xx(j)-xx2)**2  + (yy(k)-yy2)**2
     vx(:) = solnVec(VELX_VAR,pos(1),pos(2),:)
     vy(:) = solnVec(VELY_VAR,pos(1),pos(2),:)  
     do ii = 1, numCells
      ! Mass 1
        if (useMass1) then 
          dr32_1 = sqrt((zz(ii)-zz1)**2 + tmpdr32_1)
          dr32_1 = (dr32_1)**3
          g_1    = - grv_newton*grv_m1*(zz(ii)-zz1)/dr32_1
          else 
          g_1    = 0.
        end if
      ! Mass 2
        if (useMass2) then 
          dr32_2 = sqrt((zz(ii)-zz2)**2 + tmpdr32_2)
          dr32_2 = (dr32_2)**3
          g_2    = - grv_newton*grv_m2*(zz(ii)-zz2)/dr32_2
          else
          g_2    = 0.
        end if
      ! Centrifugal and Coriolis are zero on z axis
        cent = 0.
        cori = 0.
      ! Add them up
        grav(ii) = g_1 + g_2
     enddo

  endif





!==============================================================================

! Relsease the solnVector pointer, this is important
  call Grid_releaseBlkPtr(blockID, solnVec)


#ifndef FIXEDBLOCKSIZE
  deallocate(xx)
  deallocate(yy)
  deallocate(zz)
#endif

  return

end subroutine Gravity_accelOneRow
