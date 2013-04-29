!!****if* source/physics/Gravity/GravityMain/Roche/Gravity_accelAtCoords
!!
!! NAME
!!
!!  Gravity_accelAtCoords 
!!
!! SYNOPSIS
!!
!!  Gravity_accelAtCoords(integer(IN) :: numPoints,
!!                      real(IN)      :: iCoords(:),
!!                      real(IN)      :: jCoords(:),
!!                      real(IN)      :: kCoords(:),
!!                      integer(IN)   :: accelDir,
!!                      real(OUT)     :: accel(numPoints),
!!                      integer(IN)   :: blockID,
!!                      integer(IN),optional :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration in a
!!  specified direction for a vector of points given by their
!!  coordinates.
!!
!! ARGUMENTS
!!
!!  iCoords,jCoords,kCoords: coordinates of the points where the
!!                           gravitational accelation is requested.
!!                           Each of these arrays should either be
!!                           of lenght numPoints (or more), in which
!!                           case its nth value is used for the nth
!!                           point; or else of dimension 1, in which
!!                           case the value is used for all points.
!!  accelDir :    The acceleration direction:  allowed values are 
!!              IAXIS, JAXIS and IAXIS. These values are defined
!!              in constants.h.
!!  numPoints :  Number of cells to update in accel()
!!  accel     :   Array to receive results
!!  blockID  :  The local identifier of the block to work on,
!!                not applicable in pointmass gravity.
!!  potentialIndex :  optional, not applicable in pointmass gravity
!! 
!!***

subroutine Gravity_accelAtCoords (numPoints, iCoords,jCoords,kCoords, accelDir,&
     accel, blockID, &
     potentialIndex)

!=======================================================================

  use Gravity_data, ONLY: grv_m1,grv_x1,grv_y1,grv_z1,&
                          grv_m2,grv_x2,grv_y2,grv_z2,&
                          grv_ome,grv_newton,&
                          useGravity,useCentrifugal,useCoriolis

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: accelDir, numPoints
  real, dimension(:),INTENT(in) :: iCoords,jCoords,kCoords
  real, dimension(numPoints),INTENT(OUT) :: accel
  integer, intent(IN),optional :: blockID
  integer, intent(IN),optional :: potentialIndex

!==========================================================================

#ifdef FIXEDBLOCKSIZE
  real,dimension(numPoints) ::xx,yy,zz
#else
  real,allocatable,dimension(:) ::xx,yy,zz
#endif
  real :: dr32_1,tmpdr32_1,dr32_2,tmpdr32_2
  integer :: ii

!==============================================================================
  if (.NOT.useGravity) then
     accel(1:numPoints) = 0.0
     return
  end if

#ifndef FIXEDBLOCKSIZE
  allocate(xx(numPoints))
  allocate(yy(numPoints))
  allocate(zz(numPoints))
#endif



! Extracting coordinates
  zz = 0.
  yy = 0.
  if (NDIM == 3) then 
     if (size(kCoords) .GE. numPoints) then
        zz(1:numPoints) = kCoords(1:numPoints)
     else
        zz(1:numPoints) = kCoords(1)
     end if
  endif
  if (NDIM >= 2) then
     if (size(jCoords) .GE. numPoints) then
        yy(1:numPoints) = jCoords(1:numPoints)
     else
        yy(1:numPoints) = jCoords(1)
     end if
  endif
  if (size(iCoords) .GE. numPoints) then
     xx = iCoords(1:numPoints)
  else
     xx = iCoords(1)
  end if





  if (accelDir .eq. IAXIS) then                       ! x-component
     do ii = 1, numPoints
     ! Mass 1
       tmpdr32_1 = sqrt((xx(ii)-xx1)**2 + (yy(ii)-yy1)**2 + (zz(ii)-zz1)**2)
       dr32_1    = tmpdr32_1**3
       g_1       = - grv_newton*grv_m1*(xx(ii)-xx1)/dr32_1
     ! Mass 2
       tmpdr32_2 = sqrt((xx(ii)-xx2)**2 + (yy(ii)-yy2)**2 + (zz(ii)-zz2)**2)
       dr32_2    = tmpdr32_2**3
       g_2       = - grv_newton*grv_m2*(xx(ii)-xx2)/dr32_2
     ! Centrifugal
       if (useCentrifugal) then 
         cent = (grv_ome**2)*xx(ii)
         else
         cent = 0.
       end if
     ! Coriolis
       if (useCoriolis) then 
         stop "Coriolis not supported in Gravity_accelAtCoords"
         else
         cori = 0. 
       end if
     ! All together
       accel(ii) = g_1 + g_2 + cent + cori
     end do





  else if (accelDir .eq. JAXIS) then          ! y-component

     do ii = 1, numPoints
     ! Mass 1
       tmpdr32_1 = sqrt((xx(ii)-xx1)**2 + (yy(ii)-yy1)**2 + (zz(ii)-zz1)**2)
       dr32_1    = tmpdr32_1**3
       g_1       = - grv_newton*grv_m1*(yy(ii)-yy1)/dr32_1
     ! Mass 2
       tmpdr32_2 = sqrt((xx(ii)-xx2)**2 + (yy(ii)-yy2)**2 + (zz(ii)-zz2)**2)
       dr32_2    = tmpdr32_2**3
       g_2       = - grv_newton*grv_m2*(yy(ii)-yy2)/dr32_2
     ! Centrifugal
       if (useCentrifugal) then 
         cent = (grv_ome**2)*yy(ii)
         else
         cent = 0.
       end if
     ! Coriolis
       if (useCoriolis) then 
         stop "Coriolis not supported in Gravity_accelAtCoords"
         else
         cori = 0. 
       end if
     ! All together
       accel(ii) = g_1 + g_2 + cent + cori
     end do





  else if (accelDir .eq. KAXIS) then          ! z-component

     do ii = 1, numPoints
     ! Mass 1
       tmpdr32_1 = sqrt((xx(ii)-xx1)**2 + (yy(ii)-yy1)**2 + (zz(ii)-zz1)**2)
       dr32_1    = tmpdr32_1**3
       g_1       = - grv_newton*grv_m1*(zz(ii)-zz1)/dr32_1
     ! Mass 2
       tmpdr32_2 = sqrt((xx(ii)-xx2)**2 + (yy(ii)-yy2)**2 + (zz(ii)-zz2)**2)
       dr32_2    = tmpdr32_2**3
       g_2       = - grv_newton*grv_m2*(zz(ii)-zz2)/dr32_2
     ! Centrifugal
       if (useCentrifugal) then 
         cent = 0.
         else
         cent = 0.
       end if
     ! Coriolis
       if (useCoriolis) then 
         stop "Coriolis not supported in Gravity_accelAtCoords"
         else
         cori = 0. 
       end if
     ! All together
       accel(ii) = g_1 + g_2 + cent + cori 
     end do

  end if

!==============================================================================
#ifndef FIXEDBLOCKSIZE
  deallocate(xx)
  deallocate(yy)
  deallocate(zz)
#endif

  return

end subroutine Gravity_accelAtCoords
