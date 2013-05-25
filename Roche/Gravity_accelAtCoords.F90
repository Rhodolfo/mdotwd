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
                          grv_ome,grv_com,grv_newton,&
                          useGravity,useCentrifugal,useCoriolis,&
                          useMass1,useMass2

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
  real :: dr,dx
  real :: xx1,yy1,zz1,xx2,yy2,zz2,cent,cori,g_1,g_2
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
  xx1 = grv_x1
  yy1 = grv_y1
  zz1 = grv_z1
  xx2 = grv_x2
  yy2 = grv_y2
  zz2 = grv_z2
! Assume we are on plane for 2D calculations, Roche potential useless in 1D
  yy  = 0.
  zz  = 0.
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
       if (useMass1) then 
         dr  = sqrt((xx(ii)-xx1)**2 + (yy(ii)-yy1)**2 + (zz(ii)-zz1)**2)
         dx  = xx(ii)-xx1
         g_1 = - grv_newton*grv_m1*dx/(dr**3)
         else
         g_1 = 0.
       end if
     ! Mass 2
       if (useMass2) then 
         dr  = sqrt((xx(ii)-xx2)**2 + (yy(ii)-yy2)**2 + (zz(ii)-zz2)**2)
         dx  = xx(ii)-xx2
         g_2 = - grv_newton*grv_m2*dx/(dr**3)
         else
         g_2 = 0.
       end if
     ! Centrifugal
       if (useCentrifugal) then 
         cent = (grv_ome**2)*(xx(ii)-grv_com)
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
       if (useMass1) then 
         dr  = sqrt((xx(ii)-xx1)**2 + (yy(ii)-yy1)**2 + (zz(ii)-zz1)**2)
         dx  = yy(ii)-yy1
         g_1 = - grv_newton*grv_m1*dx/(dr**3)
         else
         g_1 = 0.
       end if
     ! Mass 2
       if (useMass2) then 
         dr  = sqrt((xx(ii)-xx2)**2 + (yy(ii)-yy2)**2 + (zz(ii)-zz2)**2)
         dx  = yy(ii)-yy2
         g_2 = - grv_newton*grv_m2*dx/(dr**3)
         else
         g_2 = 0.
       end if
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
       if (useMass1) then 
         dr  = sqrt((xx(ii)-xx1)**2 + (yy(ii)-yy1)**2 + (zz(ii)-zz1)**2)
         dx  = zz(ii)-zz1
         g_1 = - grv_newton*grv_m1*dx/(dr**3)
         else
         g_1 = 0.
       end if
     ! Mass 2
       if (useMass2) then 
         dr  = sqrt((xx(ii)-xx2)**2 + (yy(ii)-yy2)**2 + (zz(ii)-zz2)**2)
         dx  = zz(ii)-zz2
         g_2 = - grv_newton*grv_m2*dx/(dr**3)
         else
         g_2 = 0.
       end if 
     ! Centrifugal
       cent = 0.
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
