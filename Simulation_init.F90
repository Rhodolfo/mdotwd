!****if* source/Simulation/SimulationMain/Sod/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the spherical polytrope
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp

  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"
#include "Flash_mpi.h"



  integer :: myPE, i
  real    :: Msun, Rsun, Xhydrogen, Yhelium 
  real    :: G, mu

  ! -----------------------------------------------------------------------------------
  
#ifndef FIXEDBLOCKSIZE
  call Driver_abortFlash("SORRY, Cellular not defined for non fixed block size")
  real             rvec(blockExtent(HIGH,IAXIS)*blockExtent(HIGH,JAXIS)*blockExtent(HIGH,KAXIS))
#endif  
  G = 6.673e-8
  Msun = 1.9889225d33
  Rsun = 6.96d10
  !-----------------------------------------------------------------------------
  ! grab the parameters relevant for this problem
  !-----------------------------------------------------------------------------
  call Driver_getMype(MESH_COMM, sim_meshMe)
  myPE = sim_meshMe
  
  ! Binary parameters
  call RuntimeParameters_get('mass1',sim_acc_mass)
  call RuntimeParameters_get('npoly',sim_acc_n)
  call RuntimeParameters_get('mass2',sim_don_mass)
  call RuntimeParameters_get('npoly',sim_don_n)
  ! Fluff
  call RuntimeParameters_get('smlrho', sim_smallrho)
  call RuntimeParameters_get('smallX', sim_smallx)
  ! Domain boundaries
  call RuntimeParameters_get('xmax',sim_xmax)
  call RuntimeParameters_get('xmin',sim_xmin)
  ! Relaxation 
  call RuntimeParameters_get('sim_trelax'   ,sim_trelax)
  call RuntimeParameters_get('sim_relaxrate',sim_relaxrate)
  ! Compoistion
#if NSPECIES > 0  
  sim_xn(2:NSPECIES)     = sim_smallx
  sim_xn(1) = 1.e0 - (NSPECIES-1)*sim_smallx
#endif



  ! ========================================
  ! = Create the initial condition profile =
  ! ========================================

  mu = 2.0/(1.0 + 3.0*Xhydrogen + 0.5*Yhelium)
! Doing accretor
  call sim_lane_emden(sim_acc_n,sim_acc_mass,sim_acc_rhoc,mu,&
       sim_acc_rProf,sim_acc_radius,sim_acc_rhoProf,sim_acc_mProf,&
       sim_acc_pProf,sim_acc_vProf,sim_acc_cProf,sim_acc_c)
! Doing donor if included
  call sim_lane_emden(sim_don_n,sim_don_mass,sim_don_rhoc,mu,&
       sim_don_rProf,sim_don_radius,sim_don_rhoProf,sim_don_mProf,&
       sim_don_pProf,sim_don_vProf,sim_don_cProf,sim_don_c)

  ! =======================================
  !  (convert masses to cgs)
  ! =======================================
  sim_acc_mass = sim_acc_mass * msun
  sim_don_mass = sim_don_mass * msun



  ! print initial condition profile to file
  if(myPE .eq. 0) then
     open(45,file="initial_acc_prof.dat")
     open(46,file="initial_don_prof.dat")
     do i=1,SIM_NPROFILE
        write(45,*) i,sim_acc_rProf(i),sim_acc_rhoProf(i),sim_acc_pProf(i),sim_acc_vProf(i)
        write(46,*) i,sim_don_rProf(i),sim_don_rhoProf(i),sim_don_pProf(i),sim_don_vProf(i)
     enddo
     close(45)
     close(46)
  endif
  ! Calculating separation and stuff
  sim_ratio      = sim_don_mass / sim_acc_mass
  sim_eggle      = (0.49*(sim_ratio**(2./3.))) / &
                   (0.60*(sim_ratio**(2./3.))+log(1.+(sim_ratio**(1./3.))))  
  sim_separ      = sim_don_radius / sim_eggle
  sim_omega      = sqrt(G*(sim_acc_mass+sim_don_mass)/(sim_separ**3.)) 
  sim_acc_center = 0.
  sim_don_center = sim_acc_center-sim_separ
  sim_centmass   = (sim_acc_mass*sim_acc_center + sim_don_mass*sim_don_center) &
                 / ( sim_acc_mass + sim_don_mass)
  sim_L1         = 0.

! Writing data over here
  if (myPE .eq. 0) then
     open(unit=45,file="rho.par",status="unknown")
     i = 45
     write(i,*) "============================================================"
     write(i,*)  "Polytrope Parameters for Accretor"
     write(i,*)  "mu =",mu
     write(i,*)  "Gamma = 1+1/n = ", 1.0+1.0/sim_acc_n
     write(i,*)  "mass = ", sim_acc_mass, "  rhoc = ", sim_acc_rhoc
     write(i,*)  "Using mass, rhoc to compute polyk"
     write(i,*) "============================================================"
     write(i,*)  "Polytrope Parameters for Donor"
     write(i,*)  "mu =",mu
     write(i,*)  "Gamma = 1+1/n = ", 1.0+1.0/sim_don_n
     write(i,*)  "mass = ", sim_don_mass, "  rhoc = ", sim_don_rhoc
     write(i,*)  "Using mass, rhoc to compute polyk"
     write(i,*) "============================================================"
     write(i,*)  "Resulting Accretor Star Parameters"
     write(i,*)  "R_acc_      = ", sim_acc_radius
     write(i,*)  "sound_speed = ", sim_acc_c
     write(i,*)  "t_crossing  = ", sim_acc_radius / sim_acc_c
     write(i,*) "============================================================"
     write(i,*)  "Resulting Donor Star Parameters"
     write(i,*)  "R_acc_      = ", sim_don_radius
     write(i,*)  "sound_speed = ", sim_don_c
     write(i,*)  "t_crossing  = ", sim_don_radius / sim_don_c
     write(i,*) "============================================================" 
     write(i,*)  "Resulting Binary System Parameters"
     write(i,*)  "min_x,max_x = ", sim_xmin,sim_xmax
     write(i,*)  "Donor Loc   = ", sim_don_center
     write(i,*)  "Separation  = ", sim_separ
     write(i,*)  "Frequency   = ", sim_omega
     write(i,*)  "Period      = ", 6.28/sim_omega
     write(i,*)  "Mass Ratio  = ", sim_ratio
     write(i,*) "============================================================" 
     close(45)
   end if

! Do the data writing

  return
end subroutine Simulation_init



! adiabatic P(rho/rho0,P0), i.e. dS = 0
  real function Pres(densratio,P0)
  implicit none
  real,intent(in) :: densratio, P0
  Pres = P0 * densratio**(5.0/3.0) 
  end function Pres



  subroutine sim_lane_emden(ordern,polmass,rhoc,mu,rProf,polr,rhoProf,mProf,pProf,vProf,cProf,c)
  use Simulation_data, only: SIM_NPROFILE, sim_xmax
  implicit none
  real :: ordern,polmass,rhoc,mu,polr,c
  real :: Pres
  real :: xsurf, ypsurf, tol, polyk
  integer :: iend, ipos, mode, j
  integer :: n, nsteps, iprint 
  real, dimension(SIM_NPROFILE) :: xx, yy, yp, mass, ebind, &
       zopac, rhom, zbeta, ztemp, exact, radius, rhop, prss 
  ! R, u are FACE CENTERED range: 0,n
  real, dimension(:,:),allocatable :: R, u
  ! rho, P are CELL CENTERED range 1,n
  real, dimension(:,:), allocatable :: rho
  real, dimension(:), allocatable :: P, rho0, P0
  ! For that solver
  real, dimension(SIM_NPROFILE)    :: rProf  
  real, dimension(SIM_NPROFILE)    :: rhoProf
  real, dimension(SIM_NPROFILE)    :: pProf  
  real, dimension(SIM_NPROFILE)    :: vProf 
  real, dimension(SIM_NPROFILE)    :: cProf
  real, dimension(SIM_NPROFILE)    :: mProf
  character(len=MAX_STRING_LENGTH) :: profFile
 
  mode = 2
  tol = 1.e-10
  polyk = 1e13/(mu**(5./3.))
  ! Solve Lane-Emden
  call polytr(ordern,polmass,rhoc,polyk,mu,mode,&
       xx,yy,yp,radius,rhop,mass,prss,ebind,rhom, &
       ztemp,zbeta,exact,xsurf,ypsurf,SIM_NPROFILE,iend,ipos,tol)
  rhoc   = rhoc
  polr   = radius(ipos)
  ! NOW CREATE PROFILE FOR COLLISIONLESS COLLAPSE
  nsteps = 10000
  iprint = nsteps/100
  n = ipos
  ! ALLOCATE LOCAL rho,P,R,u arrays
  allocate(R(0:n,2))
  allocate(u(0:n,2))
  allocate(rho(1:n,2))
  allocate(P(1:n))
  allocate(P0(1:n))
  allocate(rho0(1:n)) 
  ! SET THE INITAL PROFILES
  P = 0
  u = 0
  R(0,1) = 0.0
  do j=1,n
     R(j,1) = radius(j)
  enddo
  rho(1,1) = rhop(1)
  P(1) = prss(1)
  do j=2,n
     rho(j,1) = 0.50*( rhop(j-1) + rhop(j) )
     P(j)     = 0.50*( prss(j-1) + prss(j) )
  enddo   
  rho0(:) = rho(:,1)
  P0(:)   = P(:)
  !interior of profile
  c = 0.0
  do j=1,n-1
     rProf(j)   = 0.50*(R(j,1)+R(j-1,1))
     rhoProf(j) = rho(j,1)
     mProf(j)   = mass(j)
     pProf(j)   = P(j)
     vProf(j)   = u(j,1)
     cProf(j)   = sqrt((1.0 + 1.0/ordern)*P(j)/rho(j,1))
     c          = c + cProf(j)
  end do
  c = c / n
  !exterior, continue variables... 
  if (sim_xmax .gt. rProf(n)) then
     do j=n,SIM_NPROFILE
        rProf(j) =  0.50*(R(n,1)+R(n-1,1)) &
             + (j-n)*(sim_xmax - 0.50*(R(n,1)+R(n-1,1)) )/(SIM_NPROFILE - n)
        rhoProf(j) = rho(n,1)
        mProf(j) = mass(n)
        pProf(j) = P(n)
        vProf(j) = u(n,1)
     end do
  else
     call Driver_abortFlash('Box smaller than star')
  endif
 ! CLEAN UP
  deallocate(rho)
  deallocate(P)
  deallocate(R)
  deallocate(u)
  return
  end subroutine sim_lane_emden
