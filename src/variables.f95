!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Global variables used in the main program								  !
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: February 2018                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TO DO...
! ... arrange variables in 'semantic' groups
module variables
	
implicit none
! solver-related variables:
integer :: i,k,n
real*8, dimension(:,:), allocatable :: q,q_e,q_w,q_e_bar,q_w_bar,q_as,q_old
real*8, dimension(:), allocatable :: x,dx,zb,zb_e,zb_w,hr,hl,ur,ul,hstar,ustar,sl,sr,sstar,eta,hv_e,hv_w,eta_bar
real*8 :: tt,tsim,dt,manning
real*8, dimension(:,:), allocatable :: f,fv_e,fv_w
integer, dimension(:), allocatable :: is_wet
integer :: sl_variable

! numerical parameters:
real*8 :: courant,beta,hlim
integer :: bc_w,bc_e,which_dt

! time series, if present:
real*8, dimension(:,:), allocatable :: ts_w,ts_e
character(len=1) :: ts_var_w,ts_var_e

! for local recordings (sensors):
logical :: is_sensors
integer :: nsens
integer, dimension(:), allocatable :: sens_loc
real*8 :: tw_sens

! for writing to file and screen:
real*8 :: time_start
real*8 :: tw,percent_ws
integer :: steps,j_ws
character(len=200) :: wkpath
character (len=30) :: ic_file,format_global

end module variables

!!!!!!!!!!!!!!!!!!!!! GLOSARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! beta							needed in slope limiter (0 = no slope limiter)
! bc_e							type of boundary condition at East boundary
! bc_w							type of boundary condition at West boundary
! courant						Courant number for stability criterion
! dt							time interval for integration (variable)
! dx(:)							vector with cells' widths
! eta(:)						water surface level
! eta_bar(:)					predicted value of eta
! f(k,ii)						numeric (HLLC) flux through interface ii
! format_global					format for printing global results to file
! fv_e(k,i)						bedslope-related flux following Valiani at East face of cell i
! fv_w(k,i)						bedslope-related flux following Valiani at West face of cell i
! hl(ii)						water depth at left side of interface ii
! hr(ii)						water depth at right side of interface ii
! hlim							water depth considered as zero within the code
! hstar(ii)						water depth in the star region
! hv_e(i)						Valiani depth at East face of cell i
! hv_w(i)						Valiani depth at West face of cell i
! i 							index related to cell number
! ic_file						file with initial conditions
! is_sensors					are local sensors present?
! is_wet(ii)					determines whether interface ii is dry or wet
! j_ws							index for printing-to-screen purposes
! k								index related to nunmber of conserv. law
! manning						Manning's n
! n								number of cells
! nsens							number of local sensors recording
! percent_ws					pseudo percentage of progress for print-to-screen purposes
! q(:,:)						main variable to be integrated in time [u, uh]
! q_e(:,:)						variable at East face of cell
! q_w(:,:)						variable at West face of cell
! q_e_bar(:,:)					predicted value of q_e
! q_w_bar(:,:)					predicted value of q_w
! q_as(:,:)						integrated q before including bed friction
! q_old(:,:)					records old value of q (for use in boundary conditions)
! sl(ii)						left-ward wave speed
! sr(ii)						right-ward wave speed
! sstar(ii)						wave speed in the star region
! sl_variable					determines whether slope limiter is applied to eta or h
! sens_loc(i)					location of sensor i (index of cell)
! steps							stores simulation steps or number of iterations
! tt							total time of simulation
! tw							time  interval to print global results to files
! tsim							current simulation time
! ts_e(2,:)						stores input time series imposed at East boundary
! ts_w(2,:)						stores input time series imposed at West boundary
! tw_sens						recording/sampling time interval for sensors
! time_start					calls CPU time at start of simulation
! ts_var_e						variable (h or q) to be read as time series at East boundary
! ts_var_w						variable (h or q) to be read as time series at West boundary
! ul(ii)						water velocity at left side of interface ii
! ur(ii)						water velocity at right side of interface ii
! ustar(ii)						water velocity in the star region
! wkpath						working directory
! which_dt						criterion for estimating time step (e.g. all domain is wet?)
! x(:)							x coordinates
! zb(:)							bed elevation
! zb_e(:)						bed elevation at East face of cell
! zb_w(:)						bed elevation at West face of cell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



