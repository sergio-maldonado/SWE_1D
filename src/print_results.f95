!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prints results to files												  !
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: February 2018                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module print_results
	implicit none
		
contains
	
	!========= Prints global variables ===================
	subroutine print_global_vars(n,q,eta,forma_global)
		implicit none
		integer, intent(in) :: n
		real*8, dimension(2,-1:n+2), intent(in) :: q
		real*8, dimension(n), intent(in) :: eta
		character*30, intent(in) :: forma_global
		integer :: i
				
			write(50,forma_global)(q(1,i),i=1,n)
			write(51,forma_global)(q(2,i),i=1,n)
			write(52,forma_global)(eta(i),i=1,n)
	
	end subroutine
	
	
	
	!========= Prints local variables (sensors) ===================
	subroutine print_local_vars(n,q,nsens,loc,forma_global,tsim)
		implicit none
		integer, intent(in) :: n,nsens
		real*8, dimension(2,-1:n+2), intent(in) :: q
		character*30, intent(in) :: forma_global
		integer, dimension(:), intent(in) :: loc
		real*8, intent(in) :: tsim
		integer :: i
		
			write(70,forma_global)tsim,(q(1,loc(i)),i=1,nsens),(q(2,loc(i)),i=1,nsens)
		
	end subroutine
	
	
	
	!========= Prints progrress to screen ===================
	subroutine print_progress(tsim,totaltime,steps,percent_write,j)
		use parameters
		implicit none
		real*8, intent(in) :: tsim,totaltime
		real*8, intent(inout) :: percent_write
		integer, intent(in) :: steps
		integer, intent(inout) :: j
		real*8 :: percent,mean_dt

		percent = tsim*100.0d0/totaltime					!real percent
		percent_write = tsim*100.0d0/totaltime - j*ws		!fake percent for printing purposes
		mean_dt = tsim/steps
		
		if (percent_write > ws) then
			write(*,'(f5.1,a10,4x,a13,en12.3,a1)')percent,'% complete','average dt = ',mean_dt,'s'
			write(60,'(f5.1,a10,4x,a13,en12.3,a1)')percent,'% complete','average dt = ',mean_dt,'s' !log file
			percent_write = 0.0d0
			j = j + 1
		end if
		
	end subroutine
	
	
	
	!========= Calls previous subroutines ====================
	subroutine all_print()
		use variables
		use parameters
		implicit none
		
		if (mod(tsim,tw) < dt) call print_global_vars(n,q,eta,format_global)
		
		if (is_sensors) then
			if (mod(tsim,tw_sens) < dt) call print_local_vars(n,q,nsens,sens_loc,format_global,tsim)
		end if
		
		call print_progress(tsim,tt,steps,percent_ws,j_ws)
		
	end subroutine
	
end module print_results