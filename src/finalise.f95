!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocates variables, closes files	  								  !
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: February 2018                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module finalise
	implicit none
	
	
contains
	
	!========= Deallocates memory from all variables ====================
	subroutine dealloca_vars()
		use variables
		use parameters
		implicit none
		
		deallocate(q,q_e,q_w,q_e_bar,q_w_bar,q_as,q_old)
		deallocate(f,fv_e,fv_w)
		deallocate(x,zb,zb_e,zb_w,eta,eta_bar)
		deallocate(hr,hl,ur,ul,hstar,ustar,sl,sr,sstar,is_wet)
		deallocate(hv_e,hv_w)
		
		if (is_sensors) then
			deallocate(sens_loc)
			close(70)
		end if
		
		if (bc_w == 10) deallocate(ts_w)
		if (bc_e == 10) deallocate(ts_e)
		
		close(50)
		close(51)
		close(52)
		
		close(60)
		
	end subroutine dealloca_vars
	
	
	
	!========= Footnotes at end of simulation ====================
	subroutine footnote(path,time_start)
		implicit none
		real*8, intent(in) :: time_start
		character (len=200), intent(in) :: path
		real*8 :: time_finish
		
		call cpu_time(time_finish)
		
		write(*,*)' '
		write(*,*)'Simulation finished.'
		write(*,'(A,1X,EN12.3)')'Total CPU time of simulation (s): ',time_finish-time_start
		write(*,*)'See simulation.log for simulation details.'
		write(*,*)' '
		write(*,*)'Working directory is:'
		write(*,*)trim(path)
		write(*,*)' '
		write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
		
		!log file:
		write(60,*)' '
		write(60,*)'Simulation finished.'
		write(60,'(A,1X,EN12.3)')'Total CPU time of simulation (s): ',time_finish-time_start
		
	end subroutine footnote

	
	!========= Calls all previous subroutines ====================
	subroutine end_simulation()
		use variables
		use parameters
		implicit none

		call footnote(wkpath,time_start)
		call dealloca_vars()

	end subroutine end_simulation
		
		
	
end module finalise