!TO DO...
!... write documentation, perhaps using Matlab Live script
!... document file/paper with validation cases (use SWASHES)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shallow Water Equations Solver, 1D, Finite Volume, HLLC solver		  !
!                                                                         !
! -Includes bedslope term through method of Valiani (2006)				  !
! -MUSCL-Hancock time integration										  !
! -Friction solved through semi-implicit method				 			  !
!                                                                         !
! Author: Sergio Maldonado  							                  !
! s.maldonado@soton.ac.uk                                                 !
! University of Southampton, U.K.                                         !
!                                                                         !
! Created on: January 2018                                                !
!                                                                         !
!                                                                         !
! This software has been created strictly for academic (research and      !
! teaching) purposes.                                                     !
! This code is free software; you can redistribute it and/or              !
! modify it under the terms of the GNU Lesser General Public              !
! License as published by the Free Software Foundation; either            !
! version 2.1 of the License, or (at your option) any later version.      !
!                                                                         !
! This software is distributed in the hope that it will be useful,        !
! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
! Lesser General Public License for more details.                         !
!                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program main_program
	use variables
	use parameters
	
	use initialise
	use finalise
	use solver
	use tests_bugs
	use print_results
	
	implicit none
	
	!Reads input files, initial conditions, etc.
	call prep_simulation	!in initialise
	
	
	do while (tsim < tt)
		call main_computations	!in solver
		call all_print			!in print_results
	end do
	
	!TESTS:
	!call test_select	!in tests_bugs
	
	!postprocessing:
	call end_simulation
	
end program main_program
