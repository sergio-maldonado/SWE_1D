!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters used in the main program									  !
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: February 2018                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module parameters
	
implicit none

integer, parameter :: m = 2
real*8, parameter :: g = 9.81d0
real*8, parameter :: zero = 1.0d-9
real*8, parameter :: ws = 10.0d0		

end module parameters

!!!!!!!!!!!!!!!!!!!!! GLOSARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! m			=		number of main variable dimensions (conservation laws), i.e. 2 in 1D
! g			=		gravitational acceleration
! zero		=		considered as 0 within computations
! ws		=		interval in % to print results to screen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
