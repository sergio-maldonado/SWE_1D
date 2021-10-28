!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Equations to be solved, i.e. particular form of the SWEs			   	  !
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: February 2018                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Subroutines related to the actual equations to be solved (i.e. particular form of SWE - underground, unbalanced, etc.)
!NOTES:

!Tested OK (09/01/18)

module equations
	implicit none
		
contains
	
	!========= Computes fluxes as function of conserved and primitive variables ===================
	subroutine fluxes(h,vel,is_wet,flux)
		!Standard formulation conservative form 
		use parameters
		implicit none
		real*8, dimension(:), intent(in) :: h,vel
		real*8, dimension(:,:), intent(inout) :: flux
		integer, dimension(:), intent(in) :: is_wet
		
		flux(1,:) = is_wet(:)*vel(:)*h(:)
		flux(2,:) = is_wet(:)*( vel(:)*vel(:)*h(:) + 0.5d0*g*h(:)*h(:) )
		
	end subroutine fluxes
	
	
	!========= Computes Valiani's slope-related fluxes ===================
	subroutine fluxes_valiani(h,is_wet,flux)
		use parameters
		implicit none
		real*8, dimension(:), intent(in) :: h
		real*8, dimension(:,:), intent(inout) :: flux
		integer, dimension(:), intent(in) :: is_wet
		
		flux(1,:) = 0.0d0
		flux(2,:) = is_wet(:)*(0.5d0*g*h(:)*h(:) )
		
	end subroutine fluxes_valiani
	
end module equations