!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tests to be carried out during a simulation							  !
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: February 2018                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module tests
	implicit none
	
	
contains
	
	!========= Tests for negative depths ====================
	subroutine test_neg_depth(h,text)
		use parameters
		implicit none
		real*8, dimension(:), intent(in) :: h
		character(len=*), intent(in) :: text
		integer :: i,n
		
		n = size(h)
		do i=1,n
			if(h(i) < 0.0d0) then
				write(*,*)'!!!!!!!! ERROR !!!!!!!!'
				write(*,*)'Negative water depth found!'
				write(*,'(A15,I4,2X,A9,EN12.3)')'... at index = ',i,'depth is:',h(i)
				write(*,*)text
				write(*,*)'Simulation terminated.'
				stop
			end if
		end do
		
	end subroutine test_neg_depth
	
	
	
	!========= Tests for high velocities ====================
	subroutine test_high_vel(bed,depth,disch,hlim)
		implicit none
		real*8, intent(in) :: hlim
		real*8, dimension(:), intent(in) :: bed,depth,disch
		integer :: i,n
		real*8 :: h
		logical :: switch
		
		switch = .false.
		n = size(bed)
		
		do i=1,n
			h = depth(i)
			if (h > hlim) then
				if(disch(i)/h > 10.0d0) then	!10m/s defined arbitrarily
					write(*,*)'Initial velocity higher than 10 m/s, is this correct?'
					write(*,*)'Check index',i
					write(*,'(A15,EN12.3,1X,EN12.3,A1)')'where [h,u] = [',h,disch(i)/h,']'
					write(*,'(A15,EN12.3,1X,EN12.3,A1)')'and [eta,q] = [',h+bed(i),disch(i),']'
					read(*,*)
					switch = .true.
				end if
			end if			
		end do
		
		
	end subroutine test_high_vel
	
	
	
	!========= Tests for negative or zero dx ====================
	subroutine test_dx(dx)
		use parameters
		implicit none
		real*8, dimension(:), intent(in) :: dx
		integer :: i,n
	
		n = size(dx)
		do i=1,n
			if (dabs(dx(i)) < zero) then
				write(*,*)'A negative or zero value of dx was found in bathymetry file!'
				stop
			end if
		end do
		
		
	end subroutine test_dx
	
	
	
	!========= Tests for division by zero ====================
	subroutine div_zero(number,zero,flag0)
		implicit none
		real*8, intent(in) :: number,zero
		logical, intent(out) :: flag0
		
		flag0 = .false.
		
		if (dabs(number) < zero) flag0 = .true.
		
	end subroutine div_zero
	
	
	
	!========= Panic subroutine ====================
	subroutine panic(code,text)
		implicit none
		integer, intent(in) :: code
		character(len=*), optional, intent(in) :: text
		
		if (code == 0) then
			write(*,*)'Something went wrong! Potentially division by zero. Check:'
			write(*,*)text
			write(*,*)' '
			stop
		else if (code == 1) then
			write(*,*)'Something went wrong with computation of numeric fluxes!'
			stop
		else if (code == 2) then
			write(*,*)'!!!!!!!! ERROR !!!!!!!!'
			write(*,*)'Unexpected error in:'
			write(*,*)text
			write(*,*)'Simulation aborted.'
			stop
		else if (code == 3) then
			write(*,*)'!!!!!!!! ERROR !!!!!!!!'
			write(*,*)text
			write(*,*)'Simulation aborted.'
			stop
		else
			write(*,*)'Something went wrong! No information on type of error.'
			stop
		end if
		
	end subroutine panic
	
	
end module tests