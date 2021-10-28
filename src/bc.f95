!Subroutines related to Boundary conditions
!	Variables defined as q = [h,uh]
!	3 ghost cells added to each boundary

module bc
	implicit none
	
contains
	
	!========= Reflective boundary conditions ====================
	! type 1 - Imposes wall
	!			*not very well tested
	subroutine bc_reflective(n,q,we)
		implicit none
		integer, intent(in) :: n
		real*8, dimension(2,-1:n+2), intent(inout) :: q
		character*1, intent(in) :: we
		
		select case (we)
		
			case("w")
			!q(1,1) = q(1,2)
			q(1,1) = q(1,1)
			!q(1,0) = q(1,3)
			q(1,0) = q(1,0)
			!q(1,-1) = q(1,4)
			q(1,-1) = q(1,-1)
		
			q(2,1) = -q(2,2)
			!q(2,1) = q(2,1)
			q(2,0) = -q(2,3)
			!q(2,0) = q(2,0)
			q(2,-1) = -q(2,4)
			!q(2,-1) = q(2,-1)
		
			case("e")
			q(1,n) = q(1,n-1)		
			q(1,n+1) = q(1,n-2)
			q(1,n+2) = q(1,n-3)
		
			q(2,n) = -q(2,n-1)
			q(2,n+1) = -q(2,n-2)
			q(2,n+2) = -q(2,n-3)
		
			case default
			write(*,*)'Unexpected error in bc.f95'
			write(*,*)'Simulation aborted.'
			stop
			
		end select
			
	end subroutine bc_reflective
	
	
	
	!========= Transmissive boundaries 1 ====================
	! Type 2 - fixed prescribed h, q is linearly extrapolated
	! 			*works well for uniform flow test
	subroutine bc_transm1(n,q,we)
		implicit none
		integer, intent(in) :: n
		real*8, dimension(2,-1:n+2), intent(inout) :: q
		character*1, intent(in) :: we
		
		select case (we)
		
			case("w")
			q(1,1) = q(1,1)
			q(1,0) = q(1,0)
			q(1,-1) = q(1,-1)
			
			q(2,1) = 2.0d0*q(2,2) - q(2,3)
			q(2,0) = 2.0d0*q(2,1) - q(2,2)
			q(2,-1) = 2.0d0*q(2,0) - q(2,1)
		
			case("e")
			
			q(1,n) = q(1,n)
			q(1,n+1) = q(1,n+1)
			q(1,n+2) = q(1,n+2)
			
			q(2,n) = 2.0d0*q(2,n-1) - q(2,n-2)
			q(2,n+1) = 2.0d0*q(2,n) - q(2,n-1)
			q(2,n+2) = 2.0d0*q(2,n+1) - q(2,n)
		
			case default
			write(*,*)'Unexpected error in bc.f95'
			write(*,*)'Simulation aborted.'
			stop
			
		end select

	end subroutine bc_transm1
	
	
	
	!========= For still water test ====================
	! Type 3 - fixed depth, imposed zero discharge
	!			*works well for still pond test
	subroutine bc_still(n,q,we)
		implicit none
		integer, intent(in) :: n
		real*8, dimension(2,-1:n+2), intent(inout) :: q
		character*1, intent(in) :: we
		
		select case (we)
		
			case("w")
			q(1,1) = q(1,1)
			q(1,0) = q(1,0)
			q(1,-1) = q(1,-1)
			
			q(2,1) = 0.0d0
			q(2,0) = 0.0d0
			q(2,-1) = 0.0d0
		
			case("e")
			
			q(1,n) = q(1,n)
			q(1,n+1) = q(1,n+1)
			q(1,n+2) = q(1,n+2)
			
			q(2,n) = 0.0d0
			q(2,n+1) = 0.0d0
			q(2,n+2) = 0.0d0
		
			case default
			write(*,*)'Unexpected error in bc.f95'
			write(*,*)'Simulation aborted.'
			stop
			
		end select
		
	end subroutine bc_still
	
	
	
	!========= For the hump case ====================
	! Type 4 - fixed prescribed h downstream and q upstream
	! 			*
	subroutine bc_hump(n,q,we)
		implicit none
		integer, intent(in) :: n
		real*8, dimension(2,-1:n+2), intent(inout) :: q
		character*1, intent(in) :: we
		
		select case (we)
		
			case("w")
			q(1,1) = q(1,2)
			q(1,0) = q(1,1)
			q(1,-1) = q(1,0)
			
			q(2,1) = q(2,1)
			q(2,0) = q(2,0)
			q(2,-1) = q(2,-1)
		
			case("e")
			
			q(1,n) = q(1,n)
			q(1,n+1) = q(1,n+1)
			q(1,n+2) = q(1,n+2)
			
			q(2,n) = 2.0d0*q(2,n-1) - q(2,n-2)
			q(2,n+1) = 2.0d0*q(2,n) - q(2,n-1)
			q(2,n+2) = 2.0d0*q(2,n+1) - q(2,n)
		
			case default
			write(*,*)'Unexpected error in bc.f95'
			write(*,*)'Simulation aborted.'
			stop
			
		end select

	end subroutine bc_hump
	
	
	
	!========= Transmissive Riemann fixed h  ====================
	! Type 5 - conserves Riemann invariant, uses characteristics, prescribed h
	! 			*Also seems to work well for uniform flow test
	subroutine bc_invariant_h(n,q,q_old,we)
		use parameters
		implicit none
		integer, intent(in) :: n
		real*8, dimension(2,-1:n+2), intent(inout) :: q
		real*8, dimension(2,-1:n+2), intent(in) :: q_old
		character*1, intent(in) :: we
		real*8 :: inv
		
		select case (we)
		
			case("w")
			q(1,1) = q(1,1)
			q(1,0) = q(1,0)
			q(1,-1) = q(1,-1)
			
			inv = q_old(2,2) - 2.0d0*dsqrt(g)*q_old(1,2)**1.5d0
			q(2,1) = inv + 2.0d0*dsqrt(g)*q(1,1)**1.5d0
			inv = q_old(2,1) - 2.0d0*dsqrt(g)*q_old(1,1)**1.5d0
			q(2,0) = inv + 2.0d0*dsqrt(g)*q(1,0)**1.5d0
			inv = q_old(2,0) - 2.0d0*dsqrt(g)*q_old(1,0)**1.5d0
			q(2,-1) = inv + 2.0d0*dsqrt(g)*q(1,-1)**1.5d0
			
			case("e")
			
			q(1,n) = q(1,n)
			q(1,n+1) = q(1,n+1)
			q(1,n+2) = q(1,n+2)
			
			inv = q_old(2,n-1) + 2.0d0*dsqrt(g)*q_old(1,n-1)**1.5d0
			q(2,n) = inv - 2.0d0*dsqrt(g)*q(1,n)**1.5d0
			inv = q_old(2,n) + 2.0d0*dsqrt(g)*q_old(1,n)**1.5d0
			q(2,n+1) = inv - 2.0d0*dsqrt(g)*q(1,n+1)**1.5d0
			inv = q_old(2,n+1) + 2.0d0*dsqrt(g)*q_old(1,n+1)**1.5d0
			q(2,n+2) = inv - 2.0d0*dsqrt(g)*q(1,n+2)**1.5d0
		
			case default
			write(*,*)'Unexpected error in bc.f95'
			write(*,*)'Simulation aborted.'
			stop
			
		end select
		
	end subroutine bc_invariant_h
	
	
	
	!========= Transmissive Riemann fixed q  ====================
	! Type 5 - conserves Riemann invariant, uses characteristics, prescribed q
	! 			*Also seems to work well for uniform flow test
	subroutine bc_invariant_q(n,q,q_old,we)
		use parameters
		implicit none
		integer, intent(in) :: n
		real*8, dimension(2,-1:n+2), intent(inout) :: q
		real*8, dimension(2,-1:n+2), intent(in) :: q_old
		character*1, intent(in) :: we
		real*8 :: inv
		
		select case (we)
		
			case("w")
			
			q(2,1) = q(2,1)
			q(2,0) = q(2,0)
			q(2,-1) = q(2,-1)
			
			inv = 2.0d0*dsqrt(g)*q_old(1,2)**1.5d0 - q_old(2,2)
			q(1,1) = ( (inv + q(2,1))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
			inv = 2.0d0*dsqrt(g)*q_old(1,1)**1.5d0 - q_old(2,1)
			q(1,0) = ( (inv + q(2,0))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
			inv = 2.0d0*dsqrt(g)*q_old(1,0)**1.5d0 - q_old(2,0)
			q(1,-1) = ( (inv + q(2,-1))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
			
			case("e")
			
			q(2,n) = q(2,n)
			q(2,n+1) = q(2,n+1)
			q(2,n+2) = q(2,n+2)
			
			inv = 2.0d0*dsqrt(g)*q_old(1,n-1)**1.5d0 + q_old(2,n-1)
			q(1,n) = ( (inv - q(2,n))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
			inv = 2.0d0*dsqrt(g)*q_old(1,n)**1.5d0 + q_old(2,n)
			q(1,n+1) = ( (inv - q(2,n+1))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
			inv = 2.0d0*dsqrt(g)*q_old(1,n+1)**1.5d0 + q_old(2,n+1)
			q(1,n+2) = ( (inv - q(2,n+2))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
		
			case default
			write(*,*)'Unexpected error in bc.f95'
			write(*,*)'Simulation aborted.'
			stop
			
		end select
		
	end subroutine bc_invariant_q
		
! 	!========= For uniform flow gravity driven ====================
! 	! fixed prescribed h, q is linearly extrapolated
! 	subroutine bc_uniform(n,q)
! 		implicit none
! 		integer, intent(in) :: n
! 		real*8, dimension(2,-1:n+2), intent(inout) :: q
! 		real*8 :: manning,slope,q_unif
! 		!Variables defined as q = [h,uh]
! 		!Variable 1:
! 		q(1,-1) = q(1,-1)
! 		q(1,n+2) = q(1,n+2)
!
! 		q(1,0) = q(1,0)
! 		q(1,n+1) = q(1,n+1)
!
! 		q(1,1) = q(1,1)
! 		q(1,n) = q(1,n)
!
! 		!Variable 2:
! 		q(2,1) = 2.0d0*q(2,2) - q(2,3)
! 		q(2,n) = 2.0d0*q(2,n-1) - q(2,n-2)
!
! 		q(2,0) = 2.0d0*q(2,1) - q(2,2)
! 		q(2,n+1) = 2.0d0*q(2,n) - q(2,n-1)
!
! 		q(2,-1) = 2.0d0*q(2,0) - q(2,1)
! 		q(2,n+2) = 2.0d0*q(2,n+1) - q(2,n)
		
! TEMP:	Also converges to right solution but with intial 'shock'
! 		manning = 0.04d0
! 		slope = 1.0d-3
! 		q_unif = (1.0d0/manning)*dsqrt(slope)*q(1,1)**(5.0d0/3.0d0)
! 		q(2,1) = q_unif
! 		q(2,n) = q_unif
!
! 		q(2,0) = q_unif
! 		q(2,n+1) = q_unif
!
! 		q(2,-1) = q_unif
! 		q(2,n+2) = q_unif

!	end subroutine
	
	
	
	!========= time series ===================
	!interpolates one variable, the other one computed from invariant
	subroutine bc_timeseries(n,q,q_old,bd,ts,variable,tsim)
		use parameters
		use tests
		implicit none
		integer, intent(in) :: n
		real*8, intent(in) :: tsim
		real*8, dimension(:,:), intent(in) :: ts
		real*8, dimension(2,-1:n+2), intent(inout) :: q
		real*8, dimension(2,-1:n+2), intent(in) :: q_old
		character*1, intent(in) :: bd,variable
		integer :: var,is_west,i,ind
		logical :: scan
		real*8 :: inv,t0,t1,y0,y1,y
		
		if (variable == 'h') then
			var = 1
			else if (variable == 'q') then
				var = 2
			else
				call panic(2,"bc.f95")
		end if
		
		if (bd == 'w') then
			is_west = 1
			else if (bd == 'e') then
				is_west = 0
			else
				call panic(2,"bc.f95")
		end if
		
		scan = .true.
		i=1
		do while (scan)
			if (ts(1,i) <= tsim .and. ts(1,i+1) >= tsim) scan=.false.
			i=i+1
		end do
		
		! linear interpolation:
		t0 = ts(1,i-1)
		t1 = ts(1,i)
		y0 = ts(2,i-1)
		y1 = ts(2,i)
		
		y = y0 + (tsim - t0)*(y1-y0)/(t1-t0)
		
		! assigns interpolated value to selected variable:
		do i=0,2 
			ind = n+i - (n+1)*is_west
			q(var,ind) = y
		end do
		
		
		! assigns the other variable from Riemann invariant
		select case (variable)
		
			case("h")
			
				select case (bd)
				case("w")
				inv = q_old(2,2) - 2.0d0*dsqrt(g)*q_old(1,2)**1.5d0
				q(2,1) = inv + 2.0d0*dsqrt(g)*q(1,1)**1.5d0
				inv = q_old(2,1) - 2.0d0*dsqrt(g)*q_old(1,1)**1.5d0
				q(2,0) = inv + 2.0d0*dsqrt(g)*q(1,0)**1.5d0
				inv = q_old(2,0) - 2.0d0*dsqrt(g)*q_old(1,0)**1.5d0
				q(2,-1) = inv + 2.0d0*dsqrt(g)*q(1,-1)**1.5d0
			
				case("e")
				inv = q_old(2,n-1) + 2.0d0*dsqrt(g)*q_old(1,n-1)**1.5d0
				q(2,n) = inv - 2.0d0*dsqrt(g)*q(1,n)**1.5d0
				inv = q_old(2,n) + 2.0d0*dsqrt(g)*q_old(1,n)**1.5d0
				q(2,n+1) = inv - 2.0d0*dsqrt(g)*q(1,n+1)**1.5d0
				inv = q_old(2,n+1) + 2.0d0*dsqrt(g)*q_old(1,n+1)**1.5d0
				q(2,n+2) = inv - 2.0d0*dsqrt(g)*q(1,n+2)**1.5d0
				end select
			
			case("q")
			
				select case (bd)
				case("w")
				inv = 2.0d0*dsqrt(g)*q_old(1,2)**1.5d0 - q_old(2,2)
				q(1,1) = ( (inv + q(2,1))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
				inv = 2.0d0*dsqrt(g)*q_old(1,1)**1.5d0 - q_old(2,1)
				q(1,0) = ( (inv + q(2,0))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
				inv = 2.0d0*dsqrt(g)*q_old(1,0)**1.5d0 - q_old(2,0)
				q(1,-1) = ( (inv + q(2,-1))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
			
				case("e")
				inv = 2.0d0*dsqrt(g)*q_old(1,n-1)**1.5d0 + q_old(2,n-1)
				q(1,n) = ( (inv - q(2,n))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
				inv = 2.0d0*dsqrt(g)*q_old(1,n)**1.5d0 + q_old(2,n)
				q(1,n+1) = ( (inv - q(2,n+1))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
				inv = 2.0d0*dsqrt(g)*q_old(1,n+1)**1.5d0 + q_old(2,n+1)
				q(1,n+2) = ( (inv - q(2,n+2))/2.0d0/dsqrt(g) )**(2.0d0/3.0d0)
				end select
		
			case default
			write(*,*)'Unexpected error in bc.f95'
			write(*,*)'Simulation aborted.'
			stop
		
		end select
		
	
	end subroutine
	
	
	
	!========= Wraps all previous subroutines ===================
	subroutine boundaries()
		use variables
		implicit none
		character*1 :: bound
		
		!West boundary:
		select case (bc_w)
		
			case (1)
				bound = "w"
				call bc_reflective(n,q,bound)
			case (2)
				bound = "w"
				call bc_transm1(n,q,bound)
			case (3)
				bound = "w"
				call bc_still(n,q,bound)
			case (4)
				bound = "w"
				call bc_hump(n,q,bound)
			case (5)
				bound = "w"
				call bc_invariant_h(n,q,q_old,bound)
			case (6)
				bound = "w"
				call bc_invariant_q(n,q,q_old,bound)
			case (10)
				bound = "w"
				call bc_timeseries(n,q,q_old,bound,ts_w,ts_var_w,tsim)
			case default
				write(*,*)'Unexpected error in bc.f95 (no valid bc type)'
				write(*,*)'Simulation aborted.'
				stop
		
		end select 
		
		!East boundary:
		select case (bc_e)
		
			case (1)
				bound = "e"
				call bc_reflective(n,q,bound)
			case (2)
				bound = "e"
				call bc_transm1(n,q,bound)
			case (3)
				bound = "e"
				call bc_still(n,q,bound)
			case (4)
				bound = "e"
				call bc_hump(n,q,bound)	
			case (5)
				bound = "e"
				call bc_invariant_h(n,q,q_old,bound)
			case (6)
				bound = "e"
				call bc_invariant_q(n,q,q_old,bound)
			case (10)
				bound = "e"
				call bc_timeseries(n,q,q_old,bound,ts_e,ts_var_e,tsim)
			case default
				write(*,*)'Unexpected error in bc.f95 (no valid bc type)'
				write(*,*)'Simulation aborted.'
				stop
		
		end select 
		
		
	end subroutine
		
	
end module bc