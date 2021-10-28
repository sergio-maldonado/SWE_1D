!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines related to time integration: compute dt, find new values,	  !
! includes bed friction (solution of ODE)								  !
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: February 2018                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module time_integrate
	implicit none
		
contains
	
	!========= Selects criterion to get dt ===================
	subroutine select_dt(selector,n,dx,courant,hl,hr,ul,ur,dtmin)
		implicit none
		integer, intent(in) :: n,selector
		real*8, dimension(:), intent(in) :: dx
		real*8, dimension(0:n), intent(in) :: hl,hr,ul,ur
		real*8, intent(in) :: courant
		real*8, intent(out) :: dtmin
		
		select case (selector)
			case (1)
			call min_dt_1(n,dx,courant,hl,hr,ul,ur,dtmin)
			
			case (2)
			call min_dt_2(n,dx,courant,hl,hr,ul,ur,dtmin)
			
			case default
			write(*,*)'Unexpected error in time_integrate.f95 (no valid criterion for dt)'
			write(*,*)'Simulation aborted.'
			stop
		
		end select
	
	end subroutine select_dt
	
	
	
	!========= Estimates minimum time step required, method 1/2 ===================
	subroutine min_dt_1(n,dx,courant,hl,hr,ul,ur,dtmin)
	! less conservative estimate, if all domain is wet
	! eq. (9.21) of Toro p. 158
		use parameters
		implicit none
		integer, intent(in) :: n
		real*8, dimension(0:n), intent(in) :: hl,hr,ul,ur
		real*8, dimension(:), intent(in) :: dx
		real*8, intent(in) :: courant
		real*8, intent(out) :: dtmin
		real*8 :: dt_prov,dxmin,smax
		real*8, dimension(0:n) :: al,ar,s1,s2
		integer :: ii
			
		al = dsqrt(g*hl)
		ar = dsqrt(g*hr)
			
		s1 = abs(ur) + ar
		s2 = abs(ul) + al
		
		! include only innner cell interfaces (no boundaries)	
		dtmin = 1.0d3 !a very large number
		do ii=1,n-1
			smax = max( s1(ii),s2(ii) )
			dxmin = min( dx(ii),dx(ii+1) )	
		
			if (dabs(smax) > zero) then	
				dt_prov = courant*dxmin/smax
				if (dt_prov < dtmin) dtmin = dt_prov
			end if

		end do
			
	end subroutine min_dt_1
	
	
			
	!========= Estimates minimum time step required, method 2/2 ===================
	subroutine min_dt_2(n,dx,courant,hl,hr,ul,ur,dtmin)
	! conservative estimate, considering (fast) dry front propagation
	! eq. (10.77) of Toro p. 195
		use parameters
		implicit none
		integer, intent(in) :: n
		real*8, dimension(0:n), intent(in) :: hl,hr,ul,ur
		real*8, dimension(:), intent(in) :: dx
		real*8, intent(in) :: courant
		real*8, intent(out) :: dtmin
		real*8 :: dt_prov,dxmin,smax
		real*8, dimension(0:n) :: al,ar,s1,s2
		integer :: ii
			
		al = dsqrt(g*hl)
		ar = dsqrt(g*hr)
			
		s1 = abs(ur) + 2.0d0*ar
		s2 = abs(ul) + 2.0d0*al
		
		! include only innner cell interfaces (no boundaries)	
		dtmin = 1.0d3 !a very large number
		do ii=1,n-1
			smax = max( s1(ii),s2(ii) )
			dxmin = min( dx(ii),dx(ii+1) )	
		
			if (dabs(smax) > zero) then	
				dt_prov = courant*dxmin/smax
				if (dt_prov < dtmin) dtmin = dt_prov
			end if

		end do
			
	end subroutine min_dt_2
		
		
		
	!========= Predictor step at dt/2 ===================	
	subroutine predictor_step(n,dt,dx,q_e,q_w,hr,hl,ur,ul,is_wet,q_e_bar,q_w_bar)
	! Based on method described in Toro p.207, eq. (11.24)
		use parameters
		use equations
		implicit none
		integer, intent(in) :: n
		real*8, intent(in) :: dt
		real*8, dimension(:), intent(in) :: dx
		real*8, dimension(:,:), intent(in) :: q_e,q_w
		real*8, dimension(0:n), intent(in) :: hr,hl,ur,ul
		integer, dimension(0:n), intent(in) :: is_wet
		real*8, dimension(:,:), intent(inout) :: q_e_bar,q_w_bar
		real*8, dimension(m,n) :: fluxes_l,fluxes_r
		integer :: k
		
		! obtain fluxes as function of left values:	
		call fluxes(hr(0:n-1),ur(0:n-1),is_wet(0:n-1),fluxes_l(:,:))
		! obtain fluxes as function of right values:
		call fluxes(hl(1:n),ul(1:n),is_wet(1:n),fluxes_r(:,:))
		
		! update boundary values:
		do k=1,m	
			
			q_w_bar(k,:) = q_w(k,:) + 0.5d0*dt*(fluxes_l(k,:) - fluxes_r(k,:))/dx(:)
			q_e_bar(k,:) = q_e(k,:) + 0.5d0*dt*(fluxes_l(k,:) - fluxes_r(k,:))/dx(:)
			
		end do
		
	end subroutine predictor_step
	
	
	
	!========= Predictor step at dt/2 with Valiani fluxes ===================	
	subroutine predictor_step_valiani(n,dt,dx,q_e,q_w,hr,hl,ur,ul,bed,is_wet,&
									&q_e_bar,q_w_bar,flux_slope_e,flux_slope_w,eta_bar)								
	! Based on method described in Toro p.207, eq. (11.24), but including slope-related fluxes as in Valiani (2006)
		use parameters
		use equations
		implicit none
		integer, intent(in) :: n
		real*8, intent(in) :: dt
		real*8, dimension(:), intent(in) :: dx
		real*8, dimension(:,:), intent(in) :: q_e,q_w
		real*8, dimension(0:n), intent(in) :: hr,hl,ur,ul
		real*8, dimension(:), intent(in) :: bed
		integer, dimension(0:n), intent(in) :: is_wet
		real*8, dimension(:,:), intent(inout) :: q_e_bar,q_w_bar
		real*8, dimension(:,:), intent(in) :: flux_slope_e,flux_slope_w
		real*8, dimension(:), intent(out) :: eta_bar
		real*8, dimension(m,n) :: fluxes_l,fluxes_r
		integer :: k
		
		! obtain fluxes as function of left values:	
		call fluxes(hr(0:n-1),ur(0:n-1),is_wet(0:n-1),fluxes_l(:,:))
		! obtain fluxes as function of right values:
		call fluxes(hl(1:n),ul(1:n),is_wet(1:n),fluxes_r(:,:))
		
		! update boundary values:
		do k=1,m	
			
			q_w_bar(k,:) = q_w(k,:) + 0.5d0*dt*(fluxes_l(k,:) - fluxes_r(k,:) + flux_slope_e(k,:) - flux_slope_w(k,:) )/dx(:)
			q_e_bar(k,:) = q_e(k,:) + 0.5d0*dt*(fluxes_l(k,:) - fluxes_r(k,:) + flux_slope_e(k,:) - flux_slope_w(k,:) )/dx(:)
			
		end do
		
		eta_bar(:) = 0.5d0*(q_w_bar(1,:) + q_e_bar(1,:)) + bed(:)
		
	end subroutine predictor_step_valiani	
		
		
		
	!========= Update variables at t+dt ===================
	subroutine new_value(n,dt,dx,flux,q)	
		use parameters
		implicit none
		integer, intent(in) :: n
		real*8, intent(in) :: dt
		real*8, dimension(:), intent(in) :: dx
		real*8, dimension(:,:), intent(in) :: flux
		real*8, dimension(m,-1:n+2), intent(inout) :: q
		integer :: i,k
		
		do k=1,m
			do i=2,n-1
				
				q(k,i) = q(k,i) - dt*(flux(k,i) - flux(k,i-1))/dx(i)
				
			end do
		end do
		
	end subroutine
	
	
	
	!========= Update variables at t+dt ===================
	!	It includes slope term following Valiani (2006)
	subroutine new_value_valiani(n,dt,dx,flux,q,q_pred,flux_slope_e,flux_slope_w)	!MODI
		use parameters
		implicit none
		integer, intent(in) :: n
		real*8, intent(in) :: dt
		real*8, dimension(:), intent(in) :: dx
		real*8, dimension(:,:), intent(in) :: flux,flux_slope_e,flux_slope_w
		real*8, dimension(m,-1:n+2), intent(in) :: q
		real*8, dimension(m,-1:n+2), intent(out) :: q_pred
		integer :: i,k
		
		do k=1,m
			do i=2,n-1
				
				q_pred(k,i) = q(k,i) - dt*(flux(k,i) - flux(k,i-1) - flux_slope_e(k,i) + flux_slope_w(k,i) )/dx(i)
				
			end do
		end do
		
	end subroutine
	
	
	
	!========= Includes bed friction w semi-implicit scheme ===================
	! Definition of variables is q1 = h; q2 = u h
	! Uses Manning's friction
	! Two options available for updating
	subroutine bed_friction(n,dt,hlim,mann,q_pred,q)
		use parameters
		implicit none
		integer, intent(in) :: n
		real*8, intent(in) :: dt,hlim,mann
		real*8, dimension(m,-1:n+2), intent(inout) :: q
		real*8, dimension(m,-1:n+2), intent(in) :: q_pred
		real*8 :: h,vel
		integer :: i,method
		
		!Select method to update:
		! method 1 uses vars. at n and pred. Gives right results for uniform flow test.  USE 1.
		! method 2 uses only predicted vars. Gives apporx. results for uniform flow test.
		! method 3 uses only vars. at n. Very similar to 1.
		method = 1	
		
		select case (method)
		
		case (1)
			
			do i=2,n-1
				if (q(1,i) > hlim) then
					h = q(1,i)
					vel = q(2,i)/h
					q(2,i) = q_pred(2,i)/( 1.0d0 + g*dt*mann*mann*dabs(vel)/(h**(4.0d0/3.0d0)) )
				
				elseif (q_pred(1,i) > hlim) then
					h = q_pred(1,i)
					vel = q_pred(2,i)/h
					q(2,i) = q_pred(2,i)/( 1.0d0 + g*dt*mann*mann*dabs(vel)/(h**(4.0d0/3.0d0)) )
					else
						q(2,i) = q_pred(2,i)
				end if
					q(1,i) = q_pred(1,i)
			end do
		
		case (2)
			
			do i=2,n-1
				if (q_pred(1,i) > hlim) then
					h = q_pred(1,i)
					vel = q_pred(2,i)/h
					q(2,i) = q_pred(2,i)/( 1.0d0 + g*dt*mann*mann*dabs(vel)/(h**(4.0d0/3.0d0)) )

					else
						q(2,i) = q_pred(2,i)
				end if
					q(1,i) = q_pred(1,i)
			end do
		
			case (3)
			
				do i=2,n-1
					if (q(1,i) > hlim) then
						h = q(1,i)
						vel = q(2,i)/h
						q(2,i) = q_pred(2,i)/( 1.0d0 + g*dt*mann*mann*dabs(vel)/(h**(4.0d0/3.0d0)) )
		
						else
							q(2,i) = q_pred(2,i)
					end if
						q(1,i) = q_pred(1,i)
				end do
		
		end select
		
		
	end subroutine		
	
end module time_integrate