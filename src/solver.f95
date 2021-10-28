!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Main body of the program, contains data reconstruction, Riemann		  !
!	solver, slope-fluxes (as in Valiani)								  ! 
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: February 2018                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module solver
	implicit none
	
	
contains
	
	!========= Reconstructs data within cells ===================
	! Note that beta = 0 implies piecewise constant reconstruction
	subroutine data_recons(n,beta,bed,bed_e,bed_w,w,w_e,w_w,option)
		use parameters
		implicit none
		real*8, intent(in) :: beta
		real*8 :: delta_i_left,delta_i_right,delta_bar
		real*8, dimension(m,-1:n+2), intent(in) :: w
		real*8, dimension(:,:), intent(out) :: w_e,w_w
		real*8, dimension(n), intent(in) :: bed,bed_e,bed_w
		real*8, dimension(m,n) :: wp,wp_e,wp_w
		integer :: i,k
		integer, intent(in) :: n,option
		
		select case (option)
			case (1)	!Slope limiter applied to eta (free surface level)
			wp(1,1:n) = bed(1:n) + w(1,1:n)
			case (2)	!Slope limiter applied to h
			wp(1,1:n) = w(1,1:n)
		end select
		
		wp(2,1:n) = w(2,1:n)
	
		do k=1,m
			do i=2,n-1

				delta_i_left = wp(k,i) - wp(k,i-1)
				delta_i_right = wp(k,i+1) - wp(k,i)
			
			
				if(delta_i_right > zero) then
		    		delta_bar = dmax1(0.0d0, dmin1(beta*delta_i_left,delta_i_right), dmin1(delta_i_left,beta*delta_i_right))
					else if (delta_i_right < -zero) then
		    		delta_bar = dmin1(0.0d0, dmax1(beta*delta_i_left,delta_i_right), dmax1(delta_i_left,beta*delta_i_right))
				else
					delta_bar = 0.0d0
				end if

				wp_e(k,i) = wp(k,i) + delta_bar/2.0d0
				wp_w(k,i) = wp(k,i) - delta_bar/2.0d0

			end do
			wp_e(k,1) = wp_w(k,2)
			wp_e(k,n) = wp(k,n)
			wp_w(k,1) = wp(k,1)
			wp_w(k,n) = wp_e(k,n-1)
		end do
		
		select case (option)
			case (1)
			w_e(1,:) = wp_e(1,:) - bed_e(:)
			w_w(1,:) = wp_w(1,:) - bed_w(:)
			case (2)
			w_e(1,:) = wp_e(1,:)
			w_w(1,:) = wp_w(1,:)
		end select
	
		w_e(2,:) = wp_e(2,:)
		w_w(2,:) = wp_w(2,:)
		
	end subroutine data_recons
	
	
	
	!========= Reconstructs bed (linearly) ===================
	subroutine bed_recons(n,bed,bed_e,bed_w)	
		implicit none
		real*8, dimension(:), intent(in) :: bed
		real*8, dimension(:), intent(inout) :: bed_e,bed_w
		integer, intent(in) :: n
		
		!in West faces:
		bed_w(2 : n) = (bed(1 : n-1) + bed(2 : n))/2.0d0
		bed_w(1) = 2.0d0*bed_w(2) - bed_w(3)
		!East faces equal to West faces of adjacet cell (no discontinuities in bed)
		bed_e(1 : n-1) = bed_w(2 : n)
		bed_e(n) = 2.0d0*bed_e(n-1) - bed_e(n-2)
		
	end subroutine bed_recons
	
	
	
	!========= Computes left and right states at interfaces ===================
	!		Also determines dry-dry interfaces
	subroutine left_right_states(n,hlim,q_e,q_w,hl,hr,ul,ur,is_wet)
		implicit none
		real*8, intent(in) :: hlim
		real*8, dimension(:,:), intent(in) :: q_e,q_w
		real*8, dimension(0:n), intent(inout) :: hr,hl,ur,ul
		integer, dimension(0:n), intent(inout) :: is_wet
		integer, intent(in) :: n
		integer :: ii
		
		!---------Attention to definition of q
		!	q1 = h; q2 = u h
		
		do ii=1,n-1	! ii is number of interfaces = number of cells +1 to include domain boundaries
			
			if (q_e(1,ii) > hlim) then
				hl(ii) = q_e(1,ii)
				ul(ii) = q_e(2,ii)/hl(ii)
			else
				hl(ii) = 0.0d0
				ul(ii) = 0.0d0
			end if
			
			if (q_w(1,ii+1) > hlim) then
				hr(ii) = q_w(1,ii+1)
				ur(ii) = q_w(2,ii+1)/hr(ii)
			else
				hr(ii) = 0.0d0
				ur(ii) = 0.0d0
			end if
			
			if (hl(ii) < hlim .and. hr(ii) < hlim) then
				is_wet(ii) = 0	!dry
			else
				is_wet(ii) = 1	!wet
			end if
		
		end do
		
		! values at boundaries
		if (q_w(1,1) > hlim) then
			hr(0) = q_w(1,1)
			ur(0) = q_w(2,1)/hr(0)
			is_wet(0) = 1
		else
			hr(0) = 0.0d0
			ur(0) = 0.0d0
			is_wet(0) = 0
		end if
		
		if (q_e(1,n) > hlim) then
			hl(n) = q_e(1,n)
			ul(n) = q_e(2,n)/hl(n)
			is_wet(n) = 1
		else
			hl(n) = 0.0d0
			ul(n) = 0.0d0
			is_wet(n) = 0
		end if
		! the following are irrelevant, but just to give them values:
		hr(n) = hl(n)
		ur(n) = ul(n)
		hl(0) = hr(0)
		ul(0) = ur(0)
		
	
	end subroutine left_right_states
	
	
	
	!========= Computes Valiani depths and fluxes at interfaces ===================
	subroutine valiani_depths(n,hlim,is_wet,bed_e,bed_w,eta,hv_e,hv_w,fv_e,fv_w)
		use equations
		use tests
		implicit none
		integer, intent(in) :: n
		real*8, intent(in) :: hlim
		integer, dimension(0:n), intent(in) :: is_wet
		real*8, dimension(:), intent(in) :: bed_e,bed_w,eta
		real*8, dimension(n), intent(out) :: hv_e,hv_w
		real*8, dimension(:,:), intent(out) :: fv_e,fv_w
		integer :: i
		
		do i=1,n
			hv_w(i) = eta(i) - bed_w(i)
			hv_e(i) = eta(i) - bed_e(i)
			
			if (hv_w(i) < hlim) hv_w(i) = 0.0d0
			if (hv_e(i) < hlim) hv_e(i) = 0.0d0
		end do
		
		call test_neg_depth(hv_w(:),'Check valiani_depths in solver.f95.')
		call test_neg_depth(hv_e(:),'Check valiani_depths in solver.f95.')
		
		call fluxes_valiani(hv_w(1:n),is_wet(0:n-1),fv_w(:,1:n))
		call fluxes_valiani(hv_e(1:n),is_wet(1:n),fv_e(:,1:n))
	
	end subroutine valiani_depths
	
	
	
	!========= Computes star region variables ===================
	subroutine star_region_variables(n,hl,hr,ul,ur,hst,ust)
		! h_* and u_* computed as in eq. (10.18) of Toro p.179
		! No risk of div by zero
		use parameters
		implicit none
		integer, intent(in) :: n
		real*8, dimension(0:n), intent(in) :: hr,hl,ur,ul
		real*8, dimension(0:n), intent(out) :: hst,ust
		real*8, dimension(0:n) :: al,ar
	
		al = dsqrt(g*hl)
		ar = dsqrt(g*hr)
		hst = (1.0d0/g)*(0.5d0*(al + ar) + 0.25d0*(ul - ur))**2.0d0
		ust = 0.5d0*(ul + ur) + al - ar
		
	end subroutine star_region_variables
	
	
	
	!========= Computes waves speeds ===================
	subroutine waves_speeds(n,hlim,is_wet,hl,hr,ul,ur,ust,hst,sl,sr,sst)
		! waves speeds estimates according to eqs. (10.22), (10.23). (10.27) of Toro pp.180-182
		use parameters
		implicit none
		real*8, intent(in) :: hlim
		integer, intent(in) :: n
		real*8, dimension(0:n), intent(in) :: hr,hl,ur,ul,ust,hst
		real*8, dimension(0:n), intent(inout) :: sl,sr,sst
		integer, dimension(0:n), intent(in) :: is_wet
		real*8, dimension(0:n) :: al,ar,ast
		integer :: ii
		real*8 :: deno
		
		al = dsqrt(g*hl)
		ar = dsqrt(g*hr)
		ast = dsqrt(g*hst)
	
		do ii=1,n-1	! values at 0 and n are irrelevant
			if (hl(ii) > hlim) then
				!sl(ii) = ul(ii) - al(ii)*ql(ii)
				sl(ii) = min(ul(ii) - al(ii), ust(ii) - ast(ii))
			else
				sl(ii) = ur(ii) - 2.0d0*ar(ii)
			end if

			if (hr(ii) > hlim) then
				!sr(ii) = ur(ii) + ar(ii)*qr(ii)
				sr(ii) = max(ur(ii) + ar(ii), ust(ii) + ast(ii))
			else
				sr(ii) = ul(ii) + 2.0d0*al(ii)
			end if

			deno = hr(ii)*(ur(ii) - sr(ii)) - hl(ii)*(ul(ii) - sl(ii))
			if (dabs(deno) > zero) then
				sst(ii) = ( sl(ii)*hr(ii)*(ur(ii) - sr(ii)) - sr(ii)*hl(ii)*(ul(ii) - sl(ii)) )/deno
			else
				sst(ii) = 0.0d0
			end if
			
			if (is_wet(ii) == 0) then
				sl(ii) = 0.0d0
				sr(ii) = 0.0d0
				sst(ii) = 0.0d0
			end if
			
		end do
		
		!!This does not work for dry bed:
! 		sl = ul - al*ql
! 		sr = ur + ar*qr
! 		sst = (sl*hr*(ur - sr) - sr*hl*(ul - sl) )/( hr*(ur - sr) - hl*(ul - sl) )
	
	end subroutine waves_speeds
	
	
	
	!========= Computes HLLC fluxes at intercells ===================
	subroutine fluxes_hllc(n,q_e,q_w,hr,hl,ur,ul,sl,sr,sst,is_wet,f_hllc)
		use parameters
		use equations
		use tests
		implicit none
		integer, intent(in) :: n
		real*8, dimension(:,:), intent(in) :: q_e,q_w
		real*8, dimension(0:n), intent(in) :: hr,hl,ur,ul
		real*8, dimension(0:n), intent(in) :: sl,sr,sst
		real*8, dimension(:,:), intent(inout) :: f_hllc
		integer, dimension(0:n), intent(in) :: is_wet
		real*8, dimension(m,1:n-1) :: f_hll,f_l,f_r,u_stl,u_str,f_asl,f_asr
		integer :: ii,k,zero_code,flux_code
		
		zero_code = 0
		flux_code = 1
		
		
		call fluxes(hl(1:n-1),ul(1:n-1),is_wet(1:n-1),f_l(:,:))
		call fluxes(hr(1:n-1),ur(1:n-1),is_wet(1:n-1),f_r(:,:))
		
		
		do ii=1,n-1
			if (is_wet(ii) == 1) then
					
				if (dabs(sl(ii) - sst(ii)) > zero) then  
					! Toro Eq. (10.26) p.182	
					u_stl(1,ii) = hl(ii)*(sl(ii) - ul(ii))/(sl(ii) - sst(ii))
					u_stl(2,ii) = sst(ii)*hl(ii)*(sl(ii) - ul(ii))/(sl(ii) - sst(ii))
				else
					u_stl(1,ii) = q_e(1,ii)
					u_stl(2,ii) = q_e(2,ii)
				end if
				
				if (dabs(sr(ii) - sst(ii)) > zero)then
					u_str(1,ii) = hr(ii)*(sr(ii) - ur(ii))/(sr(ii) - sst(ii))
					u_str(2,ii) = sst(ii)*hr(ii)*(sr(ii) - ur(ii))/(sr(ii) - sst(ii))			
				else
					u_str(1,ii) = q_w(1,ii+1)
					u_str(2,ii) = q_w(2,ii+1)
				end if
			
			else	
				u_stl(1,ii) = 0.0d0
				u_stl(2,ii) = 0.0d0
				u_str(1,ii) = 0.0d0
				u_str(2,ii) = 0.0d0	
			end if
		end do				
	
		
		
		do k=1,m
			do ii=1,n-1
				if (is_wet(ii) == 1) then
					
					f_asl(k,ii) = f_l(k,ii) + sl(ii)*(u_stl(k,ii) - q_e(k,ii))			
					f_asr(k,ii) = f_r(k,ii) + sr(ii)*(u_str(k,ii) - q_w(k,ii+1))
					
					f_hll(k,ii) = (sr(ii)*f_l(k,ii) - sl(ii)*f_r(k,ii) + &
								&sr(ii)*sl(ii)*(q_w(k,ii+1) - q_e(k,ii)))/&
								&(sr(ii) - sl(ii))
				else
					f_asl(k,ii) = 0.0d0
					f_asr(k,ii) = 0.0d0
					
					f_hll(k,ii) = 0.0d0
					
				end if
				
				
				if (sl(ii) >= 0.0d0) then
					f_hllc(k,ii) = f_l(k,ii)
					elseif (sl(ii) <= 0.0d0 .and. sst(ii) >= 0.0d0) then
						f_hllc(k,ii) = f_asl(k,ii)
						!f_hllc(k,ii) = f_hll(k,ii)
						elseif (sr(ii) >= 0.0d0 .and. sst(ii) <= 0.0d0) then
							f_hllc(k,ii) = f_asr(k,ii)
							!f_hllc(k,ii) = f_hll(k,ii)
							elseif (sr(ii) <= 0.0d0) then
								f_hllc(k,ii) = f_r(k,ii)
							else
								call panic(flux_code)
				end if
				
			end do				
		end do
		
	end subroutine fluxes_hllc
	
	
		
	!========= Wraps all previous subroutines ===================
	subroutine main_computations
		use variables
		use parameters
		use time_integrate
		use bc

		implicit none
		
		!~~~~~~~~~~For second order integration in time:
		call bed_recons(n,zb,zb_e,zb_w)
		call data_recons(n,beta,zb,zb_e,zb_w,q,q_e,q_w,sl_variable)
		call left_right_states(n,hlim,q_e,q_w,hl,hr,ul,ur,is_wet)
		call valiani_depths(n,hlim,is_wet,zb_e,zb_w,eta,hv_e,hv_w,fv_e,fv_w)
		call star_region_variables(n,hl,hr,ul,ur,hstar,ustar)
		call waves_speeds(n,hlim,is_wet,hl,hr,ul,ur,ustar,hstar,sl,sr,sstar)
		!call min_dt_2(n,dx,courant,hl,hr,ul,ur,dt)
		call select_dt(which_dt,n,dx,courant,hl,hr,ul,ur,dt)
		call predictor_step_valiani(n,dt,dx,q_e,q_w,hr,hl,ur,ul,zb,is_wet,q_e_bar,q_w_bar,fv_e,fv_w,eta_bar)
		call valiani_depths(n,hlim,is_wet,zb_e,zb_w,eta_bar,hv_e,hv_w,fv_e,fv_w)
		call left_right_states(n,hlim,q_e_bar,q_w_bar,hl,hr,ul,ur,is_wet)
		call star_region_variables(n,hl,hr,ul,ur,hstar,ustar)
		call waves_speeds(n,hlim,is_wet,hl,hr,ul,ur,ustar,hstar,sl,sr,sstar)
		call fluxes_hllc(n,q_e_bar,q_w_bar,hr,hl,ur,ul,sl,sr,sstar,is_wet,f)
		call new_value_valiani(n,dt,dx,f,q,q_as,fv_e,fv_w)	!MODI
		q_old(:,:) = q(:,:)
		call bed_friction(n,dt,hlim,manning,q_as,q)
		tsim = tsim + dt
		if (tsim > tt) tsim = tt
		call boundaries()
		
		!~~~~~~~~~~For first order integration in time: (NOT TESTED - May need modification)
! 	 	call bed_recons(n,zb,zb_e,zb_w)
! 		call data_recons(n,beta,zb,zb_e,zb_w,q,q_e,q_w)
! 		call left_right_states(n,q_e,q_w,hl,hr,ul,ur,is_wet)
! 		call valiani_depths(n,is_wet,zb_e,zb_w,eta,hv_e,hv_w,fv_e,fv_w)
! 		call star_region_variables(n,hl,hr,ul,ur,hstar,ustar)
! 		call waves_speeds(n,is_wet,hl,hr,ul,ur,ustar,hstar,sl,sr,sstar)
! 		call min_dt_2(n,dx,courant,hl,hr,ul,ur,dt)
! 		call fluxes_hllc(n,q_e,q_w,hr,hl,ur,ul,sl,sr,sstar,is_wet,f)
! 		call new_value_valiani(n,dt,dx,f,q,fv_e,fv_w)
! 		!call bc_reflective(n,q)
! 		call bc_fixed(n,q)		
! 		tsim = tsim + dt

		steps = steps + 1
		eta(1:n) = zb(1:n) + q(1,1:n)
		
	end subroutine main_computations
		
		
	
end module solver