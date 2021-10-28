!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different tests for each subroutine									  !
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: January 2018                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module tests_bugs
	implicit none
	
	
contains
	
	!========= writes files with reconstructed E and W data ===================
	subroutine test_data_recons(n,q_e,q_w)
		implicit none
		real*8, dimension(2,n), intent(in) :: q_e,q_w
		integer, intent(in) :: n
		integer :: i
		
		open(10,file = 'test_data_recons.dat',status='unknown')
		write(10,*)'% i, qE_1, qE_2, qW_1, qW_2'
		do i=1,n
			write(10,'(i3,4(2x,es10.3))') i,q_e(1,i),q_e(2,i),q_w(1,i),q_w(2,i)
		end do
		close(10)
		
	end subroutine test_data_recons 
	
	
	!========= writes files with reconstructed E and W bed ===================
	subroutine test_bed_recons(n,bed_e,bed_w)
		implicit none
		real*8, dimension(:), intent(in) :: bed_e,bed_w
		integer, intent(in) :: n
		integer :: i
		
		open(11,file = 'test_bed_recons.dat',status='unknown')
		write(11,*)'% i, zb_E, zb_W'
		do i=1,n
			write(11,'(i3,2(2x,es10.3))') i,bed_e(i),bed_w(i)
		end do
		close(11)
		
	end subroutine test_bed_recons
	
	
	!========= writes files with reconstructed left and right states ===================
	subroutine test_left_right(n,hl,hr,ul,ur)
		implicit none
		real*8, dimension(:), intent(in) :: hr,hl,ur,ul
		integer, intent(in) :: n
		integer :: ii
		
		open(12,file = 'test_left_right.dat',status='unknown')
		write(12,*)'% ii, hL, hR, uL, uR'
		do ii=0,n
			write(12,'(i3,4(2x,es10.3))') ii,hl(ii),hr(ii),ul(ii),ur(ii)
		end do
		close(12)
		
	end subroutine test_left_right
	
	
	!========= writes files with reconstructed star region variables ===================
	subroutine test_star_vars(n,hst,ust)
		implicit none
		real*8, dimension(:), intent(in) :: hst,ust
		integer, intent(in) :: n
		integer :: ii
		
		open(13,file = 'test_star_vars.dat',status='unknown')
		write(13,*)'% ii, h*, u*'
		do ii=1,n-1
			write(13,'(i3,2(2x,es10.3))') ii,hst(ii),ust(ii)
		end do
		close(13)
		
	end subroutine test_star_vars
	
	
	!========= writes files with speeds estimates at intercells ===================
	subroutine test_waves_speeds(n,sl,sr,sst)
		implicit none
		real*8, dimension(:), intent(in) :: sl,sr,sst
		integer, intent(in) :: n
		integer :: ii
		
		open(14,file = 'test_waves_speeds.dat',status='unknown')
		write(14,*)'% ii, sL, sR, s_*'
		do ii=1,n-1
			write(14,'(i3,3(2x,es10.3))') ii,sl(ii),sr(ii),sst(ii)
		end do
		close(14)
		
	end subroutine test_waves_speeds
	
	
	!========= writes files with fluxes ===================
	subroutine test_fluxes(f)
		implicit none
		real*8, dimension(:,:), intent(in) :: f
		!integer, intent(in) :: n
		integer :: i,n
		
		n = size(f(1,:))
		
		open(15,file = 'test_fluxes.dat',status='unknown')
		write(15,*)'% i, flux1,flux2'
		do i=1,n
			write(15,'(i3,2(2x,es10.3))') i,f(1,i),f(2,i)
		end do
		close(15)
		
	end subroutine test_fluxes
	
	
	!========= writes files with predicted boundary values ===================
	subroutine test_predicted_values(q_e_pred,q_w_pred)
		implicit none
		real*8, dimension(:,:), intent(in) :: q_e_pred,q_w_pred
		integer :: i,n
		
		n = size(q_e_pred(1,:))
		
		open(16,file = 'test_predicted_vals.dat',status='unknown')
		write(16,*)'% i, pred q_w(1), q_w(2), q_e(1), q_e(2)'
		do i=1,n
			write(16,'(i3,4(2x,es10.3))') i,q_w_pred(1,i),q_w_pred(2,i),q_e_pred(1,i),q_e_pred(2,i)
		end do
		close(16)
		
	end subroutine test_predicted_values
	
	
	!========= writes files with updated variables ===================
	subroutine test_new_value(q)
		implicit none
		real*8, dimension(:,:), intent(in) :: q
		integer :: i,n
		
		n = size(q(1,:))
		
		open(17,file = 'test_new_values.dat',status='unknown')
		write(17,*)'% i, q1, q2'
		do i=1,n
			write(17,'(i3,2(2x,es10.3))') i,q(1,i),q(2,i)
		end do
		close(17)
		
	end subroutine test_new_value
	
	
	!========= select tests to be performed ===================	
	subroutine test_select
		! (un)comment lines as required
		use variables
		implicit none
		
		!TESTS:
		!call test_bed_recons(n,zb_e,zb_w)
		!call test_data_recons(n,q_e,q_w)
		!call test_left_right(n,hl,hr,ul,ur)
		!call test_star_vars(n,hstar,ustar)
		!call test_waves_speeds(n,sl,sr,sstar)
		call test_fluxes(f)
		!call test_predicted_values(q_e_bar,q_w_bar)
		!call test_new_value(q)
		
	end subroutine test_select
	
end module tests_bugs