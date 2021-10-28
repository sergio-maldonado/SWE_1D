!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reads input data and initilialise variables							  !
!                                                                         !
! Author: Sergio Maldonado (s.maldonado@soton.ac.uk) 	                  !
!                                                                         !
! Created on: January 2018                                                !
! Modified: February 2018                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module initialise
	implicit none
	
	
contains
	
	
	!========= starts screen ====================
	subroutine welcome_screen(time_start)
		implicit none
		real*8, intent(inout) :: time_start
		
		call cpu_time(time_start)
		
		write(*,*)' '
		write(*,*)'************************************************************'
		write(*,*)' '
		write(*,*)'         1D Shallow Water Equations Solver'
		write(*,*)'            (FV HLLC approx. solver)'
		write(*,*)' '
		write(*,*)'************************************************************'
		write(*,*)' '
		write(*,*)'Starting simulation...'
		
	end subroutine
	
	
	
	!========= Reads file with working directory ====================
	subroutine get_wkdir(path)
		implicit none
		character*1 :: rub
		character(len=200), intent(out) :: path
		
		write(*,*)'Getting working directory...'
		
		open(1,file = 'working_directory.inp',status='old')
		read(1,*)rub
		read(1,*)path
		close(1)
		
	end subroutine
		
		
		
	!========= Reads input file with simulation-related settings ====================
	subroutine read_settings_sim(n,file_name,format_global,total_time,twrite_global,manning,is_sensors,n_sensors,path)
		implicit none
		integer, intent(out) :: n
		character (len=30), intent(out) :: file_name
		real*8, intent(out) :: total_time,twrite_global,manning
		character*30, intent(out) :: format_global
		character(len=200), intent(in) :: path
		character*1 :: rub,dig_cells
		character*7 :: ncells,cell
		integer, intent(out) :: n_sensors
		logical, intent(out) :: is_sensors

		
		! reads input settings file:
		write(*,*)'Reading settings_sim.inp...'
		
		open(1,file = trim(path)//'settings_sim.inp',status='old')
		read(1,*)rub
		read(1,*)rub
		read(1,*)rub
		read(1,*)rub
		read(1,*)n
		read(1,*)rub
		read(1,*)file_name
		read(1,*)rub
		read(1,*)total_time
		read(1,*)rub
		read(1,*)twrite_global
		read(1,*)rub
		read(1,*)manning
		read(1,*)rub
		read(1,*)n_sensors
		close(1)
		
		if (n_sensors == 0) then
			is_sensors = .false.
		else
			is_sensors = .true.
		end if
		
		! format to print global variables:
		write(dig_cells,'(I1)')int(log10(real(n))) + 1
		ncells = '(I'//dig_cells//')'
		write(cell,ncells)n
		format_global = '('//trim(cell)//'(EN12.3,2X))'
		format_global = trim(format_global)
		
	end subroutine
	
	
	
	!========= Reads input file with numerics-related settings ====================
	subroutine read_settings_num(beta,courant,hlim,sl_variable,bc_west_type,bc_east_type,dt_crit,path)
		use tests
		implicit none
		real*8, intent(out) :: beta,courant,hlim
		integer, intent(out) :: sl_variable,bc_west_type,bc_east_type,dt_crit
		character(len=200), intent(in) :: path
		character(len=5) :: bc_west_key,bc_east_key
		character*1 :: rub
		
		! reads input numerics-related settings file:
		write(*,*)'Reading settings_num.inp...'
		
		open(1,file = trim(path)//'settings_num.inp',status='old')
		read(1,*)rub
		read(1,*)rub
		read(1,*)rub
		read(1,*)rub
		read(1,*)beta
		read(1,*)rub
		read(1,*)courant
		read(1,*)rub
		read(1,*)hlim
		read(1,*)rub
		read(1,*)sl_variable
		read(1,*)rub
		read(1,*)bc_west_key
		read(1,*)rub
		read(1,*)bc_east_key
		read(1,*)rub
		read(1,*)dt_crit
		close(1)
		
		call bc_select(bc_west_key,bc_west_type)
		call bc_select(bc_east_key,bc_east_type)
		
		! basic tests:
		if (hlim < 0.0d0) call panic(3,"Depth threshold must be positive.")
		if (courant < 0.0d0 .or. courant > 1.0d0) call panic(3,"Check Courant number (0<C<1)")
		if (sl_variable.ne.1 .and. sl_variable.ne.2) call panic(3,"Slope limiter variable should be 1 or 2.")
		if (dt_crit.ne.1 .and. dt_crit.ne.2) call panic(3,"Criterion for dt should be 1 or 2.")
		
	end subroutine
	
	
	
	!========= Assigns type of boundary condition ====================
	subroutine bc_select(keyword,bc_type)
		implicit none
		character*5, intent(in) :: keyword
		integer, intent(out) :: bc_type
		
		select case (keyword)
			case ("refle")
				bc_type = 1
			case ("trlnq")
				bc_type = 2
			case ("still")
				bc_type = 3
			case ("hump1")
				bc_type = 4	
			case ("trinh")
				bc_type = 5
			case ("trinq")
				bc_type = 6
			case ("times")
				bc_type = 10	
			case default
				write(*,*)'!!!!!!!! ERROR !!!!!!!!'
				write(*,*)'This option for boundary condition does not exist!'
				write(*,*)'Check settings_num.inp file.'
				write(*,*)'Simulation aborted.'
				stop
		end select
		
	end subroutine
	
	
	
	!========= Allocates memory to variables ====================
	subroutine alloca_vars()
		use variables
		use parameters
		implicit none
		
		allocate(q(m,-1:n+2),q_e(m,n),q_w(m,n),q_e_bar(m,n),q_w_bar(m,n),q_as(m,-1:n+2),q_old(m,-1:n+2))
		allocate(f(m,n-1),fv_e(m,n),fv_w(m,n))
		allocate(x(n),dx(n),zb(n),zb_e(n),zb_w(n),eta(n),eta_bar(n))
		allocate(hr(0:n),hl(0:n),ur(0:n),ul(0:n),hstar(0:n),ustar(0:n),sl(0:n),sr(0:n),sstar(0:n),is_wet(0:n))
		allocate(hv_e(n),hv_w(n))
		
		if (is_sensors) then
			allocate(sens_loc(nsens))
			sens_loc = 0
		end if
			
		! initialise values:
		tsim = 0.0d0
		steps = 0
		percent_ws = 0.0d0
		j_ws = 0
		q = 0.0d0
		q_e = 0.0d0
		q_w = 0.0d0
		q_e_bar = 0.0d0
		q_w_bar = 0.0d0
		f = 0.0d0
		x = 0.d0
		dx = 0.0d0
		zb = 0.0d0
		zb_e = 0.0d0
		zb_w = 0.0d0
		hr = 0.0d0
		hl = 0.0d0
		ur = 0.0d0
		ul = 0.0d0
		hstar = 0.0d0
		ustar = 0.0d0
		sl = 0.0d0
		sr = 0.0d0
		sstar = 0.0d0
		is_wet = 1
		eta = 0.0d0
		hv_e = 0.0d0
		hv_w = 0.0d0
		fv_e = 0.0d0
		fv_w = 0.0d0
		eta_bar = 0.0d0
		q_as = 0.0d0
		q_old = 0.0d0
		
	end subroutine
	
	
	
	!========= Reads file with initial conditions ====================
	!	To standarise, input files are in terms of eta, not h
	subroutine read_inp_file(path,file_name,n,x,bed,variables,eta,hlim)
		use parameters
		use tests
		implicit none
		integer, intent(in) :: n
		real*8, intent(in) :: hlim
		character (len=16), intent(in) :: file_name
		real*8, dimension(m,-1:n+2), intent(inout) :: variables
		real*8, dimension(:), intent(inout) :: x,bed
		real*8, dimension(n), intent(inout) :: eta
		character(len=200), intent(in) :: path
		integer :: i
		character*3 :: rub
		
		open(2,file = trim(path)//file_name,status='old')
		read(2,*)rub
		
			do i=1,n
				read(2,*)x(i),bed(i),eta(i),variables(2,i)
				variables(1,i) = eta(i) - bed(i)
			end do
		
		close(2)
		variables(1,0) = variables(1,1)
		variables(1,n+1) = variables(1,n)
		variables(2,0) = variables(2,1)
		variables(2,n+1) = variables(2,n)
		
		variables(1,-1) = variables(1,1)
		variables(1,n+2) = variables(1,n)
		variables(2,-1) = variables(2,1)
		variables(2,n+2) = variables(2,n)
		
		! check that no negative depth is present:
		call test_neg_depth(variables(1,:),'check read_inp_file in initialise.f95.')
		! test for suspiciously high initial velocities:
		call test_high_vel(bed,variables(1,:),variables(2,:),hlim)
		
		
	end subroutine
	
	
	
	!========= opens files to write results ====================
	subroutine create_files(path,n,q,eta,forma_global,is_sensors,nsens,loc)
		use parameters
		implicit none
		integer, intent(in) :: n
		real*8, dimension(n), intent(in) :: eta
		real*8, dimension(m,-1:n+2), intent(in) :: q
		character*30, intent(in) :: forma_global
		character (len=200), intent(in) :: path
		logical, intent(in) :: is_sensors
		integer, dimension(:), optional, intent(in) :: loc
		integer, optional, intent(in) :: nsens
		integer :: i
		real :: t
		
		open(50,file= trim(path)//'h.dat',status='unknown')
		open(51,file= trim(path)//'q.dat',status='unknown')
		open(52,file= trim(path)//'eta.dat',status='unknown')
		
		! writes intial conditions
		write(50,forma_global)(q(1,i),i=1,n)
		write(51,forma_global)(q(2,i),i=1,n)
		write(52,forma_global)(eta(i),i=1,n)
		
		if (is_sensors) then
			t = 0.0
			open(70,file= trim(path)//'sensors.dat',status='unknown')
			write(70,'(A)')'% t(s)		(h)_i (i=1,#sensors)		(hu)_i (i=1,#sensors)'
			write(70,forma_global)t,(q(1,loc(i)),i=1,nsens),(q(2,loc(i)),i=1,nsens)
		end if
		
	end subroutine
	
	
	
	!========= obtains width of each cell ====================
	subroutine cells_width(n,x,cell_w)
		use tests
		implicit none
		integer, intent(in) :: n
		real*8, dimension(:), intent(in) :: x
		real*8, dimension(:), intent(out) :: cell_w
		integer :: i
		
		do i=2,n-1
			cell_w(i) = (x(i+1) - x(i-1))/2.0d0
			! check:
			if (cell_w(i) < 0.0d0) call panic(3,"Only ascending order in x-axis allowed.")
		end do
		cell_w(1) = x(2)-x(1)
		cell_w(n) = x(n)-x(n-1)
		
	end subroutine
	
	
	
	!========= Reads sensors information if present ====================
	subroutine sensors(path,nsens,sens_loc,tw_sens)
		use tests
		implicit none
		integer, intent(in) :: nsens
		character (len=200), intent(in) :: path
		integer, dimension(:), intent(out) :: sens_loc
		real*8, intent(out) :: tw_sens
		integer :: num_sensors,i
		character*1 :: rub
		
		write(*,*)'Reading local_sensors.inp...'
		
		open(1,file = trim(path)//'local_sensors.inp',status='old')
		read(1,*)rub
		read(1,*)rub
		read(1,*)rub
		read(1,*)rub
		read(1,*)num_sensors
		if (num_sensors .ne. nsens) call panic(3,"Number of sensors in local_sensors.inp and settings_sim.inp do not match")
		read(1,*)rub
		read(1,*)(sens_loc(i),i=1,nsens)
		read(1,*)rub
		read(1,*)tw_sens
		close(1)
		
	end subroutine
	
	
	
	!========= Reads time series ===================
	! has access to global variables
	! returns time series vector [time, variable] and variable read (h or q)
	subroutine read_timeseries(boundary)
		use tests
		use parameters
		use variables
		implicit none
		character*1, intent(in) :: boundary
		integer :: n_points,ind
		
		select case (boundary)
		
		case ("w")
		
		write(*,*)'Reading time series for West boundary...'
		open(1,file = trim(wkpath)//'ts_west.inp',status='old')
		read(1,*)n_points
			if (n_points < 2) call panic(3,'At least 2 points required for time series!')
				
		read(1,'(A1)')ts_var_w
			if (ts_var_w.ne."h".and.ts_var_w.ne."q") call panic(3,'Not valid option in ts_west.inp!')
			
			allocate(ts_w(2,n_points))
			do ind=1,n_points
				read(1,*)ts_w(1,ind),ts_w(2,ind)
			end do
		
		close(1)
		
		! check that time series make sense:
		if (ts_w(1,1) > zero) call panic(3,"time 1 in time series should be >0")
		if (ts_w(1,n_points) < tt) call panic(3,"time series should be > simulation time")
		
		case ("e")
		
		write(*,*)'Reading time series for East boundary...'
		open(1,file = trim(wkpath)//'ts_east.inp',status='old')
		read(1,*)n_points
			if (n_points < 2) call panic(3,'At least 2 points required for time series!')
			
		read(1,'(A1)')ts_var_e
			if (ts_var_e.ne."h".and.ts_var_e.ne."q") call panic(3,'Not valid option in ts_east.inp!')
			
			allocate(ts_e(2,n_points))
			do ind=1,n_points
				read(1,*)ts_e(1,ind),ts_e(2,ind)
			end do
		
		close(1)
		
		! check that time series make sense:
		if (ts_e(1,1) > zero) call panic(3,"time 1 in time series should be >0")
		if (ts_e(1,n_points) < tt) call panic(3,"time series should be > simulation time")
		
		case default
		call panic(2,"initialise.f95, read_timeseries")
		
		end select
		
		
	end subroutine
	
	
	
	!========= Translates bc type to string ===================
	subroutine transl_bc(code,type)
		implicit none
		integer, intent(in) :: code
		character(len=45), intent(out) :: type
		
		! default:
		type = 'type not defined!'
		
		if (code == 1) type = "reflective wall"
		if (code == 2) type = "transmissive (h prescr. q extrap.)"
		if (code == 3) type = "prescr. depth zero discharge"
		if (code == 4) type = "custom. for hump test"
		if (code == 5) type = "conserv. Riemann invariant, prescr. h"
		if (code == 6) type = "conserv. Riemann invariant, prescr. q"
		if (code == 10) type = "time series prescribed"
		
	end subroutine
	
	
	
	!========= Prints log file ===================
	subroutine logfile()
		use variables	
		implicit none
		character :: date*8,time*4
		character(len=30) :: limited_variable
		character(len=45) :: bcw_type,bce_type,dt_crit_s
		
		if (sl_variable == 1) limited_variable = 'water surface level'
		if (sl_variable == 2) limited_variable = 'water depth'
		
		if (which_dt == 1) dt_crit_s = 'for wet domains'
		if (which_dt == 2) dt_crit_s = 'for wet-dry fronts'
		
		call transl_bc(bc_w,bcw_type)
		call transl_bc(bc_e,bce_type)
		
		call date_and_time(date,time)
		
		open(60,file= trim(wkpath)//'simulation.log',status='unknown')
		write(60,*)'************************************************************'
		write(60,*)' '
		write(60,*)'         1D Shallow Water Equations Solver'
		write(60,*)'            (FV HLLC approx. solver)'
		write(60,*)'            	version XXX.XXX'
		write(60,*)' '
		write(60,*)'************************************************************'
		write(60,*)' '
		write(60,*)'Simulation started: YYYYMMDD   hhmm'
		write(60,*)'                    ',date,'   ',time
		write(60,*)' '
		write(60,*)'Working directory:'
		write(60,*)trim(wkpath)
		write(60,*)' '
		write(60,*)'Parameters read from settings_sim.inp and settings_num.inp'
		write(60,*)' '	
		write(60,*)'~~~~~Input data:'
		write(60,'(3X,A40,2X,EN12.3)')'Total time of simulation (s):',tt
		write(60,'(3X,A40,2X,A)')'File with initial conditions:',ic_file
		write(60,'(3X,A40,2X,I5)')'Number of cells:',n
		write(60,'(3X,A40,2X,EN12.3)')'Global results printed every (s):',tw
		write(60,'(3X,A40,2X,EN12.3)')'Mannings n:',manning
		write(60,'(3X,A40,2X,I5)')'Number of local sensors recording:',nsens
		if (is_sensors) write(60,'(3X,A)')'*Sensors input read from local_sensors.inp'
		write(60,*)' '	
		write(60,*)'~~~~~Advanced parameters:'
		write(60,'(3X,A40,2X,EN12.3)')'Slope limiter, beta:',beta
		write(60,'(3X,A40,2X,A)')'Slope limiter applied to:',trim(limited_variable)
		write(60,'(3X,A40,2X,EN12.3)')'Courant number:',courant
		write(60,'(3X,A40,2X,EN12.3)')'Zero water depth threshold (m):',hlim
		write(60,'(3X,A40,2X,A)')'Boundary condition at West:',trim(bcw_type)
		write(60,'(3X,A40,2X,A)')'Boundary condition at East:',trim(bce_type)
		write(60,'(3X,A40,2X,A)')'Criterion for computing dt:',trim(dt_crit_s)
		write(60,*)' '
		write(60,*)'Starting simulation... '
		write(60,*)' '	
		
		!For screen:
		write(*,*)' '
		write(*,'(a,1x,en12.3)')'	Total tme of simulation [s]:',tt
		write(*,*)' '
		
	end subroutine
		
	
	!========= Calls all previous subroutines ====================
	subroutine prep_simulation
		use variables
		use parameters
		implicit none
		
		call welcome_screen(time_start)
		call get_wkdir(wkpath)
		call read_settings_num(beta,courant,hlim,sl_variable,bc_w,bc_e,which_dt,wkpath)
		call read_settings_sim(n,ic_file,format_global,tt,tw,manning,is_sensors,nsens,wkpath)
		call alloca_vars()
			if (bc_w == 10) call read_timeseries('w')
			if (bc_e == 10) call read_timeseries('e')
			if (is_sensors) call sensors(wkpath,nsens,sens_loc,tw_sens)
		call read_inp_file(wkpath,ic_file,n,x,zb,q,eta,hlim)
		call create_files(wkpath,n,q,eta,format_global,is_sensors,nsens,sens_loc)
		call cells_width(n,x,dx)
		call logfile()
		
	end subroutine
		
		
	
end module initialise