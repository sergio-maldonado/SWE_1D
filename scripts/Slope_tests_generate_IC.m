% Select test case for SWE1D FV code
% produces initial conditions to test inclusion of slope-related source
% terms

clear variables
clc

%==========================================================================
% write files to (your working directory):
dire = '/Your/local/path/SWE_1D_HLL/runs/';
% Select case to be run:
testcase = 4;
% Options are:
%
%   0 = Dam break over horizontal wet bed
%   1 = A pond with still water
%   2 = Dam break over a slope (wet bed)
%   3 = Dam break over a slope (dry slope)
%   4 = Uniform flow (to test bed friction)
%   5 = Dam break over triangular hump (Brufau et al. 2002)
%   6 = Supercritical steady flow over a hump
%   7 = Subrcritical steady flow over a hump
%   8 = A ramp to test time series (a tide)
%
% Run automatically the SWE solver?
runsolver = 1;
path2solver = '/Your/loca/path//SWE_1D_HLL/src/';
solver = 'SWE1D_HLLC';  % add .exe for Windows
%==========================================================================

% ----NOTES----------------------------------------------
% Options for boundary conditions are (5 char long)
%  refle = reflective walls
%  trlnq = transmissive (h fixed, q linearly extrapolated)
%  still = fixed h and zero q, for still water test 
%  hump1 = h/q fixed downstream/upstream, for hump test
%  trinh = conserves Riemann invariant, fixed intial h
%  trinq = conserves Riemann invariant, fixed intial q
%  times = reads external time series
%
% Criteria for minimum dt:
%  1 = considers all domain wet; S = |u| + a
%  2 = dry front at some point; S = |u| + 2a
% -------------------------------------------------------

%date and time:
tnow = datetime;
tnows = datestr(tnow);
thisf = mfilename;              %extracts name of this file with no extension
thisfile = strcat(thisf,'.m');  %adds extension

file_pref = 'Slope_test_';

switch testcase
    
    % ************************************
    case 0
        l = 50.0;
        eta_l = 1.0;
        eta_r = 0.2;
        x0 = 25.0;     %position of dam
        %number of cells:
        n = 100;
        dx = l/n;
        x = dx/2.0 : dx :l-dx/2.0;
    
        zb_f = zeros(1,n);
        zb = zb_f;
        
        eta = eta_l.*(x<=x0) + eta_r.*(x>x0);

        q = 0.0;
        
        tout = 5.0;
        twrite = 0.2;
        manning = 0.0;
        
        %numerics:
        beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 1;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'refle';     %type of BC at West (see notes on top)
        bc_e = 'refle';     %type of BC at East (see notes on top)
        dt_c = 1;           %criterion for min dt (see notes on top)
        
        %activating local sensors:
        fg_sensor = false;
        sens_loc = [5,20,35];
        tw_local = 20.0;
        
        n_sens = fg_sensor*length(sens_loc);
        
   % ************************************     
    case 1
        l = 1000;
        n = 40;
        dx = l/n;
        x = dx/2.0 : dx :l-dx/2.0;
        x2 = l./2;
        hs = 5.0 - 4.0*((x-x2)./x2).^2;
        eta = 6.0*ones(1,n);
        zb = eta - hs;
        q = 0.0*ones(1,n);
        
        tout = 1000.0;
        twrite = 50.0;
        manning = 0.0;
        
        %numerics:
        beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 1;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'still';     %type of BC at West (see notes on top)
        bc_e = 'still';     %type of BC at East (see notes on top)
        dt_c = 1;           %criterion for min dt (see notes on top)
        
        %activating local sensors:
        fg_sensor = true;
        sens_loc = [20];
        tw_local = 5.0;
        
        n_sens = fg_sensor*length(sens_loc);
    
    % ************************************
    case 2
        l = 50.0;
        x_r = 30.0;     %ramp starts at
        l_r = l - x_r;  %length of ramp
        z_r = 0.1;      %height of ramp
        eta_l = 1.0;
        eta_r = 0.2;
        x0 = 25.0;     %position of dam
        %number of cells:
        n = 100;
        dx = l/n;
        x = dx/2.0 : dx :l-dx/2.0;
    
        zb_f = zeros(1,n);
        zb_s = (z_r/l_r)*x - (z_r/l_r)*x_r;
        zb = zb_f.*(x<=x_r) + zb_s.*(x>x_r);
        
        eta = eta_l.*(x<=x0) + eta_r.*(x>x0);

        q = 0.0;
        
        tout = 5.0;
        twrite = 0.2;
        manning = 0.0;
        
        %numerics:
        beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 1;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'refle';     %type of BC at West (see notes on top)
        bc_e = 'refle';     %type of BC at East (see notes on top)
        dt_c = 1;           %criterion for min dt (see notes on top)
        
        %activating local sensors:
        fg_sensor = false;
        sens_loc = [5,20,35];
        tw_local = 20.0;
        
        n_sens = fg_sensor*length(sens_loc);
        
   % ************************************     
    case 3
        l = 50.0;
        x_r = 30.0;     %ramp starts at
        l_r = l - x_r;  %length of ramp
        z_r = 0.1;      %height of ramp
        eta_l = 1.0;
        %eta_r = 0.0;
        x0 = 25.0;     %position of dam
        %number of cells:
        n = 100;
        dx = l/n;
        x = dx/2.0 : dx :l-dx/2.0;
    
        zb_f = zeros(1,n);
        zb_s = (z_r/l_r)*x - (z_r/l_r)*x_r;
        zb = zb_f.*(x<=x_r) + zb_s.*(x>x_r);
        
        eta_r = zb;
        eta = eta_l.*(x<=x0) + eta_r.*(x>x0);

        q = 0.0;
        
        tout = 5.0;
        twrite = 0.2;
        manning = 0.0;
        
        %numerics:
        beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 1;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'refle';     %type of BC at West (see notes on top)
        bc_e = 'refle';     %type of BC at East (see notes on top)
        dt_c = 2;           %criterion for min dt (see notes on top)
        
        %activating local sensors:
        fg_sensor = true;
        sens_loc = [80 90];
        tw_local = 0.5;
        
        n_sens = fg_sensor*length(sens_loc);
        
   % ************************************     
    case 4
        l = 1000.0;         %length of channel
        slope = 1./1000;
        datum = 1.0;
        %number of cells:
        n = 40;
        dx = l/n;
        x = dx/2.0 : dx :l-dx/2.0;
        
        zb = -slope*x + slope*l + datum;
        
        h = 5.0;
        q = 0.0;
        eta = zb + h;
        
        tout = 10000.0;
        twrite = 500.0;
        manning = 0.04;
        
        %numerics:
        beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 2;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'trinh';     %type of BC at West (see notes on top)
        bc_e = 'trinh';     %type of BC at East (see notes on top)
        dt_c = 1;           %criterion for min dt (see notes on top)
        
        %activating local sensors:
        fg_sensor = true;
        sens_loc = [5,20,35];
        tw_local = 20.0;
        
        n_sens = fg_sensor*length(sens_loc);
    
    % ************************************    
    case 5
        l = 38.0;
        zs = 25.5;      %hump starts at
        zl = 6.0;       %length of hump
        zh = 0.4;       %height of hump
        
        hd = 0.75;      %water depth in reservoir
        xd = 15.5;      %position of dam
        
        %number of cells:
        n = 380;
        dx = l/n;
        x = dx/2.0 : dx :l-dx/2.0;
    
        zb = zeros(1,n);
        zb_uh = (2.*zh/zl)*x - (2.*zh/zl)*zs;
        zb_dh = (-2.*zh/zl)*x + (2.*zh/zl)*(zs + zl);
        
        x_ch = zs + zl/2;   %coordinate of center of hump
        x_eh = zs + zl;     %coordinate of end of hump
        
        for i=1:n
            if (x(i)<=zs)
                zb(i) = 0.0;
            elseif (x(i)>zs && x(i)<=x_ch)
                zb(i) = zb_uh(i);
            elseif (x(i)>x_ch && x(i)<x_eh)
                zb(i) = zb_dh(i);
            else
                zb(i) = 0.0;
            end
        end
                    
        eta_r = zb;
        eta_l = hd;
        
        eta = eta_l.*(x<=xd) + eta_r.*(x>xd);  
        q = 0.0;
        
        tout = 10.0;
        twrite = .2;
        manning = 0.0125;
        
        %numerics:
        beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 1;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'refle';     %type of BC at West (see notes on top)
        bc_e = 'refle';     %type of BC at East (see notes on top)
        dt_c = 2;           %criterion for min dt (see notes on top)
        
        %activating local sensors:
        fg_sensor = false;
        sens_loc = [5,20,35];
        tw_local = 20.0;
        
        n_sens = fg_sensor*length(sens_loc);
        
       % ************************************ 
        case 6
        l = 25.0;
  
        %number of cells:
        n = 256;
        dx = l/n;
        x = dx/2.0 : dx :l-dx/2.0;
    
        zb = zeros(1,n);
        zb_h = 0.2 - 0.05*(x - 10.0).^2;

        for i=1:n
            if (x(i)<=8.0)
                zb(i) = 0.0;
            elseif (x(i)>8.0 && x(i)<12.0)
                zb(i) = zb_h(i);
            else
                zb(i) = 0.0;
            end
        end
                            
        eta = 0.33;  
        q = 0.18;
        
        tout = 500.0;
        twrite = 10.0;
        manning = 0.0;
        
        %numerics:
        beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 1;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'trinq';     %type of BC at West (see notes on top)
        bc_e = 'trinh';     %type of BC at East (see notes on top)
        dt_c = 1;           %criterion for min dt (see notes on top)
        
        %activating local sensors:
        fg_sensor = false;
        sens_loc = [5,20,35];
        tw_local = 20.0;
        
        n_sens = fg_sensor*length(sens_loc);
        
       % ************************************ 
        case 7
        l = 25.0;
  
        %number of cells:
        n = 256;
        dx = l/n;
        x = dx/2.0 : dx :l-dx/2.0;
    
        zb = zeros(1,n);
        zb_h = 0.2 - 0.05*(x - 10.0).^2;

        for i=1:n
            if (x(i)<=8.0)
                zb(i) = 0.0;
            elseif (x(i)>8.0 && x(i)<12.0)
                zb(i) = zb_h(i);
            else
                zb(i) = 0.0;
            end
        end
                            
        eta = 2.0;  
        q = 4.42;
        
        tout = 500.0;
        twrite = 10.0;
        manning = 0.0;
        
        %numerics:
        beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 1;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'trinq';     %type of BC at West (see notes on top)
        bc_e = 'trinh';     %type of BC at East (see notes on top)
        dt_c = 1;           %criterion for min dt (see notes on top)
        
        %activating local sensors:
        fg_sensor = false;
        sens_loc = [5,20,35];
        tw_local = 20.0;
        
        n_sens = fg_sensor*length(sens_loc);
        
         % ************************************     
        case 8
        l = 50.0;
        x_r = 30.0;     %ramp starts at
        l_r = l - x_r;  %length of ramp
        z_r = 1.0;      %height of ramp
        eta_l = 0.5;
        %eta_r = 0.0;
        %x0 = 25.0;     %position of dam
        %number of cells:
        n = 100;
        dx = l/n;
        x = dx/2.0 : dx :l-dx/2.0;
    
        zb_f = zeros(1,n);
        zb_s = (z_r/l_r)*x - (z_r/l_r)*x_r;
        zb = zb_f.*(x<=x_r) + zb_s.*(x>x_r);
        
        eta_r = zb;
        eta = eta_l.*(zb<=eta_l) + eta_r.*(zb>eta_l);

        q = 0.0;
        
        tout = 7200.0;
        twrite = 100;
        manning = 0.02;
        
        %numerics:
        beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 1;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'times';     %type of BC at West (see notes on top)
        bc_e = 'refle';     %type of BC at East (see notes on top)
        dt_c = 2;           %criterion for min dt (see notes on top)
        
        %activating local sensors:
        fg_sensor = true;
        sens_loc = [5 50];
        tw_local = 10.;
        
        n_sens = fg_sensor*length(sens_loc);
        
    otherwise
        error('This test case is not an option.')
end
% ****************************************************

filename = strcat(file_pref,num2str(testcase),'.ic');

%plotting:
figure
yyaxis left
plot(x,zb,'-k')
hold on
plot(x,eta,'x')
ylim([0 1.25*max(eta)])
ylabel('free surface')

yyaxis right
plot(x,q,'*')
ylim([min(-0.1,1.25*min(q)) max(0.1,1.25*max(q))])
ylabel('discharge')


% Writing file
Wr(1,:) = x;
Wr(2,:) = zb;
Wr(3,:) = eta;
Wr(4,:) = q;

header = '%% x(m)   z_b(m)  eta(m)  q(m^2/s)\n';
format = '%6.4e %6.4e %6.4e %6.4e\n';
foot = ['%% File created by ' thisfile ' on ' tnows];

% Calls function to print variables with initial conditions:
printfile(header,foot,filename,format,Wr,dire)

% --- Prints file with simulation settings:
disp('Writing file with simulation settings (settings_sim.inp)...')
h1 = '%% ====================================\n';
sign = ['\n-------------------\n%% File created by ' thisfile ' on ' tnows '\n'];

fid = fopen('settings_sim.inp','w');
fprintf(fid, h1);
fprintf(fid, '%%  Main settings for the simulation\n');
fprintf(fid, h1);
fprintf(fid, '%% Number of cells:\n');
fprintf(fid, '%i\n', n);
fprintf(fid, '%% File with initial conditions:\n');
fprintf(fid, '%s\n', filename);
fprintf(fid, '%% Total time of simulation (s):\n');
fprintf(fid, '%.3g\n', tout);
fprintf(fid, '%% Print global variables every (s):\n');
fprintf(fid, '%.3f\n', twrite);
fprintf(fid, '%% Mannings coefficient n:\n');
fprintf(fid, '%.5f\n', manning);
fprintf(fid, '%% Amount of local recordings (sensors):\n');
fprintf(fid, '%i\n', n_sens);
fprintf(fid, sign);
fclose(fid);


% --- Prints file with numerical settings:
disp('Writing file with numerics settings (settings_num.inp)...')
fid = fopen('settings_num.inp','w');
fprintf(fid, h1);
fprintf(fid, '%%    Numerics-related settings\n');
fprintf(fid, h1);
fprintf(fid, '%% Slope limiter, beta:\n');
fprintf(fid, '%.3f\n', beta_sl);
fprintf(fid, '%% Courant number:\n');
fprintf(fid, '%.3f\n', cour);
fprintf(fid, '%% Zero depth threshold (m):\n');
fprintf(fid, '%.3g\n', zeroh);
fprintf(fid, '%% Slope limiter applied to (1)free surface or (2)water depth?\n');
fprintf(fid, '%i\n', sl_var);
fprintf(fid, '%% Boundary condition at West (5 char keyword):\n');
fprintf(fid, '%s\n', bc_w);
fprintf(fid, '%% Boundary condition at East (5 char keyword):\n');
fprintf(fid, '%s\n', bc_e);
fprintf(fid, '%% Less/More conservative criterion for min dt (1/2):\n');
fprintf(fid, '%i\n', dt_c);
fprintf(fid, sign);
fclose(fid);

beta_sl = 1.0;      %slope limiter
        cour = 0.8;         %Courant number
        zeroh = 1.0e-3;     %water depth threshold
        sl_var = 2;         %slope limiter applied to (1)eta or (2)h
        bc_w = 'tran1';     %type of BC at West (see notes on top)
        bc_e = 'tran1';    


% --- Prints file with sensors settings if present:
if fg_sensor
    disp('Writing file with information on sensors (local_sensors.inp)...')
    
    fid = fopen('local_sensors.inp','w');
    fprintf(fid, h1);
    fprintf(fid, '%% Information on local sensors\n');
    fprintf(fid, h1);
    fprintf(fid, '%% Number of sensors:\n');
    fprintf(fid, '%i\n', n_sens);
    fprintf(fid, '%% Location (index) of sensors:\n');
    fprintf(fid, '%i ', sens_loc);
    fprintf(fid, '\n%% Output results every (s):\n');
    fprintf(fid, '%.3f\n', tw_local);
    fprintf(fid, sign);
    fclose(fid);
    
    yyaxis left
    plot(x(sens_loc),zb(sens_loc),'dk','MarkerSize',8,'MarkerFaceColor','k') 
end


disp('Done.')
disp('Output files located in:')
disp(pwd)
disp(' ')
        

% Running the simulation from this code: 
% (comment accordingly below Unix/Mac OS or Windows) 
if (runsolver)
    disp('Running executable...')
    % for Unix/Mac OS:
    command = char(strcat('cd',{' '},path2solver,' && ./',solver,{' '},'&'));
    % for Windows:
    % command = char(strcat('cd',{' '},path2solver,' && ',solver,{' '}));
    system(command);
    disp('Simulation ended, see results in working directory.')
end
    
