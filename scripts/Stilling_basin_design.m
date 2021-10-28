% Generates intial conditions for solving the Stilling Basin problem with
% variable stilling basin geometries (user input). It also runs the HLL
% solver if the option is activated.
%
% Objective: to illustrate the pros and cons of solving a problem like this
% with a numerical model solving the unsteady 1D SWEs, in contrast to
% approaching the problem via simple considerations based on Specific
% Energy theory.
% 
% Sergio Maldonado

clear variables
clc
close all

%==========================================================================
% write files to:
dire = '/Users/sergio/Documents/My_codes/SWE_1D/SWE_1D_HLLC/v2/runs';
%
% Run automatically the SWE solver?
runsolver = true;
path2solver = '/Users/sergio/Documents/My_codes/SWE_1D/SWE_1D_HLLC/v2/src/';
solver = 'SWE1D_HLLC'; % add .exe for Windows

% Input for Stilling Basin geometry:
l = 200.0;          % length of domain (m)
n = 200;            % number of cells
dx = l/n;
x = dx/2.0 : dx :l-dx/2.0;

h = 1.4;            % fixed downstream water depth (m)
u = 1.6;            % fixed downstream velocity (m/s)
q = u*h;

slope = 1/5;        % slope of chute spillway
xch = 50.0;         % chute starts at xch
xsb = 150.00;        % stilling basin ends at xsb

hweir = 0.0;        % height of upstream weir m)
lweir = 3.0;        % length of weir (m)

dz = 2.0;           % difference in elevation between channels (m)
s = 0.1;           % height ot step (m)

% Basic numeric input:
tout = 100.0;       % total time of simulation (s)
twrite = 1.0;       % print results every twrite seconds
manning = 0.0;      % Manning's n

% activating local sensors:
fg_sensor = false;
sens_loc = [5 190]; % position of sensors (in vector index, not m)
tw_local = 1.0;     % write local results every tw_local s
n_sens = fg_sensor*length(sens_loc);

% Advanced numeric parameters/options:
beta_sl = 1.0;      % slope limiter
cour = 0.8;         % Courant (CFL) number
zeroh = 1.0e-3;     % water depth threshold (m)
sl_var = 1;         % slope limiter applied to (1)eta or (2)h
bc_w = 'trinq';     % type of BC at West (see notes below)
bc_e = 'trinq';     % type of BC at East (see notes below)
dt_c = 1;           % criterion for min dt (see notes below)

%==========================================================================

% ----NOTES----------------------------------------------
% Current options for boundary conditions are (must be 5 char long):
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


% date and time:
tnow = datetime;
tnows = datestr(tnow);
thisf = mfilename;              %extracts name of this file with no extension
thisfile = strcat(thisf,'.m');  %adds extension

file_name = 'Stillbas_test';
filename = strcat(file_name,'.ic');

% construct the bathymetry:
m = -slope;
b = (s+dz) - m*xch;
lch = (s+dz)/tan(slope);
xendch = xch + lch;

% datum taken at z=0
zb = (s+dz)*(x<=xch) + s*(x>=xsb) + (m*x+b).*(x>xch & x<xendch) + hweir*(x<=xch & x>=(xch-lweir));
eta = zb+h;

%cplotting initilal conditions
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

% Printing to file
Wr(1,:) = x;
Wr(2,:) = zb;
Wr(3,:) = eta;
Wr(4,:) = q;

header = '%% x(m)   z_b(m)  eta(m)  q(m^2/s)\n';
format = '%6.4e %6.4e %6.4e %6.4e\n';
foot = ['%% File created by ' thisfile ' on ' tnows];

% Calls function to print variables with initial conditions:
%(make sure printfile.m is in the same folder)
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
    