%Program for results visualisation for the SWE1D FV code

clear variables
clc

% Options ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% files located in:
dire = '/Your/local/path/SWE_1D_HLL/runs/';

plot_anim = 1;          %plot free surface evolution animation? (1/0 = Y/N)
plot_v = 0;             %plot evolution of velocity?
plot_q = 0;             %plot evolution of discharge? 
plot_sensors = 0;       %plot sensors if present?
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Default: --------
tpausa = 0.1;
% -----------------

cd(dire)

% Obtain information from settings_sim.inp:
fid = fopen('settings_sim.inp');
for i=1:6       %headers/not necessary info
    tline = fgetl(fid);
end
initial_cond_file = fgetl(fid); %initial conditions file
for i=1:3
    tline = fgetl(fid);
end
st_tw = fgetl(fid);
dt =  str2double(st_tw);   %results interval
for i=1:3
    tline = fgetl(fid);
end
st_ns = fgetl(fid);
nsens = str2double(st_ns);
fclose(fid);

% Loading data files:
disp('Loading files...')
eta = load('eta.dat');
ic = load(initial_cond_file);
load q.dat
load h.dat
disp('loaded!')

[m,n] = size(q);

x = ic(:, 1);
ni = 0:n-1;
zb = ic(:, 2);
dx = x(length(x)) - x(length(x) - 1);
xf = x(length(x)) + dx;
eta_up = max(max(eta)) + 0.2*max(max(eta));

% plot free sruface evolution:
if (plot_anim)

    disp('plotting...')
    tsim = 0;
    figure
    for i = 1:m
    plot(x,zb,'-k','LineWidth',2)
    tsim = tsim + dt;
    tit = strcat('Free surface evolution, t = ',num2str(tsim),' s');
    title(tit)
    hold on
    plot(x,eta(i,:),'+b','LineWidth',2)
    hold off
    set(gca, 'YLim',[0 eta_up],'XLim',[x(1) xf])
    xlabel('x (m)')
    ylabel('z (m)')

     if (i==1)
        pause
        elseif(i>1)
        pause(tpausa) 
        %pause
     end
    end

end

% plotting velocity:
if (plot_v)
    disp('plotting velocities...')
    v = zeros(m,n);
    for j=1:m
        for i=1:n
            if (h(j,i) > 0.0)
                v(j,i) = q(j,i)/h(j,i);
            end
        end
    end
    v_up = max(max(v)) + 0.2*max(max(v));
    v_min = min(min(v)) - 0.2*min(min(v));
    
    figure
    tsim = 0;
    for j = 1:m
    plot(x,v(j,:),'-r','LineWidth',2)
    tsim = tsim + dt;
    tit = strcat('Velocity, t = ',num2str(tsim),' s');
    title(tit)
    hold on
    set(gca, 'YLim',[v_min v_up],'XLim',[x(1) xf])
    xlabel('x (m)')
    ylabel('velocity (m/s)')
    
        if (j==1)
        pause
        elseif(j>1)
        pause(tpausa)    
        end
    end
    plot(x,v(1,:),'--k','LineWidth',2)
    plot(x,v(m,:),'-b','LineWidth',2)
    
end

% plotting discharge:
if (plot_q)
    disp('plotting discharge...')
    
    q_up = max(max(q)) + 0.2*max(max(q));
    q_min = min(min(q)) - 0.2*min(min(q));
    
    figure
    tsim = 0;
    for j = 1:m
    %hold on    
    plot(x,q(j,:),'*m','LineWidth',2)
    tsim = tsim + dt;
    tit = strcat('Discharge, t = ',num2str(tsim),' s');
    title(tit)
    hold off
    set(gca, 'YLim',[q_min q_up],'XLim',[x(1) xf])
    xlabel('x (m)')
    ylabel('discharge (m^2/s)')
    
        if (j==1)
        pause
        elseif(j>1)
        pause(tpausa)    
        end
    end
    
end



% Plotting sensors:
if (plot_sensors)
    if (nsens == 0) 
        error('No sensors recording according to input file')
    else
        disp('plotting sensors time series...')
        s = load('sensors.dat');
        labels = strings(1,nsens);
        
        figure
        subplot(1,2,1)
        for i=2:nsens+1
            plot(s(:,1),s(:,i))
            hold on
        end
        hold off
        title('Water depth')
        xlabel('time (s)')
        ylabel('water depth (m)')
        
        subplot(1,2,2)
        j=1;
        for i=2+nsens:nsens+nsens+1
            plot(s(:,1),s(:,i))
            hold on
            labels{j} = strcat('sensor ',num2str(j));
            j=j+1;
        end
        title('Discharge')
        xlabel('time (s)')
        ylabel('discharge (m^2/s)')
        legend(labels)
        
    end
end