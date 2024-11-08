% check out rotary spectra of near-surface current (saildrone)
clear; close all; clc
path_cur = '~/Documents/projects/sd-windpower/SD_currents/';
cd(path_cur)
restoredefaultpath
addpath('~/Documents/MATLAB/COARE3.6/')
% addpath(genpath('~/Documents/MATLAB/mymatlab/'))
addpath(genpath('~/Documents/MATLAB/nanmatlab/'))
% settings
year = '2022';
platf_num = '1031';
path_tau = '~/Documents/projects/sd-windpower/SD_windstress/';
path = '~/Documents/projects/sd-windpower/plot_SD/';
fn = ['adcp-1min-merge-' year '-SD' platf_num '.nc'];
fn_wind = ['airsea-1min-merge-' year '-SD' platf_num '.nc'];
pathout = [path_cur 'fig_' year '-SD' platf_num '/'];
if exist(pathout,'dir') == 0
    mkdir(pathout)
end
% === load adcp merged file
vars_load = {'time','latitude','longitude','depth',...
    'vel_east','vel_north'};

for i = 1:numel(vars_load)
    tmp = ncread([path fn], vars_load{i});
    eval([vars_load{i} '= tmp;'])
end
dt_adcp = datetime(str2double(year),1,1)+seconds(time);
time_adcp = time;
clear tmp time

% === load wind data from airsea merged file
vars_load = {'time','UWND_MEAN','VWND_MEAN'};
for i = 1:numel(vars_load)
    tmp = ncread([path fn_wind], vars_load{i});
    eval([vars_load{i} '= tmp;'])
end
dt_airsea = datetime(str2double(year),1,1)+seconds(time);
time_airsea = time;
clear tmp time

%% interpolate adcp near-surface current in time
% 2022-SD1031 mostly on 10 min intervals
% dt_10min = datetime(2022,7,16):minutes(10):datetime(2022,10,31,23,50,0);
iz = 1; % depth index to be focused here
dt_10min = datetime(2022,7,16):minutes(10):datetime(2022,10,18);
time_10min = nan*ones( size(dt_10min) );
for i = 1:numel(time_10min)
    time_10min(i) = seconds(dt_10min(i)-datetime(2022,1,1)+1);
end
vel_east_intp = nan*ones( numel(depth), numel(dt_10min) );
vel_north_intp = vel_east_intp;
for i = 1:numel(depth)
    nanx = isnan(vel_east(i,:));
    vel_east_intp(i,:) = interp1(dt_adcp(~nanx),vel_east(i,~nanx), dt_10min);
    vel_north_intp(i,:) = interp1(dt_adcp(~nanx),vel_north(i,~nanx), dt_10min);
end
latitude_intp = interp1(dt_adcp,latitude, dt_10min);
longitude_intp = interp1(dt_adcp,longitude, dt_10min);

clear time_adcp vel_east vel_north longitud latitude

% put wind to 10-min averaged
uwnd_10min = nan*ones(size(dt_10min));
vwnd_10min = uwnd_10min;
sec_10min = 60*10;
for i = 1:numel(time_10min)
    it_wnd = find( (time_airsea >= time_10min(i)-0.5*sec_10min) & ...
        (time_airsea<= time_10min(i)+0.5*sec_10min));
    nanx = isnan( UWND_MEAN(it_wnd) );
    if sum(~nanx) > 0
        uwnd_10min(i) = mean( UWND_MEAN(it_wnd) );
        vwnd_10min(i) = mean( VWND_MEAN(it_wnd) );
    end
end

clear time_airsea UWND_MEAN VWND_MEAN
%% plot wind
close
plot(dt_10min, uwnd_10min, dt_10min, vwnd_10min); hold on
plot(dt_10min, sqrt(uwnd_10min.^2+vwnd_10min.^2),'k')
legend('uwnd','vwnd','wdsp')
%% compute wind stress - use formula: tau=Cd*rhoa*|u|*u
u_correct = uwnd_10min - vel_east_intp(iz,:);
v_correct = vwnd_10min - vel_north_intp(iz,:);
[taux,tauy] = ra_windstr_nc( uwnd_10min,vwnd_10min,1.15 );
[taux_c,tauy_c] = ra_windstr_nc( u_correct,v_correct,1.15 );

close; clc
subplot(3,1,1)
plot(dt_10min, uwnd_10min); hold on
plot(dt_10min, vwnd_10min)
legend('uwnd','vwnd')
subplot(3,1,2)
plot( dt_10min, taux, dt_10min, taux_c,'.'); hold on
plot( dt_10min, tauy, dt_10min, tauy_c,'.')
plot( dt_10min, sqrt(taux.^2+tauy.^2), dt_10min, sqrt(taux_c.^2+tauy_c.^2),'.')
legend('tau_x','tau_x(c)','tau_y','tau_y(c)')
subplot(3,1,3)
% plot( dt_10min, (1-taux_c./taux)*100,'x'); hold on
% plot( dt_10min, (1-tauy_c./tauy)*100,'x')
histogram( (1-taux_c./taux)*100,-500:10:500); hold on
histogram( (1-tauy_c./tauy)*100,-500:10:500)
xlim([-100 100]); grid on;
legend('reduction in tau_x','reduction in tau_y','Orientation','horizontal',Location='best')
ylabel('(1-\tau_{correct}/\tau)*100 (%)')
saveas(gcf,[pathout 'taux_tauy_correc(' num2str(depth(iz)) 'm)_' ...
    year '-SD' platf_num ...
    datestr(dt_10min(1),'yyyymmdd-') datestr(dt_10min(end),'yyyymmdd') '.png'])

%% wavelet analysis
Fs = 60*24/10; % cpd
f_inert = coriolisf(mean(latitude_intp))/2/pi*86400; % cpd
f_m2 = 1/((12+25/60)/24); % cpd
close; clc
str_wt_out = {'u','v'};
for v = 1:numel( str_wt_out )
    if str_wt_out{v} == 'u'
        [wt,f,coi] = cwt(vel_east_intp(1,:),Fs);
    elseif str_wt_out{v} == 'v'
        [wt,f,coi] = cwt(vel_north_intp(1,:),Fs);
    else
    end
    % plot wavelet analysis result
    [c,h] = contourf( time_10min/86400+1, f, sqrt(real(wt).^2 + imag(wt).^2),...
        0:0.05:0.25);
    hold on; grid on
    set(h,'Linecolor','none')
    yscale('log'); ylim([min(f),max(f)])
    colorbar()
    cmap = flip(hot);
    colormap(cmap)
    plot(time_10min/86400+1, coi,'--','color',[.8 .8 .8],'LineWidth',2.5)
    ttick = 200:10:300;
    tticklabel = {};
    for i = 1:numel(ttick)
        tticklabel{i} = datestr(datetime(2022,1,1)+days(ttick-1),'mm/dd');
    end
    xticks(ttick); xticklabels(tticklabel);
    ylabel('cpd')
    % add M2 & inertial frequency
    hinert = plot(time_10min([1,end])/86400+1,[f_inert, f_inert],'--','color','blue','LineWidth',2);
    hm2 = plot(time_10min([1,end])/86400+1,[f_m2, f_m2],'--','color','m','LineWidth',2);
    legend([hinert, hm2],'f_{inert}','f_{M2}')

    saveas(gcf,[pathout 'wavelet-' str_wt_out{v} '(' num2str(depth(iz)) 'm)_' ...
        year '-SD' platf_num '_' ...
        datestr(dt_10min(1),'yyyymmdd-') datestr(dt_10min(end),'mmdd') '.png'])
end

%% select a period of time close to FIONA and bandpass inertial period.
% Then compute wind power to inertial wind power input.
% itime = ( dt_10min>= datetime(str2double(year),9,15) ) & ...
%     ( dt_10min<= datetime(str2double(year),10,5) );
% bandpass: 0.8f~1.2f. f(16N) = 4.01e-5 rad/s.
n_f = {[1.2,0.8],[1.3,0.7],[1.4,0.6]};
ubp_n_f = {};
vbp_n_f = {};
nsubplot = 3;
rgb_ncol = lines(numel(n_f));
% plot several bandpass results
close; clc;
% yd_here = time_10min(itime)/86400+1;
subplot(nsubplot,1,1)
plot(dt_10min, vel_east_intp(iz,:),'k.'); hold on
ylabel(['u_{adcp}(' depth(iz) 'm)'])
subplot(nsubplot,1,2)
plot(dt_10min, vel_north_intp(iz,:),'k.'); hold on
ylabel(['v_{adcp}(' depth(iz) 'm)'])

for i = 1:numel(n_f)
    f_itime = coriolisf(mean(latitude_intp(:))); % rad/s
    disp(['mean inertial frequency=',num2str(f_itime),' rad/s. ',...
        'inertial period=',num2str( 1/(f_itime/2/pi*86400)*24 ),' hr'])
    flow_cutt = n_f{i}(1)*f_itime/2/pi; % Hz
    fhi_cutt =  n_f{i}(2)*f_itime/2/pi;
    nlow_cutt = round( (1/flow_cutt) );
    nhi_cutt = round( (1/fhi_cutt) );
    % clc
    disp('low and high cutt off frequency: (Hz)')
    disp([flow_cutt,fhi_cutt])
    % bandpass
    [ubp_n_f{i},d] = bandpass(vel_east_intp(iz,:),sort([flow_cutt,fhi_cutt]),1./600);
    [vbp_n_f{i},d] = bandpass(vel_north_intp(iz,:),sort([flow_cutt,fhi_cutt]),1./600);
    subplot(nsubplot,1,1)
    plot(dt_10min, ubp_n_f{i},'-',color=rgb_ncol(i,:))
    text(dt_10min(round((i*0.15-0.1)*numel(dt_10min))),0.45,...
        [num2str(n_f{i}(2)) 'f~' num2str(n_f{i}(1)) 'f'],color=rgb_ncol(i,:),fontsize=14)
    subplot(nsubplot,1,2)
    plot(dt_10min, vbp_n_f{i},'-',color=rgb_ncol(i,:))
    % compute wind power
    tau_dot_sfcvel = taux_c .* ubp_n_f{i} + tauy_c .*vbp_n_f{i};
    subplot(nsubplot,1,3)
    plot(dt_10min, tau_dot_sfcvel*1000,'-',color=rgb_ncol(i,:)); hold on
    ylabel('tau \cdot u_{inertial} (mW/m^{-2})')
end

% figure settings
for i = 1:nsubplot
    subplot(nsubplot,1,i)
    xlim([min(dt_10min) max(dt_10min)])
    if i <= 2
        ylim([-0.55 0.55])
    end
    if i == 3
        ylim([-50 50])
    end
    grid on
end

% compute energy input
nanx = isnan( tau_dot_sfcvel );
disp( [num2str(sum(nanx)/numel(nanx)*100) '% nan data in energy input'] )
energy_input = nan*ones( size(tau_dot_sfcvel) );
energy_input(1) = 0;
for i = 2:numel(tau_dot_sfcvel)
    nanx = isnan(tau_dot_sfcvel(1:i));
    energy_input(i) = trapz( time_10min(~nanx),tau_dot_sfcvel(~nanx) );
end

yyaxis right
plot(dt_10min, energy_input*1e-3,'k.')
grid on; ylabel('\int\Pi dt (kJ/m^2)')

saveas(gcf,[pathout 'uv(' num2str(depth(iz)) 'm)-adcp(raw-bpass)_wind-power-input_' ...
    year '-SD' platf_num '_' ...
    datestr(dt_10min(1),'yyyymmdd-') '-' datestr(dt_10min(end),'mmdd') '.png'])