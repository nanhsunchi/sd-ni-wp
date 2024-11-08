%% compute wind stress (saildrone) using COARE3.6
% 2024/07/23
clear; close all; clc
cd ~/Documents/projects/sd-windpower/SD_windstress/
addpath('~/Documents/MATLAB/COARE3.6/')
addpath(genpath('~/Documents/MATLAB/mymatlab/'))
addpath(genpath('~/Documents/MATLAB/nanmatlab/'))
year = '2022';
platf_num = '1031';
path = '~/Documents/projects/sd-windpower/plot_SD/';
fn = ['airsea-1min-merge-' year '-SD' platf_num '.nc'];
fn_adcp = ['adcp-1min-merge-' year '-SD' platf_num '.nc'];
% ncdisp([path fn]);
% === load air-sea merged file
vars_load = {'time','latitude','longitude','UWND_MEAN','VWND_MEAN',...
'WIND_MEASUREMENT_HEIGHT_MEAN','WAVE_SIGNIFICANT_HEIGHT',...
'TEMP_AIR_MEAN','RH_MEAN','BARO_PRES_MEAN','PAR_AIR_MEAN',...
'WIND_MEASUREMENT_HEIGHT_MEAN','TEMP_SBE37_MEAN','SAL_SBE37_MEAN',...
'WATER_CURRENT_SPEED_MEAN','WATER_CURRENT_DIRECTION_MEAN'};

for i = 1:numel(vars_load)
    tmp = ncread([path fn], vars_load{i});
    eval([vars_load{i} '= tmp;'])
end
time_airsea = time;
clear time
%% load adcp merged file
vars_adcp_load = {'depth','time','vel_east','vel_north'};
for i = 1:numel(vars_adcp_load)
    tmp = ncread([path fn_adcp], vars_adcp_load{i});
    eval([vars_adcp_load{i} '= tmp;'])
end
time_adcp = time;
clear time

% interpolate adcp near-surface current to match time_airsea
vel_east_intp = nan*ones(numel(depth), numel(time_airsea));
vel_north_intp = vel_east_intp;
for i = 1:numel(depth)
    nanx = isnan(vel_east(i,:));
    vel_east_intp(i,:) = interp1(time_adcp(~nanx),vel_east(i,~nanx), time_airsea);
    vel_north_intp(i,:) = interp1(time_adcp(~nanx),vel_north(i,~nanx), time_airsea);
end

%% compare water current in airsea file vs. adcp current
vel_mag = sqrt(vel_north_intp.^2+vel_east_intp.^2);
vel_dir = wrapTo360( atan2d(vel_north_intp,-vel_east_intp)-90 ); % direction to
U_WATER_CURRENT = WATER_CURRENT_SPEED_MEAN.*sind( WATER_CURRENT_DIRECTION_MEAN );
V_WATER_CURRENT = WATER_CURRENT_SPEED_MEAN.*cosd( WATER_CURRENT_DIRECTION_MEAN );

close; wygiwys; wysiwyg;
subplot(4,1,1)
plot( time_airsea, vel_mag(1,:),'.' ); hold on
plot( time_airsea, WATER_CURRENT_SPEED_MEAN,'.' )
ylabel('m/s')
subplot(4,1,2)
plot( time_airsea, vel_dir(1,:),'.' ); hold on
plot( time_airsea, WATER_CURRENT_DIRECTION_MEAN,'.' );
legend('adcp direction','water current direction')
ylabel('degrees')
subplot(4,1,3)
plot( time_airsea, vel_east_intp(1,:), '.'); hold on
plot( time_airsea, U_WATER_CURRENT, '.')
ylabel('m/s')
subplot(4,1,4)
plot( time_airsea, vel_north_intp(1,:), '.'); hold on
plot( time_airsea, V_WATER_CURRENT, '.')
ylabel('m/s')

%% convert to match the coare input
% U_corrected = sqrt( (UWND_MEAN-vel_east_intp(1,:)').^2 + ...
%     (VWND_MEAN-vel_north_intp(1,:)').^2 );
str_var_app = 'vel_corrected_east'; % !!! change this!!!
U_corrected = UWND_MEAN-vel_east_intp(1,:)';
% str_var_app = 'vel_corrected_east'; % !!! change this!!!
% U_corrected = VWND_MEAN-vel_north_intp(1,:)';
close; plot(time_airsea, U_corrected,'.')

zu = 5.2; % m
zt = 2.3; 
zq = 2.3;
SW_dn = PAR_AIR_MEAN./2.4;
ts_depth = depth(1);

%%
clc
% B=coare36vnWarm_et(Jd,U,zu,Tair,zt,RH,zq,P,Tsea,SW_dn,LW_dn,Lat,Lon,zi,
% Rainrate,ts_depth,Ss,cp,sigH,zrf_u,zrf_t,zrf_q)
% B=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq Cd Ch Ce  L  zeta dT_skinx dq_skinx dz_skin Urf Trf Qrf RHrf UrfN TrfN QrfN  lw_net sw_net Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 hrain Qs Evap T10 T10N Q10 Q10N  RH10 P10 rhoa10 gust wc_frac Edis];
  %   1   2   3   4   5   6    7      8   9  10  11  12  13 14 15  16  17  18       19        20     21  22  23   24   25   26   27     28      29  30  31  32 33   34    35     36   37      38   39  40  41  42   43   44   45    46   47    48   49      50   

% it = 1:numel(time_airsea);
Input = cat(2,time_airsea(it)./86400+1,U_corrected(it),TEMP_AIR_MEAN(it),...
    RH_MEAN(it),BARO_PRES_MEAN(it),TEMP_SBE37_MEAN(it),SW_dn(it),500*ones(numel(it),1),...
    latitude(it),longitude(it),...
    SAL_SBE37_MEAN(it),WAVE_SIGNIFICANT_HEIGHT(it));
Res = coare36vnWarm_et(Input(:,1),Input(:,2),zu,Input(:,3),zt,...
    Input(:,4),zq,Input(:,5),Input(:,6),Input(:,7),Input(:,8),...
    Input(:,9),Input(:,10),600,...
    zeros(numel(it),1),ts_depth,Input(:,11),nan*ones(size(it)),...
    nan*ones(size(it)),10,10,10);

%% plot COARE flux output
iRes_plot = [2,3,4,5,13,21,31];
varnm_plot = {'tau','hsb','hlb','hbb','Cd','Urf','rhoa'};

close
wysiwyg; wygiwys; 
for i = 1:numel(iRes_plot)
    subplot(numel(iRes_plot),1,i)
    plot( datetime(2022,1,1)+days(Input(:,1)), Res(:,iRes_plot(i)),'.',MarkerSize=2 )
    ylabel(varnm_plot{i},FontSize=13)
    grid on
end

saveas(gcf,['select_coare_output_' year '-' platf_num '_' str_var_app '.png'])
%% save to file
SELECT_COARE_OUTPUT = Res(:,[1 iRes_plot]);
SELECT_COARE_OUTPUT_varsnm = ['Jd' varnm_plot];
TIME = time_airsea;

README = {'"SELECT_COARE_OUTPUT(time,variables)" is the results from coare3.6 algorithm.',
    'The variable names are in varnm_plot.'};
save(['tau_SD_coare3.6_' year '-SD' platf_num '_' str_var_app '.mat'],...
    'TIME','SELECT_COARE_OUTPUT_varsnm','SELECT_COARE_OUTPUT')
