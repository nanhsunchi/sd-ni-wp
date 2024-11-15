%% read saildrone wind data (Hurricane 2022 - try 1031)
clear; close all; clc
year = '2022';
platf_num = '1031';
path = ['~/Documents/Data/Hurricane' year '/delayed/'];
outpath = ['~/Documents/projects/sd-windpower/plot_SD/fig_' platf_num];
fn = ['wind-1min-merge-' year '-SD' platf_num '.nc'];
ncdisp([path fn]);

varnms = {'time','longitude','latitude','UWND_MEAN','VWND_MEAN',...
    'GUST_WND_MEAN','WIND_MEASUREMENT_HEIGHT_MEAN'};
% time is seconds since 1970/1/1 0Z
% load selected variables from file
for i = 1:numel(varnms)
    str_eval = strcat(varnms{i},'= ncread([path fn],',"'",varnms{i},"'",');');
    % disp(str_eval)
    eval(str_eval)
end

%% plot wind
close
dtime = datetime(str2double(year),1,1) + seconds(time);
subplot(2,1,1)
plot(dtime, UWND_MEAN, dtime, VWND_MEAN, ...
    dtime, sqrt(UWND_MEAN.*UWND_MEAN+VWND_MEAN.*VWND_MEAN), ...
    dtime, GUST_WND_MEAN)
legend('u','v','wdsp','gust wdsp','Orientation','horizontal')
grid on
ylabel('(m s^{-1})')
subplot(2,1,2)
plot(dtime, WIND_MEASUREMENT_HEIGHT_MEAN)
grid on
ylabel('(m)')

filename = fullfile(outpath,['wind-' year '-SD' platf_num '.png']);
saveas(gcf,filename)
