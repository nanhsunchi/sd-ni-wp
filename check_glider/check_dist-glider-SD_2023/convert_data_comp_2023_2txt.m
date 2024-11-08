% Check out Lev's glider data_comp_2023.mat and convert it to text file
% 2024/10/25
clear; clc; close all;
fn_mat = 'data_comp_2023.mat';
dataIN = load(fn_mat);
dataIN = dataIN.data_comp;
fn_csv_out = 'SD-glider_2023.csv';

col_nm = {'distance between SD and glider (km) ',...
    'Time of glider profile (datenum)',...
    'Glider Name',...
    'Glider Lat',...
    'Glider Lon',...
    'SD Name',...
    'SD Wind Speed (kts)'};
%%
clc
dist_SD_glider = str2double( dataIN(:,1) );
% plot(dist_SD_glider)
dtime_glider_profile = datetime( str2double(dataIN(:,2)),'ConvertFrom','datenum' );
% disp(dtime_glider_profile)
% === table
T = table(dist_SD_glider, dtime_glider_profile, dataIN(:,3),...
    str2double( dataIN(:,4) ), str2double( dataIN(:,5) ),...
    dataIN(:,6), str2double( dataIN(:,7) ),'VariableNames',col_nm );
% === write table to csv file
writetable(T, fn_csv_out, 'WriteVariableNames', true)
