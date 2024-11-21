%% Compute velocity's spectra using RC's multitaper method
% 2024/10/14
clear; close all; clc
cd('~/Documents/projects/sd-ni-wp/test_method_buoy-SD/')
restoredefaultpath
addpath(genpath('~/Documents/MATLAB/mymatlab/'))
addpath('~/Documents/projects/sd-ni-wp/data_manipulate/')

% data paths/ settings
year = '2023';
platf_num = '1042';
path_cur_SD = '/Users/chi/Documents/projects/sd-ni-wp/data_manipulate/data_merge_adcp/';
fname = ['adcp-raw-merge-' year '-SD' platf_num '.nc'];
path_tide_SD = '/Users/chi/Documents/projects/sd-ni-wp/data_manipulate/data_SD_uv_tide/';
fname_tide = ['timeseries_uv-tide_' year '-SD' platf_num '.txt'];

%% === Read SD adcp data
tlim_plot = [datetime(str2double(year),10,19), datetime(str2double(year),11,15,12,0,0)];
tticks = tlim_plot(1):days(1):tlim_plot(2);
% === load variables
dtime = datetime(str2double(year),1,1,0,0,ncread([path_cur_SD fname],'time'));
u = ncread([path_cur_SD fname], 'vel_east');
v = ncread([path_cur_SD fname], 'vel_north');
depth = ncread([path_cur_SD fname],'depth');
data_flag = ncread([path_cur_SD fname],'data_flag');
%% ===  put flagged data to NaN
u(data_flag ~= 0) = NaN;
v(data_flag ~= 0) = NaN;

%% === Linear interpolate the velocity to 15 min grid before compute the spectra
clc
dtime_intp = tlim_plot(1):minutes(15):tlim_plot(2);
u_intp = NaN*ones( length(depth),length(dtime_intp) );
v_intp = u_intp;
for z = 1:length(depth)
    isok = ~isnan(u(z,:));
    if sum(isok) > 2
        u_intp(z,:) = interp1( dtime(isok), squeeze(u(z,isok)), dtime_intp );
        v_intp(z,:) = interp1( dtime(isok), squeeze(v(z,isok)), dtime_intp );
    end
end

%% plot SD velocity
close; clc
u_lim = [-40,40]; cint = 3;
speclim = [0,50];
for i = 1:2
    subplot(2,1,i)
    if i == 1
        [C,h]= contourf(datenum(dtime_intp), depth, u_intp*100,...
            u_lim(1):cint:u_lim(2));
    else
        [C,h]= contourf(datenum(dtime_intp), depth, v_intp*100,...
            u_lim(1):cint:u_lim(2));
    end
    set(h,'LineColor','none')
    colormap(redblue)
    caxis(u_lim);
    datetick('x','mm/dd')
    set(gca, 'YDir','reverse','YLim',speclim,'XLim',datenum(tlim_plot))
    colorbar
end

%% compute spectra
col = ['b','r'];
col_shade = {[204,255,255]/255,[255,240,245]/255};
varnms_plot = {'u_intp','v_intp'};
f_coriolis = 0.672e-4;
f_coriolis_cpd = f_coriolis*86400/(2*pi);
f_M2_cpd = 2/((24+50/60)/24);
close all
for i = 1:2
    data = eval([varnms_plot{i} '(1,:)']);
    %%%
    nt = numel(data);
    nbin = 64; % no of frequency bands in request.
    bins = newcrtbins( nt,2,2,nbin );
    dof = 2*( diff(bins')+1 );
    
    [R,f,C] = multitapersp...
        (demean(data), demean(data), 2/(2*numel(data)),1,'R');
    mf = savg(f,bins); mC = mysavg(C,bins);
    sp1 = mC(:,1); sp2 = mC(:,2); 
    cosp = mC(:,3); qdsp = mC(:,4);
    [mcoh,mpha,msig,muppha,mlwpha] = ...
        cohphaconf(95,cosp-sqrt(-1)*qdsp,sp1,sp2,dof');
    [up,lw] = conospec(dof,0.95);

    % plot spectrum
    mf_cpd = mf/(15*60)*86400;
    h1p = plainshaded(mf_cpd, mf_cpd.*sp1'.*up, mf_cpd, mf_cpd.*sp1'.*lw, ...
        col_shade{i}); alpha(h1p,0.5)
    hold on
    h1 = plot(mf_cpd,mf_cpd.*sp1','LineWidth',2,'Color',col(i)); grid on;
    uistack(h1,'top')
    uistack(h1p,'bottom')
    set(gca,'XScale','log','YScale','log')
    xlabel('cpd'); ylabel('\Phi_u (m^2 s^{-2})')
    speclim = get(gca,'YLim');
end

plot( f_coriolis_cpd*[1 1],speclim,'b--' ); 
text(f_coriolis_cpd, speclim(2)*1.3,'f','Color','b')
plot( f_M2_cpd*[1 1],speclim,'b--' ); 
text(f_M2_cpd, speclim(2)*1.3,'M2','Color','b')

% save figure
saveas(gcf,append('spectra_uv(QC)_',year,'-SD',platf_num,'_',...
    string(datetime(tlim_plot(1),'format','yMMd')),'_', ...
    string(datetime(tlim_plot(2),'format','MMd')),'.png'))

%% Read tidal current - AVISO FES 2014 model
table_tide = readtable([path_tide_SD fname_tide]);
[tnum, ~] = size(table_tide);
dtime_tide = datetime(zeros(tnum,1),0,0);
for i = 1:tnum
    dtime_tide(i) = table2array(table_tide(i,1))+...
        table2array(table_tide(i,2));
end
% linear interpolate in time
u_tide = table2array(table_tide(:,5));
v_tide = table2array(table_tide(:,6));
isok = ~isnan(u_tide);
u_tide_intp = interp1( dtime_tide(isok), squeeze(u_tide(isok)), dtime_intp );
v_tide_intp = interp1( dtime_tide(isok), squeeze(v_tide(isok)), dtime_intp );

%% plot the detide SD velocity
close; clc
u_lim = [-40,40]; cint = 3;
zlim = [0,50];
varnms_plot = {'u','u-u_{tide}','v','v_{tide}'};
for i = 1:4
    subplot(4,1,i)
    if i == 1
        [C,h]= contourf(datenum(dtime_intp), depth, u_intp*100,...
            u_lim(1):cint:u_lim(2));
    elseif i == 2
        [C,h]= contourf(datenum(dtime_intp), depth, u_intp*100-u_tide_intp,...
            u_lim(1):cint:u_lim(2));
    elseif i == 3
        [C,h]= contourf(datenum(dtime_intp), depth, v_intp*100,...
            u_lim(1):cint:u_lim(2));
    elseif i == 4
        [C,h]= contourf(datenum(dtime_intp), depth, v_intp*100-v_tide_intp,...
            u_lim(1):cint:u_lim(2));
    end
    set(h,'LineColor','none')
    colormap(redblue)
    caxis(u_lim);
    xticks(datenum(tticks))
    xticklabels(string(tticks,'M/dd'))
    yticks(zlim(1):10:zlim(2))
    % datetick('x','mm/dd')
    set(gca, 'YDir','reverse','YLim',zlim,'XLim',datenum(tlim_plot))
    colorbar
    title(varnms_plot{i},'HorizontalAlignment','left')
end
% save figure
saveas(gcf,append('time-depth_uv(QC)_detide_',year,'-SD',platf_num,'_',...
    string(datetime(tlim_plot(1),'format','yMMd')),'_', ...
    string(datetime(tlim_plot(2),'format','MMd')),'.png'))

%% compute spectra - after removing barotropic tide
linstyle = {'-','-','--','--'};
linewidth = [1.5,1.5,3,3];
col = {'b','r','c','m'};
% col_shade = {'','',[204,255,255]/255,[255,240,245]/255};
col_shade = {'','','c','m'};
varnms_txt = {'u','v','u-u_{tide}','v-v_{tide}'};
varnms_plot = {'u_intp','v_intp','u_intp','v_intp'};
varnms_tide_plot = {'0','0','u_tide_intp','v_tide_intp'};
f_coriolis = 0.672e-4;
f_coriolis_cpd = f_coriolis*86400/(2*pi);
f_M2_cpd = 2/((24+50/60)/24);
hl = [];
speclim = [4e-3,1e1];

close all
for i = 1:numel(varnms_plot)
    % u-u_tide in m/s
    data = eval([varnms_plot{i} '(1,:)'])-0.01*eval(varnms_tide_plot{i});
    %%%
    nt = numel(data);
    nbin = 64; % no of frequency bands in request.
    bins = newcrtbins( nt,2,2,nbin );
    dof = 2*( diff(bins')+1 );
    
    [R,f,C] = multitapersp...
        (demean(data), demean(data), 2/(2*numel(data)),1,'R');
    mf = savg(f,bins); mC = mysavg(C,bins);
    sp1 = mC(:,1); sp2 = mC(:,2); 
    cosp = mC(:,3); qdsp = mC(:,4);
    [mcoh,mpha,msig,muppha,mlwpha] = ...
        cohphaconf(95,cosp-sqrt(-1)*qdsp,sp1,sp2,dof');
    [up,lw] = conospec(dof,0.95); % confidence of interval

    % plot spectrum
    clear h1p h1
    mf_cpd = mf/(15*60)*86400;
    if ~isempty(col_shade{i})
        h1p = plainshaded(mf_cpd, mf_cpd.*sp1'.*up, mf_cpd, mf_cpd.*sp1'.*lw, ...
            col_shade{i}); alpha(h1p,0.25)
    end
    hold on
    h1 = plot(mf_cpd,mf_cpd.*sp1','LineStyle',linstyle{i},'Color',col{i},'LineWidth',linewidth(i)); 
    hl = [hl; h1];
    grid on;
    uistack(h1,'top')
    if exist('h1p')== 1
        uistack(h1p,'bottom')
    end
    set(gca,'XScale','log','YScale','log')
    xlabel('cpd'); ylabel('\Phi_u (m^2 s^{-2})')
    ylim(speclim)
    % speclim = get(gca,'YLim');
end
plot( f_coriolis_cpd*[1 1],speclim,'--','Color',[.2,.2,.2] ); 
text(f_coriolis_cpd, speclim(2)*1.3,'f','Color',[.2,.2,.2],'Fontsize',14)
plot( f_M2_cpd*[1 1],speclim,'--','Color',[.2,.2,.2] ); 
text(f_M2_cpd, speclim(2)*1.3,'M2','Color',[.2,.2,.2],'Fontsize',14)
legend(hl,{'u','v','u-u_{tide}','v-v_{tide}'})

% save figure
saveas(gcf,append('spectra(nbin=',num2str(nbin),')_uv(QC)-detide_',year,'-SD',platf_num,'_',...
    string(datetime(tlim_plot(1),'format','yMMd')),'_', ...
    string(datetime(tlim_plot(2),'format','MMd')),'.png'))

%%
figure
loglog(mf_cpd, dof,'-o'); grid on