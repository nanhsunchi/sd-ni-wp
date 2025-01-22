%% Compute velocity's spectra using RC's multitaper method
% 2025/01/20
clear; close all; clc
cd('~/Documents/projects/sd-ni-wp/check_model/check_RTOFS/')
restoredefaultpath
addpath(genpath('~/Documents/MATLAB/mymatlab/'))
addpath('~/Documents/projects/sd-ni-wp/check_model/check_GFS/')

% data paths/ settings
year = '2023';
platf_num = 'RTOFS';
path_cur_SD = '/Users/chi/Documents/projects/sd-ni-wp/check_model/check_RTOFS/';
fname = 'RTOFS.merge.surface.now.20230905-20230922.nc';
%% === Read SD adcp data
tlim_plot = [datetime(str2double(year),9,5), datetime(str2double(year),9,22)];
tticks = tlim_plot(1):days(1):tlim_plot(2);
% === load variables
dtime = datetime(str2double(year),1,1,0,0,ncread([path_cur_SD fname],'time'));
u = squeeze( ncread([path_cur_SD fname], 'U_VELOCITY') );
v = squeeze( ncread([path_cur_SD fname], 'V_VELOCITY') );
u(u> 1e3) = NaN;
v(v> 1e3) = NaN;
lon = ncread( [path_cur_SD fname],'LONGITUDE' );
lon(lon> 180) = lon(lon> 180)-360;
lat = ncread( [path_cur_SD fname],'LATITUDE' );

%% test plot current velocity
close; clc
lonlim = [-65, -62.5];
latlim = [22.5, 25];
ilon = find( (lon(:,1)>= lonlim(1)) & (lon(:,1)< lonlim(2)) );
ilat = find( (lat(1,:)>= latlim(1)) & (lat(1,:)< latlim(2)) );
disp('lon & lat size:')
disp([size(ilon),size(ilat)])
for j = ilon
    for k = ilat
        subplot(3,1,1)
        plot(datenum(dtime), squeeze(u(j,k,:))*100); hold on
        subplot(3,1,2)
        plot(datenum(dtime), squeeze(v(j,k,:))*100); hold on
    end
end
for i = 1:2
    subplot(3,1,i)
    datetick('x','mm/dd')
    ylabel('(cm/s)')
end
subplot(3,1,3)
plot(lon(:,1)); hold on
plot(lat(1,:));

%% compute spectra for one timeseries within the lon/ lat limits
linstyle = {'-','-'};
linewidth = [1.5,1.5,3,3];
col_shade = {'c','m'};
varnms_plot = {'u','v'};
f_coriolis = [0.557e-4, 0.615e-4];
f_coriolis_cpd = f_coriolis*86400/(2*pi);
f_M2_cpd = 2/((24+50/60)/24);
f_diurnal_cpd = 1/((24+50/60)/24);
hl = [];
speclim = [4e-3,1e1];

close all
for i = 1:numel(varnms_plot)
    % u in m/s
    data = squeeze( eval([varnms_plot{i} '(ilon(1),ilat(1),:)']) );
    % plot(data)
    %%%
    nt = numel(data);
    nbin = 64; % no of frequency bands in request.
    bins = newcrtbins( nt,2,2,nbin );
    dof = 2*( diff(bins')+1 );

    [R,f,C] = multitapersp(demean(data), demean(data), 2/(2*numel(data)),1,'R');
    mf = savg(f,bins); mC = mysavg(C,bins);
    sp1 = mC(:,1); sp2 = mC(:,2); 
    cosp = mC(:,3); qdsp = mC(:,4);
    [mcoh,mpha,msig,muppha,mlwpha] = ...
        cohphaconf(95,cosp-sqrt(-1)*qdsp,sp1,sp2,dof');
    [up,lw] = conospec(dof,0.95); % confidence of interval

    % plot spectrum
    clear h1p h1
    mf_cpd = mf/(60*60)*86400;
    if ~isempty(col_shade{i})
        h1p = plainshaded(mf_cpd, mf_cpd.*sp1'.*up, mf_cpd, mf_cpd.*sp1'.*lw, ...
            col_shade{i}); alpha(h1p,0.1)
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
    speclim = get(gca,'YLim');
end
%%% add frequency labels on top of the plot
plot( f_M2_cpd*[1 1],speclim,'--','Color',[.2,.2,.2] ); 
text(f_M2_cpd, speclim(2)*1.3,'M2','Color',[.2,.2,.2],'Fontsize',14)
plot( f_diurnal_cpd*[1 1],speclim,'--','Color',[.2,.2,.2] );
text(f_diurnal_cpd, speclim(2)*1.3,'Diurnal','Color',[.2,.2,.2],'Fontsize',14)
plot( f_coriolis_cpd(1)*[1 1],speclim,'--','Color',[.2,.2,.2] ); 
text(f_coriolis_cpd(1), speclim(2)*2,'f(22.5N)','Color',[.2,.2,.2],'Fontsize',14,'Rotation',30)
plot( f_coriolis_cpd(2)*[1 1],speclim,'--','Color',[.2,.2,.2] ); 
text(f_coriolis_cpd(2), speclim(2)*2,'f(25N)','Color',[.2,.2,.2],'Fontsize',14,'Rotation',30)

legend(hl,{'u','v'})
% ylim([1e-4,1e1])

% save figure
% saveas(gcf,append('spectra(nbin=',num2str(nbin),')_uv(QC_',num2str(zlim(1)),'m-',num2str(zlim(2)),...
%     'm)-detide_',year,'-SD',platf_num,'_',...
%     string(datetime(tlim_plot(1),'format','yMMd')),'_', ...
%     string(datetime(tlim_plot(2),'format','MMd')),'.png'))

%% compute spectra for each timeseries within the lon/ lat limits
linstyle = '-';
linewidth = 1;
col = 'b';
col_shade = 'c';
varnms_txt = 'u';
varnms_plot = 'u';
f_coriolis = [0.557e-4, 0.615e-4];
f_coriolis_cpd = f_coriolis*86400/(2*pi);
f_M2_cpd = 2/((24+50/60)/24);
f_diurnal_cpd = 1/((24+50/60)/24);
hl = [];
speclim = [4e-3,1e1];


data_pick = squeeze( eval(varnms_plot) );
close all
for i = 1:numel(ilon)
    for j = 1:numel(ilat)
        % u/v in m/s
        data = data_pick(ilon(i),ilat(1),:);
        %%%
        nt = numel(data);
        nbin = 64; % no of frequency bands in request.
        bins = newcrtbins( nt,2,2,nbin );
        dof = 2*( diff(bins')+1 );
    
        [R,f,C] = multitapersp(demean(data), demean(data), 2/(2*numel(data)),1,'R');
        mf = savg(f,bins); mC = mysavg(C,bins);
        sp1 = mC(:,1); sp2 = mC(:,2); 
        
        % plot spectrum
        clear h1p h1
        mf_cpd = mf/(60*60)*86400;
        h1 = plot(mf_cpd,mf_cpd.*sp1','LineStyle',linstyle,'Color',col,'LineWidth',linewidth);
        hold on
    end
end
set(gca,'XScale','log','YScale','log')
xlabel('cpd'); ylabel('\Phi_u (m^2 s^{-2})')        
speclim = get(gca,'YLim');
grid on
%%% add frequency labels on top of the plot
plot( f_M2_cpd*[1 1],speclim,'--','Color',[.2,.2,.2] ); 
text(f_M2_cpd, speclim(2)*1.3,'M2','Color',[.2,.2,.2],'Fontsize',14)
plot( f_diurnal_cpd*[1 1],speclim,'--','Color',[.2,.2,.2] );
text(f_diurnal_cpd, speclim(2)*1.3,'Diurnal','Color',[.2,.2,.2],'Fontsize',14)
plot( f_coriolis_cpd(1)*[1 1],speclim,'--','Color',[.2,.2,.2] ); 
text(f_coriolis_cpd(1), speclim(2)*2,'f(22.5N)','Color',[.2,.2,.2],'Fontsize',14,'Rotation',30)
plot( f_coriolis_cpd(2)*[1 1],speclim,'--','Color',[.2,.2,.2] ); 
text(f_coriolis_cpd(2), speclim(2)*2,'f(25N)','Color',[.2,.2,.2],'Fontsize',14,'Rotation',30)
% add frequnecy lines 0.9xf(lowest latitude) and 1.1xf(highest latitude)
plot( min(f_coriolis_cpd)*[0.9 0.9],speclim,'--','Color',[.5 0 0] ); 
plot( max(f_coriolis_cpd)*[1.1 1.1],speclim,'--','Color',[.5 0 0] ); 
%%% save figure
saveas(gcf,append('spectra(nbin=',num2str(nbin),')_',varnms_txt,'(',fname,')_',...
    num2str(latlim(1)),'-',num2str(latlim(2)),'N_',num2str(lonlim(1)),'_',num2str(lonlim(2)),'.png'))

%%
