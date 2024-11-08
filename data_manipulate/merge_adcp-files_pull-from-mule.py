### Merge ADCP 5-min files: directly from Mule
import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import datetime
import os
import math
import cftime
import sys

###
year = '2024'
print('This program merges all listed SD adcp for year '+year+'. One SD for a file and plot time series of variables categorized by \n'+\
      '1D (ex: pos, nav, bottom track), 2D (ex: velocity, error vel, percent_good), 3D (ex: echo intensity, correlation, bt_range)')
### list of sd-number for each year
if year == '2021':
    sd_year = ['1031','1040','1045','1048','1060']
elif year == '2022':
    sd_year = ['1031','1032','1040','1059','1078','1083','1084']
elif year == '2023':
    sd_year = ['1031','1036','1040','1041','1042','1045','1057','1064','1065','1068','1069','1083']
elif year == '2024':
    sd_year = ['1030','1031','1036','1040','1041','1042','1045','1057','1068','1069','1083','1091']

###
### folder name that hosts the adcp data
# for importing saildrone data from mule: make sure it is mounted (mule.pmel.noaa.gov)    
for platf_num in sd_year:
    if year == '2021':
        path = '/Volumes/disk3/projects/sdig-external/sdig/saildrone/hurricane_monitoring_2021/adcp/'+platf_num+'/'
    elif year == '2022':
        path = '/Volumes/disk3/projects/sdig-external/sdig/saildrone/high_resolution/hurricane_monitoring_2022/2022_atlantic_hurricane_mission/sd-'+platf_num+'/adcp/5min/'
    elif year == '2023':
        path = '/Volumes/disk3/projects/sdig-external/sdig/saildrone/hurricane_monitoring_2023/adcp/'+platf_num+'/'
    elif year == '2024':
        path = '/Volumes/disk3/projects/sdig-external/sdig/saildrone/noaa_hurricane_2024/adcp/'+platf_num+'/'
    print('========== '+year+'-SD'+platf_num+' ==========')
    ### get the files in the directory 
    filenames_all = np.sort( os.listdir(path) )
    print('Number of Files in for this directory on Mule:',len(filenames_all))
    pathout = os.getcwd()+'/fig_merge_raw_adcp/'
    if os.path.isdir(pathout):
        pass
    else:
        os.mkdir(pathout)
    ### retain the filenames wanted only
    filenames = np.array([])
    for file in filenames_all:
        if ('saildrone' in file) & ('.nc' in file):
            filenames = np.append(filenames, file)
    print('Files have "saildrone" & ".nc" in the sub-directory on Mule:',len(filenames))
    print(filenames)
    ### the list of variables I do not want for now
    vars_no = ['trajectory','nav_start_time','nav_end_time','nav_start_latitude','nav_end_latitude',\
            'nav_start_longitude','nav_end_longitude','beam']
    ### the list of variables that do not have time dimension (depth only)
    vars_nostack = ['cell_depth']

    ### try open the nc files and see how many can be opened
    print('check numbers of files can be open or not.')
    numOK = 0
    numBAD = 0
    filenames_canopen = []
    for f in range( len(filenames) ):
        try: # block raising an exception
            ds = nc.Dataset(path+'/'+filenames[f])
            test = ds.variables['cell_depth'][:]
            ds.close()
            numOK = numOK + 1
            filenames_canopen.append(filenames[f])
            # print('Done',filenames[f])
        except: # doing nothing on exception
            numBAD = numBAD + 1
            # print('cannot open file:',filenames[f])
            pass
    print(numOK,'can be open.')
    print(numBAD,'cannot be open.')

    ###
    ### skip if there is no data in the folder
    if len(filenames_canopen) == 0:
        print('No data in the adcp folder of '+year+'-'+platf_num+'. Go to next folder listed.')
    else:
        ### create a dictionary
        vars_dic = {}
        vars_dic_attr = {}
        attr = ['long_name','units']
        ds = nc.Dataset(path+'/'+filenames_canopen[0])
        # print(ds)
        nz = len( ds.variables['cell_depth'][:] )
        varnms = list( ds.variables.keys() )
        ### add empty items to the dictionary 
        for i in range( len(varnms) ):
            vkey = varnms[i]
            if vkey not in vars_no:
                item = ds.variables[vkey][:]
                vars_dic[vkey] = np.empty( item.shape )
                ### record attributes
                for a in range( len(attr) ):
                    ds.variables[vkey].units
                    str_eval = "ds.variables['" + vkey + "']." + attr[a]
                    vars_dic_attr[vkey+'-'+attr[a]] = eval(str_eval)
                # print(vkey, item.shape, vars_dic[vkey].shape)
        print(len(varnms),'variables in nc file &', len(vars_dic),'variables are selected to append')
        ### go through each (mostly daily) nc file & append the selected variables 
        for f in range( len(filenames_canopen) ):
            try: # block raising an exception
                ds = nc.Dataset(path+'/'+filenames_canopen[f])
                for i in range( len(varnms) ):
                    vkey = varnms[i]
                    if (vkey not in vars_no) & (vkey not in vars_nostack):
                        # print(vkey)
                        vkey = varnms[i]
                        item_old = vars_dic[vkey]
                        item_app = np.squeeze( ds.variables[vkey][:] ) # add np.squeeze for 2021
                        ### append in time dimension only
                        if f == 0:
                            vars_dic[vkey] = item_app
                        else:
                            vars_dic[vkey] = np.concatenate( (item_old, item_app),axis=0 )
                        if (f == 0):# | (f == len(filenames)-1):
                            vars_dic[vars_nostack[0]] = ds.variables[vars_nostack[0]][:]
                            # print(ds.variables[vkey]._FillValue)
                            # print(vkey,item_old.shape, item_app.shape, vars_dic[vkey].shape)
                ds.close()
                print('Done',filenames[f])
            except: # doing nothing on exception
                pass
            # cur_spd = np.sqrt( np.square(u)+np.square(v) )*100  # cm/s
        print(vars_dic['time'].shape)

        ### 
        ### compute the time as seconds since year-01-01
        dtime = np.array([datetime.datetime(1970,1,1)+datetime.timedelta(seconds= vars_dic['time'][i]) for i in range(len(vars_dic['time']))])
        time_out = np.array([(dtime[i]-datetime.datetime(int(year),1,1)).total_seconds() for i in range(len(dtime))])
        print(dtime.shape, time_out.shape)
        ### check if time is monotonic and unique
        def isMonotonic(A):
            return (all(A[i] <= A[i + 1] for i in range(len(A) - 1)) or
                    all(A[i] >= A[i + 1] for i in range(len(A) - 1)))
        print('is time_out monotonic? -->', isMonotonic(time_out))

        ###
        ### sort data so that the time is monotonically increasing --> then find unique time
        it_sort = np.argsort(time_out)
        time_sort = time_out[it_sort]
        time_uniq, it_uniq = np.unique(time_sort, return_index=True)
        plt.plot(time_uniq, time_sort[it_uniq],'.')

        vkey_all = vars_dic.keys()
        print(vkey_all)
        for vkey in vkey_all:
            if vkey != 'cell_depth':
                if np.ndim(vars_dic[vkey]) == 1:
                    data = vars_dic[vkey][it_sort]
                    data_uniq = data[it_uniq]            
                elif np.ndim(vars_dic[vkey]) == 2:
                    data = vars_dic[vkey][it_sort,:]
                    data_uniq = data[it_uniq,:]
                elif np.ndim(vars_dic[vkey]) == 3:
                    data = vars_dic[vkey][it_sort,:,:]
                    data_uniq = data[it_uniq,:,:]
                ## replace
                vars_dic[vkey] = data_uniq
                print(vkey,vars_dic[vkey].shape)

        ###
        ### write merged adcp data to netcdf file: https://unidata.github.io/python-training/workshop/Bonus/netcdf-writing/
        try: ncfile.close()  # just to be safe, make sure dataset is not already open.
        except: pass
        ncfname_out = 'adcp-raw-merge-'+year+'-SD'+platf_num+'.nc'
        ncfile = nc.Dataset(ncfname_out,mode='w',format='NETCDF4_CLASSIC') 
        print(ncfile)
        ### creating dimensions
        depth_dim = ncfile.createDimension('depth', len(vars_dic['cell_depth'])) # depth axis
        time_dim = ncfile.createDimension('time', len(time_uniq)) # unlimited axis (can be appended to).
        beam_dim = ncfile.createDimension('beam', 4)
        for dim in ncfile.dimensions.items():
            print(dim)
        ### creating attributes
        ncfile.title='Merged files for '+year+' SD-'+platf_num
        print(ncfile.title)
        ncfile.subtitle="Only selected variables for adcp measurements are here. Temporal resolution is ~10-minute."
        print(ncfile.subtitle)
        print(ncfile)
        ### move the writing data to new variables
        depth_out = vars_dic['cell_depth']

        ###
        ### Creating variables
        depth = ncfile.createVariable('depth', np.float64, ('depth',))
        depth.units = 'meter'
        depth.long_name = 'depth'
        time = ncfile.createVariable('time', np.float64, ('time',))
        time.units = 'seconds since '+year+'-01-01'
        time.long_name = 'time'
        ### 1D variables
        vars_names = ['longitude','latitude','pitch','roll','heading','vehicle_vel_east','vehicle_vel_north','vehicle_vel_up','bt_vel_east','bt_vel_north','bt_vel_up',\
                'vel_east','vel_north','vel_up','error_vel','percent_good_4_beam','percent_good_3_beam','percent_good',\
                    'bt_range','bt_cor','bt_amp','bt_percent_good',\
                        'echo_intensity','correlation']
        ### 1: (time,) 2: (time, cell_depth) 3: (time, beam) 4: (time, beam, cell_depth)
        dim_cat = [1,1,1,1,1,1,1,1,1,1,1,\
            2,2,2,2,2,2,2,\
                3,3,3,3,\
                    4,4]
        for i in range( len(vars_names) ):
            ### create variables with respective dimensions categorized by dim_cat
            if dim_cat[i] == 1:
                str_exec = vars_names[i] + "= ncfile.createVariable('" + vars_names[i] + "', np.float64, ('" + 'time' + "',))"
                exec(str_exec)
                print(i,str_exec)
            elif dim_cat[i] == 2:
                str_exec = vars_names[i] + "= ncfile.createVariable('" + vars_names[i] + "', np.float64, ('" + 'time' + "','"+ 'depth' +"'))"
                exec(str_exec)
                print(i,str_exec)
            elif dim_cat[i] == 3:
                str_exec = vars_names[i] + "= ncfile.createVariable('" + vars_names[i] + "', np.float64, ('" + 'time' + "','"+ 'beam' +"'))"
                exec(str_exec)
                print(i,str_exec)
            else: 
                str_exec = vars_names[i] + "= ncfile.createVariable('" + vars_names[i] + \
                    "', np.float64, ('" + 'time' + "','"+ 'beam' + "','" + 'depth' +"'))"
                exec(str_exec)
                print(i,str_exec)
            # ### add attributes
            for a in range( len(attr) ):
                str_exec = vars_names[i] + "."+ attr[a] + " = '"+ vars_dic_attr[vars_names[i]+'-'+attr[a]] + "'"
                exec(str_exec)
            
        ###
        ### writing data
        # Note: the ":" is necessary in these "write" statements
        depth[:] = depth_out
        time[:] = time_uniq
        for i in range( len(vars_names) ):
            if dim_cat[i] == 1:
                str_exec = vars_names[i] + "[:]= vars_dic['" + vars_names[i] + "']"
            elif (dim_cat[i] == 2 ) | (dim_cat[i] == 3):
                str_exec = vars_names[i] + "[:,:]= vars_dic['" + vars_names[i] + "']"
            else:
                str_exec = vars_names[i] + "[:,:,:]= vars_dic['" + vars_names[i] + "']"
            exec(str_exec)    
        print(ncfile)
        # close the Dataset.
        ncfile.close(); print('Dataset is closed!')

        ### 
        ### plot time intervals (time(i+1)-time(i)) and print time gaps
        fig, axes = plt.subplots(1)
        fig.set_size_inches(22,3)
        test = vars_dic['time']
        dtime = np.array([datetime.datetime(1970,1,1)+datetime.timedelta(seconds= vars_dic['time'][i]) for i in range(len(vars_dic['time']))])
        dtest = np.array([(test[i+1]-test[i])/3600 for i in range( len(test)-1 )]) # 
        plt.semilogy(dtime[:-1],dtest)
        ddays = (np.max(dtime)-np.min(dtime)).days
        xticks = [np.min(dtime)+datetime.timedelta(days=5*i) for i in range(ddays)]
        xtickslabel = [xticks[i].strftime('%m/%d') for i in range(ddays)]
        plt.xticks(xticks)
        plt.gca().set_xticklabels(xtickslabel, rotation=90)
        plt.ylabel('time difference (hours)')
        plt.xlim([dtime[0],dtime[-1]])
        plt.grid()
        print('min time difference=',np.min(dtest)*3600,'sec')
        print('1st & last time in file:',dtime[0].strftime('%m/%d %H:%MZ'),dtime[-1].strftime('%m/%d %H:%MZ'))
        ### find data gaps
        d5min = (datetime.datetime(int(year),1,1,0,5,0)-datetime.datetime(int(year),1,1)).total_seconds()/3600
        d10min = (datetime.datetime(int(year),1,1,0,10,0)-datetime.datetime(int(year),1,1)).total_seconds()/3600
        d30min = (datetime.datetime(int(year),1,1,0,30,0)-datetime.datetime(int(year),1,1)).total_seconds()/3600
        d60min = (datetime.datetime(int(year),1,1,1,0,0)-datetime.datetime(int(year),1,1)).total_seconds()/3600
        # i_tgap = np.where( dtest > d30min )[0]
        i_tgap = np.where( dtest > d30min )[0]
        print(len(i_tgap),'time gaps that are > 30 min')
        print('The corresponding time gaps (in mins)',dtest[i_tgap]*60)
        for i in range( len(i_tgap) ):
            print(dtime[i_tgap[i]].strftime('%m/%d %H:%MZ'),'-',dtime[i_tgap[i]+1].strftime('%m/%d %H:%MZ'))

        plt.suptitle('Time intervals (t(i+1)-t(i)): Hurricane '+year+' - '+platf_num)
        ### save figure
        fig.savefig(pathout+'time_intervals_'+year+'-'+platf_num+'-adcp.png', dpi=300,bbox_inches='tight')
        
        ###
        ### plot merged raw data (1D data)
        print(list(vars_dic.keys()))
        vars_plot = ['longitude','latitude','pitch','roll','heading','vehicle_vel_east','vehicle_vel_north','vehicle_vel_up','bt_vel_east','bt_vel_north','bt_vel_up']
        nrow = len(vars_plot)
        tlim = [dtime[0],dtime[-1]]
        ### initiate plot
        plt.clf()
        plt.rcParams.update({'font.size': 14})
        fig, axes = plt.subplots(nrow,1)
        fig.set_size_inches(20,16)
        for i in range( nrow ):
            vkey = vars_plot[i]
            print(vkey)
            data = vars_dic[vkey]
            h = plt.subplot(nrow,1,i+1)
            plt.plot(dtime,data,'.',ms=1)
            plt.xlim(tlim)
            plt.grid()
            plt.legend(h,labels=[vkey],loc='upper left')
            if i == 0:
                plt.gca().xaxis.set_ticks_position('top')
            elif i < nrow-1:
                plt.gca().set_xticklabels('')
        plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.1,hspace=0.05)
        plt.suptitle('Merged Raw Timeseries: Hurricane '+year+' - '+platf_num)
        fig.savefig(pathout+'merge_raw_1D_'+year+'-'+platf_num+'-adcp.png', dpi=300,bbox_inches='tight')
        
        ### 
        ### plot merged raw data (2D data)
        print(list(vars_dic.keys()))
        vars_plot = ['vel_east','vel_north','vel_up','error_vel','percent_good_4_beam','percent_good_3_beam','percent_good']
        nrow = len(vars_plot)
        tlim = [dtime[0],dtime[-1]]
        zlim = [0,np.max(vars_dic['cell_depth'])]
        ### initiate plot
        plt.clf()
        plt.rcParams.update({'font.size': 14})
        fig, axes = plt.subplots(nrow,1)
        fig.set_size_inches(20,10.2)
        for i in range( nrow ):
            vkey = vars_plot[i]
            print(vkey)
            data = vars_dic[vkey]
            plt.subplot(nrow,1,i+1)
            if i < 4: # velocities
                # dmin = np.nanmin(data)
                # dmax = np.nanmax(data)
                dlim = [-1,1]#np.max([abs(dmin),abs(dmax)])*np.array([-1,1])
                cmap = 'seismic'
            else: # percentage
                dlim = [0,100]
                cmap = 'jet'
            cs = plt.pcolormesh(dtime,vars_dic['cell_depth'],data.transpose(),vmin=dlim[0],vmax=dlim[1],cmap=cmap)
            plt.xlim(tlim)
            plt.ylim(zlim)
            plt.gca().invert_yaxis()
            plt.yticks(np.arange(10,110,20))
            plt.grid()
            # plt.legend(cs,labels=[vkey],loc='upper left')
            plt.gca().set_title(vkey, x=0.02, y=1, pad=-12, fontsize=12, color='gray', fontweight='bold', horizontalalignment='left')
            ### figure settings
            if i == 0:
                plt.gca().xaxis.set_ticks_position('top')
            elif i < nrow-1:
                plt.gca().set_xticklabels('')

            ### add colorbar - m/s
            if i == 0:
                cbar_ax = fig.add_axes([.91,.5,.005,0.3])
                axf = plt.colorbar(cs,orientation='vertical',cax=cbar_ax,extend='both')
                axf.set_label('(m/s)',fontsize=11)
            elif i == 4:
                cbar_ax = fig.add_axes([.91,.12,.005,0.3])
                axf = plt.colorbar(cs,orientation='vertical',cax=cbar_ax,extend='both')
                axf.set_label('(%)',fontsize=11)
        plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.1,hspace=0.05)
        plt.suptitle(pathout+'Merged Raw Timeseries: Hurricane '+year+' - '+platf_num)
        fig.savefig(pathout+'merge_raw_2D_'+year+'-'+platf_num+'-adcp.png', dpi=300,bbox_inches='tight')

        ###
        ### plot merged raw data (3D data)
        print(list(vars_dic.keys()))
        vars_plot = ['echo_intensity','correlation','bt_range']
        nbeam = 4
        nrow = len(vars_plot)*(nbeam-1)
        tlim = [dtime[0],dtime[-1]]
        zlim = [0,np.max(vars_dic['cell_depth'])]
        ### initiate plot
        plt.clf()
        plt.rcParams.update({'font.size': 14})
        fig, axes = plt.subplots(nrow,1)
        fig.set_size_inches(20,14.5)
        for i in range( len(vars_plot) ):
            vkey = vars_plot[i]
            print(vkey)
            data = vars_dic[vkey]
            cmap = 'jet'
            if i == 0: # echo intensity
                # dmin = np.nanmin(data)
                # dmax = np.nanmax(data)
                # dlim = [dmin,dmax]
                dlim = [30,200]
                cmap = 'nipy_spectral'
            elif i == 1: # correlation
                dlim = [0,100]
                cmap = 'jet'
            elif i == 2: # bt_range
                dmin = np.nanmin(data)
                dmax = np.nanmax(data)
                dlim = [dmin,dmax]
            ### plot 4 beams
            for b in range(nbeam):
                if i == 2: ### bt_range is 1D: plot in one plot
                    plt.subplot(nrow,1,i*4+1)
                    data_plot = data[:,b]
                    h = plt.plot(dtime,data_plot,label=str(b+1))
                    plt.legend()
                    plt.gca().set_title(vkey+'(cm)', x=0.02, y=1, pad=-12, fontsize=12, \
                                        color='gray', fontweight='bold', horizontalalignment='left')
                else:
                    plt.subplot(nrow,1,i*4+b+1)
                    data_plot = data[:,b,:].transpose()
                    cs = plt.pcolormesh(dtime,vars_dic['cell_depth'],data_plot,vmin=dlim[0],vmax=dlim[1],cmap=cmap)
                    plt.ylim(zlim)
                    plt.gca().invert_yaxis()
                    plt.yticks(np.arange(10,110,20))
                    plt.gca().set_title(vkey+'('+str(b+1)+')', x=0.02, y=1, pad=-12, fontsize=12, \
                                        color='gray', fontweight='bold', horizontalalignment='left')
                    if (i == 0) & (b == 0):
                        plt.gca().xaxis.set_ticks_position('top')
                    else:
                        plt.gca().set_xticklabels('')
            ### axes settings    
            plt.xlim(tlim)
            plt.grid()
            ### add colorbar - m/s
            if i == 0: # echo intensity (counts)
                cbar_ax = fig.add_axes([.91,.57,.005,0.3])
                axf = plt.colorbar(cs,orientation='vertical',cax=cbar_ax,extend='both')
                axf.set_label('(counts)',fontsize=11)
            elif i == 1:
                cbar_ax = fig.add_axes([.91,.21,.005,0.3])
                axf = plt.colorbar(cs,orientation='vertical',cax=cbar_ax,extend='both')
                axf.set_label('(%)',fontsize=11)
        plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.1,hspace=0.02)
        plt.suptitle('Merged Raw Timeseries: Hurricane '+year+' - '+platf_num)
        fig.savefig(pathout+'merge_raw_4beams_'+year+'-'+platf_num+'-adcp.png', dpi=300,bbox_inches='tight')
        plt.close('all')