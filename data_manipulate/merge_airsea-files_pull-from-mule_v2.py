import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import datetime
import os
import math
import cftime
import sys
### The main difference from previous version is only go through files once. Faster.
###
year = '2022'
print('This program merges all listed SD airsea for year '+year+'. One SD for a file and plot time series of variables categorized by \n'+\
      '1D (ex: pos, nav, bottom track), 2D (ex: velocity, error vel, percent_good), 3D (ex: echo intensity, correlation, bt_range)')
### list of sd-number for each year
if year == '2021':
    sd_year = ['1031','1040','1045','1048','1060']
elif year == '2022':
    sd_year = ['1031','1032','1040','1059','1078','1083','1084']
elif year == '2023':
    sd_year = ['1031','1036','1040','1041','1042','1045','1057','1064','1065','1068','1069','1083']
elif year == '2024':
    sd_year = ['1057','1068','1069','1083','1091'] #'1030','1031','1036','1040','1041','1042','1045',

###
### folder name that hosts the airsea data
# for importing saildrone data from mule: make sure it is mounted (mule.pmel.noaa.gov)    
for platf_num in sd_year:
    if year == '2021':
        path = '/Volumes/disk3/projects/sdig-external/sdig/saildrone/hurricane_monitoring_2021_updated/daily_files/'+platf_num+'/'
    elif year == '2022':
        path = '/Volumes/disk3/projects/sdig-external/sdig/saildrone/hurricane_monitoring_2022/delayed/'+platf_num+'/'
    elif year == '2023':
        path = '/Volumes/disk3/projects/sdig-external/sdig/saildrone/hurricane_monitoring_2023/delayed_post_mission/sd-'+platf_num+'/'
    elif year == '2024':
        path = '/Volumes/disk3/projects/sdig-external/sdig/saildrone/noaa_hurricane_2024/delayed/'+platf_num+'/'

    ### get the files in the directory 
    filenames_all = np.sort( os.listdir(path) )
    print('Number of Files in for this directory on Mule:',len(filenames_all))
    pathout = os.getcwd()+'/fig_merge_raw_airsea/'
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
    ### the list of variables I want
    vars_yes = ['latitude','longitude','time','UWND_MEAN','UWND_STDDEV','VWND_MEAN','VWND_STDDEV','GUST_WND_MEAN','GUST_WND_STDDEV',\
                'WIND_MEASUREMENT_HEIGHT_MEAN','TEMP_AIR_MEAN','RH_MEAN','BARO_PRES_MEAN','PAR_AIR_MEAN',\
                    'WAVE_DOMINANT_PERIOD','WAVE_SIGNIFICANT_HEIGHT',\
                        'TEMP_DEPTH_HALFMETER_MEAN','TEMP_SBE37_MEAN','TEMP_SBE37_STDDEV','SAL_SBE37_MEAN','SAL_SBE37_STDDEV',\
                            'WATER_CURRENT_SPEED_MEAN','WATER_CURRENT_DIRECTION_MEAN']
    
    ### create a dictionary
    vars_dic = {}
    vars_dic_attr = {}
    attr = ['long_name','units','installed_height']
    vars_no_installed_height = ['latitude','longitude','time']
    ### find the first uncorrupted file
    for f in range( len(filenames) ):
        try: # block raising an exception
            ds = nc.Dataset(path+'/'+filenames[0])
            varnms = list( ds.variables.keys() )
            break
        except: ### do nothing on exception
            pass

    ### add empty items to the dictionary 
    for i in range( len(varnms) ):
        vkey = varnms[i]
        if vkey in vars_yes:
            item = ds.variables[vkey][:]
            vars_dic[vkey] = np.empty( item.shape )
            attr_in_vkey = ds.variables[vkey].ncattrs()
            ### record attributes
            for a in range( len(attr) ):
                if attr[a] in attr_in_vkey: # if attribute of this vkey exist, write out
                    if (attr[a] == 'installed_height') & (vkey in vars_no_installed_height):
                        pass
                    else:
                        str_eval = "ds.variables['" + vkey + "']." + attr[a]
                        vars_dic_attr[vkey+'-'+attr[a]] = eval(str_eval)
            # print(vkey, item.shape, vars_dic[vkey].shape)
    ds.close()
    print(len(varnms),'variables in nc file &', len(vars_dic),'variables are selected to append')

    ### go through each daily nc file & append the selected variables 
    print('Start going through',len(filenames),'nc file & append selected variables')
    filenames_canopen = []
    numOK = 0
    numBAD = 0
    vars_exists = []
    for f in range( len(filenames) ):
        try: # block raising an exception
            ds = nc.Dataset(path+'/'+filenames[f])
            filenames_canopen.append(filenames[f])
            numOK = numOK + 1
            for i in range( len(varnms) ):
                vkey = varnms[i]
                if (vkey in vars_yes): #& (vkey not in vars_nostack):
                    # print(vkey)
                    vars_exists.append(vkey)
                    vkey = varnms[i]
                    item_old = vars_dic[vkey]
                    item_app = np.squeeze( ds.variables[vkey][:] )
                    ### append in time dimension only
                    if numOK == 1:
                        vars_dic[vkey] = item_app
                    else:
                        vars_dic[vkey] = np.concatenate( (item_old, item_app),axis=0 )
            print('Done',filenames[f])
            ds.close()
        except:
            numBAD = numBAD + 1
            print(filenames[f],'cannot be open. skip')
    print(numOK,'can be open.')
    print(numBAD,'cannot be open.')
    print('----- Done going through nc files & append selected variables -----')
    
    print('----- Start plotting: time gaps -----')
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
    d5min = (datetime.datetime(2019,1,1,0,5,0)-datetime.datetime(2019,1,1)).total_seconds()/3600
    d10min = (datetime.datetime(2019,1,1,0,10,0)-datetime.datetime(2019,1,1)).total_seconds()/3600
    d30min = (datetime.datetime(2019,1,1,0,30,0)-datetime.datetime(2019,1,1)).total_seconds()/3600
    d60min = (datetime.datetime(2019,1,1,1,0,0)-datetime.datetime(2019,1,1)).total_seconds()/3600
    i_tgap = np.where( dtest > d30min )[0]
    print(len(i_tgap),'time gaps that are > 30 min')
    print('The corresponding time gaps (in mins)',dtest[i_tgap]*60)
    for i in range( len(i_tgap) ):
        print(dtime[i_tgap[i]].strftime('%m/%d %H:%MZ'),'-',dtime[i_tgap[i]+1].strftime('%m/%d %H:%MZ'))

    plt.suptitle('Time intervals (t(i+1)-t(i)): Hurricane '+year+' - '+platf_num)
    ### save figure
    fig.savefig(pathout+'time_intervals_'+year+'-'+platf_num+'.png', dpi=300,bbox_inches='tight')

    print('----- Start plotting: wind data (1D data) -----')
    ###
    ### plot merged raw wind data (1D data)
    print(list(vars_dic.keys()))
    vars_plot = ['latitude','longitude','UWND_MEAN','UWND_STDDEV','VWND_MEAN','VWND_STDDEV','GUST_WND_MEAN','GUST_WND_STDDEV','WIND_MEASUREMENT_HEIGHT_MEAN']
    nrow = len(vars_plot)
    tlim = [dtime[0],dtime[-1]]
    print('keys in vars_dic:',vars_dic.keys())
    ### initiate plot
    plt.clf()
    plt.rcParams.update({'font.size': 14})
    fig, axes = plt.subplots(nrow,1)
    fig.set_size_inches(20,16)
    for i in range( nrow ):
        vkey = vars_plot[i]
        if vkey in vars_dic.keys():
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
    fig.savefig(pathout+'merge_raw_wind_1D_'+year+'-'+platf_num+'.png', dpi=300,bbox_inches='tight')

    print('----- Start plotting: met data (1D data) -----')
    ###
    ### plot merged raw met data (1D data)
    print(list(vars_dic.keys()))
    vars_plot = ['TEMP_AIR_MEAN','RH_MEAN','BARO_PRES_MEAN','PAR_AIR_MEAN','WAVE_DOMINANT_PERIOD','WAVE_SIGNIFICANT_HEIGHT']#,\
                #  'TEMP_DEPTH_HALFMETER_MEAN','TEMP_SBE37_MEAN','TEMP_SBE37_STDDEV','SAL_SBE37_MEAN','SAL_SBE37_STDDEV',\
                #     'WATER_CURRENT_SPEED_MEAN','WATER_CURRENT_DIRECTION_MEAN']
    nrow = len(vars_plot)
    tlim = [dtime[0],dtime[-1]]
    ### initiate plot
    plt.clf()
    plt.rcParams.update({'font.size': 14})
    fig, axes = plt.subplots(nrow,1)
    fig.set_size_inches(20,16)
    for i in range( nrow ):
        vkey = vars_plot[i]
        if vkey in vars_dic.keys():
            print(vkey)
            data = vars_dic[vkey]
            h = plt.subplot(nrow,1,i+1)
            plt.plot(dtime,data,'.',ms=1)
            plt.xlim(tlim)
            plt.grid()
            plt.legend(h,labels=[vkey],loc='upper left')
    plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.1,hspace=0.05)
    plt.suptitle('Merged Raw Timeseries: Hurricane '+year+' - '+platf_num)
    fig.savefig(pathout+'merge_raw_met_1D_'+year+'-'+platf_num+'.png', dpi=300,bbox_inches='tight')

    print('----- Start plotting: ocean data (1D data) -----')
    ###
    ### plot merged raw ocean data (1D data)
    print(list(vars_dic.keys()))
    vars_plot = ['TEMP_DEPTH_HALFMETER_MEAN','TEMP_SBE37_MEAN','TEMP_SBE37_STDDEV','SAL_SBE37_MEAN','SAL_SBE37_STDDEV',\
                'WATER_CURRENT_SPEED_MEAN','WATER_CURRENT_DIRECTION_MEAN']
    nrow = len(vars_plot)
    tlim = [dtime[0],dtime[-1]]
    ### initiate plot
    plt.clf()
    plt.rcParams.update({'font.size': 14})
    fig, axes = plt.subplots(nrow,1)
    fig.set_size_inches(20,16)
    for i in range( nrow ):
        vkey = vars_plot[i]
        if vkey in vars_dic.keys():
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
    fig.savefig(pathout+'merge_raw_sea_1D_'+year+'-'+platf_num+'.png', dpi=300,bbox_inches='tight')
    plt.close('all')

    print('----- Start write the merged data to one netcdf file -----')
    ###
    ### write the merged data to netcdf file: https://unidata.github.io/python-training/workshop/Bonus/netcdf-writing/
    try: ncfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    ncfname_out = 'airsea-raw-merge-'+year+'-SD'+platf_num+'.nc'
    ncfile = nc.Dataset(ncfname_out,mode='w',format='NETCDF4_CLASSIC') 
    print(ncfile)
    ### creating dimensions
    time_dim = ncfile.createDimension('time', len(dtime)) # unlimited axis (can be appended to).
    for dim in ncfile.dimensions.items():
        print(dim)
    ### creating attributes
    ncfile.title='Merged files for '+year+' SD-'+platf_num+' from '+path
    print(ncfile.title)
    ncfile.subtitle="Only selected variables for COARE flux calculations are here. Temporal resolution is 1-minute."
    print(ncfile.subtitle)
    print(ncfile)
    ### move the writing data to new variables
    ### compute the time as seconds since 2022-01-01
    time_out = np.array([(dtime[i]-datetime.datetime(int(year),1,1)).total_seconds() for i in range(len(dtime))])

    ###
    def isMonotonic(A):
        return (all(A[i] <= A[i + 1] for i in range(len(A) - 1)) or
                all(A[i] >= A[i + 1] for i in range(len(A) - 1)))
    print('is time_out monotonic? -->', isMonotonic(time_out))

    ###
    ### Creating variables: 'UWND_MEAN', 'UWND_STDDEV', 'VWND_MEAN', 'VWND_STDDEV', 'GUST_WND_MEAN', 'GUST_WND_STDDEV', 'WIND_MEASUREMENT_HEIGHT_MEAN'
    time = ncfile.createVariable('time', np.float64, ('time',))
    time.units = 'seconds since '+year+'-01-01'
    time.long_name = 'time'
    ### 1D variables
    ### 1: (time,)
    vkey_attr = vars_dic_attr.keys()
    for vkey in vars_dic.keys():
        if vkey != 'time':
            ### create variables
            str_exec = vkey + "= ncfile.createVariable('" + vkey + "', np.float64, ('" + 'time' + "',))"
            exec(str_exec)
            # print(i,str_exec)
            vkey_attr_4_this_vkey = [x for x in vkey_attr if vkey in x]
            ### add attributes
            for a in vkey_attr_4_this_vkey:
                ind_ = a.find('-')
                attr_this_round = a[ind_+1:]
                # if attr[a] in attr_in_vkey: # if attribute of this vkey exist, write out
                #     if (attr[a] == 'installed_height') & (vkey in vars_no_installed_height):
                #         pass
                #     else:
                str_exec = vkey + "."+ attr_this_round + " = '"+ str( vars_dic_attr[a] ) + "'"
                exec(str_exec)
                print(a,str_exec)
    
    ###
    ### writing data
    # Note: the ":" is necessary in these "write" statements
    time[:] = time_out
    for vkey in vars_dic.keys():
        if vkey != 'time':
            str_exec = vkey + "[:]= vars_dic['" + vkey + "']"
            exec(str_exec)    
    print(ncfile)
    # close the Dataset.
    ncfile.close(); print('Dataset is closed!')    
