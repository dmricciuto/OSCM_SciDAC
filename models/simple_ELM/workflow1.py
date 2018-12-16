#!/usr/bin/env python

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os
import numpy

import utils
from common import oscm_dir,dim


sys.path.append(os.environ['WW']+'/run/xuq')
import workflows as wf
import myutils as mu

###################################################################
###################################################################
###################################################################

def get_ydata(xdata_full):

    qoi='gpp'
    twelve = 12



    obsdata = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
    sites_info = utils.read_sites(obsdata)
    nsites = len(sites_info)

    ndata_all = xdata_full.shape[0]

    sim_locs=np.unique(xdata_full[:,0:2],axis=0)
    obs_locs = np.array([sites_info[j][1:3] for j in range(nsites)])
    site_xind =  utils.pick_sites2(sim_locs,obs_locs)
    sitenames = [sites_info[j][0] for j in range(nsites)]

    ydata_mean = np.empty((0,))
    ydata_std = np.empty((0,))
    ydata_all=np.empty((0,0))
    ind_data_all=[]

    stat_ind = 0 # mean
    sitenames_select = []
    for site_id, lon,lat in site_xind:

        if qoi=='gpp':
            qoi_ind=0
            var_daily = obsdata.variables['GPP'][site_id,:]*24*3600*1000

            ydata = utils.daily_to_monthly(var_daily).reshape(-1,twelve).T

            var_monthly=np.nanmean(ydata,axis=1)
            var_monthly_std=np.nanstd(ydata,axis=1)
            ydata_mean   = np.append(ydata_mean, var_monthly)
            ydata_std   = np.append(ydata_std, var_monthly_std)
            if ydata_all.shape[1]==0:
                ydata_all = ydata.copy()
            else:
                ydata_all = np.append(ydata_all, ydata,axis=0)

            ind_data=utils.pick_ind(xdata_full,[lon,lat,qoi_ind,None,stat_ind])
            ind_data_all.extend(ind_data)
            sitenames_select.append(sitenames[int(site_id)])

    print ind_data_all
    if stat_ind==0:
        ydata=ydata_mean.copy()

    # ind_calib=np.zeros((ndata_all,))
    # ind_calib[ind_data_all]=1
    # np.savetxt('ind_calib.dat', ind_calib, fmt='%d')

    # xdata_calib = xdata_full[np.array(ind_data,dtype=int),:]
    # np.savetxt('xdata_calib.txt',  xdata_calib,  fmt='%.2f')



    np.savetxt('ydata_all_mean.txt', ydata_mean, fmt='%.12f')
    np.savetxt('ydata_all_std.txt', ydata_std, fmt='%.12f')
    np.savetxt('ydata_all_sam.txt', ydata_all, fmt='%.12f')


    return sitenames_select,site_xind, ydata, ind_data_all, ydata_all, ydata_mean, ydata_std



###################################################################
###################################################################
###################################################################


ncfile=sys.argv[1]
print("Reading "+ncfile)

simdata = Dataset(ncfile)


#pnames, prange, ptrain, qtrain = utils.read_simdata_input(simdata)
#xdata_all, outnames, ytrain = utils.read_simdata_ytrain(simdata)
#mu.savepk((pnames, prange,  ptrain, qtrain, ytrain, xdata_all, outnames), 'sim')
pnames, prange,  ptrain, qtrain, ytrain, xdata_all, outnames = mu.loadpk('sim')

#print ytrain.shape
#np.savetxt('ytrain_all.dat',ytrain)
print "Read %d ensemble members of %d outputs" % (ytrain.shape[0],ytrain.shape[1])

np.savetxt('xdata_all.txt', xdata_all,fmt='%.2f')

###################################################################

sitenames_select, site_xind, ydata, ind_data_all, ydata_all, ydata_mean, ydata_std = get_ydata(xdata_all)
np.savetxt('ydata.txt', ydata_mean)
np.savetxt('sitenames_select.dat',sitenames_select,fmt='%s')
np.savetxt('sitecoords_select.dat',site_xind[:,1:],fmt='%.2f')
###################################################################


nout = xdata_all.shape[0]
ind_surr = np.array(ind_data_all) #np.arange(nout) #np.arange(2,93) #for the sake of checking #np.loadtxt('ind_surr.dat',dtype=int)

#print("Computing surrogates")
#results_all=wf.multifit(qtrain, ytrain, xdata_all, pnames, outnames, ind_surr, out_threshold=0.02, order=2, method='bcs')
#mu.savepk(results_all,'results_all')

print("Reading surrogates")
results_all=mu.loadpk('results_all')

assert(ind_surr.shape[0]==len(results_all))
xdata_surr=xdata_all[ind_surr,:]
ytrain_surr=ytrain[:,ind_surr]

np.savetxt('xdata.txt', xdata_surr,fmt='%.2f')
np.savetxt('ytrain.txt', ytrain_surr)


###################################################################

print("Saving surrogate info")

nout_surr = ind_surr.shape[0]
surr_err=np.empty((nout_surr,))
allsens=np.empty((nout_surr,dim))
for iout in range(nout_surr):
    np.savetxt('mindexp.'+str(iout)+'_pred.dat', results_all[iout]['mindex'],fmt='%d')
    np.savetxt('pccfp.'+str(iout)+'_pred.dat', results_all[iout]['cfs'])
    surr_err[iout]=results_all[iout]['rmse_val'] #errval[i]*np.linalg.norm(ydata[:,i])/np.sqrt(ndata)
    #mainsens,totsens,jointsens=mu.pce_sens('LU',results_all[iout]['mindex'],results_all[iout]['cfs'])

    allsens[iout,:]=results_all[iout]['sens'][1]


np.savetxt('surr_errors.dat',surr_err)


###################################################################

nout_calib = nout_surr#len(ind_data_all) #.shape[0]
dataerr_calib=np.empty((nout_calib,))
for iout in range(nout_calib):
    np.savetxt('mindexp.'+str(iout)+'.dat', results_all[iout]['mindex'],fmt='%d')
    np.savetxt('pccfp.'+str(iout)+'.dat', results_all[iout]['cfs'])
    dataerr_calib[iout]=results_all[iout]['rmse_val'] #errval[i]*np.linalg.norm(ydata[:,i])/np.sqrt(ndata)

np.savetxt('dataerr_calib.dat',dataerr_calib)

###################################################################
ind_calib=range(nout_calib)
ytrain_calib=ytrain_surr[:,ind_calib]
xdata_calib=xdata_surr[ind_calib,:]

##################################################################


allsens_ave=np.mean(allsens,axis=0)
np.savetxt('allsens_ave.dat',allsens_ave)

dim_infer=4
ind=np.argsort(allsens_ave)[:-dim_infer]
fixindnom=np.zeros((dim-dim_infer,2))
fixindnom[:,0]=ind
#wf.infer(xdata_surr, ind_calib, ydata, dataerr_calib,dim,\
#        nmcmc=1000,fixindnom=fixindnom,gamma=0.5)

de_params=np.argsort(allsens_ave)[-2:]
wf.infer(xdata_surr, ind_calib, ydata, dataerr_calib,dim,\
        nmcmc=1000,fixindnom=fixindnom,gamma=0.5, order=2, de_params=de_params)

#model_inf -l classical -o 0 -i uniform -t xdata_all.txt -f pcs -c LU -s pci  -u 5 -a -1 -b 1 -m 1000 -g 0.1 -d 3 -z -x xdata_calib.txt -y ydata_calib.txt  -v fixindnom.dat





