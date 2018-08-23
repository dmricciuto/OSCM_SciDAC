from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os
import numpy
import utils


if os.environ['USER']=='ksargsy' or os.environ['USER']=='root':
  print('Hello Khachik')
  oscm_dir=os.environ['HOME']+'/research/OSCM_SciDAC/'
elif os.environ['USER']=='csafta':
  print('Hello Cosmin')
  oscm_dir=os.environ['HOME']+'/Projects/OSCM_SciDAC.dmr/'
else:
  oscm_dir='../../'

pdict={'gdd_crit': 500.0, 'crit_dayl': 39300., 'ndays_on':30, 'ndays_off': 15,          \
                      'nue': 15.0, 'slatop':0.03,                                                      \
                      'livewdcn': 50, 'leafcn': 25, 'frootcn': 42,                                     \
                      'fstor2tran': 0.5, 'stem_leaf': 2.7, 'croot_stem': 0.3, 'f_livewd':0.1,          \
                      'froot_leaf': 1.0,                                                               \
                      'rg_frac': 0.3, 'br_mr': 2.52e-6, 'q10_mr': 1.5, 'cstor_tau':3.0,                \
                      'r_mort': 0.02, 'lwtop_ann': 0.7, 'leaf_long': 3.0, 'froot_long': 3.0,           \
                      'q10_hr': 1.5, 'k_l1': 1.2039728, 'k_l2':0.0725707, 'k_l3':0.0140989,            \
                      'k_s1':0.0725707, 'k_s2':0.0140989244, 'k_s3':0.00140098, 'k_s4':0.0001,         \
                      'k_frag':0.0010005, 'rf_l1s1':0.39, 'rf_l2s2':0.55, 'rf_l3s3':0.29,              \
                      'rf_s1s2':0.28, 'rf_s2s3':0.46, 'rf_s3s4':0.55, 'soil4ci':1000.,                 \
                      'cwd_flig':0.24, 'fr_flig':0.25, 'lf_flig':0.25,'fr_flab':0.25, 'lf_flab':0.25,  \
                      'fpi':0.1, 'fpg':0.9}


pmin = {}
pmax = {}

for p in pdict:
    if (p == 'crit_dayl'):
        pmin[p] = 36000.
        pmax[p] = 43000.
    elif (p == 'gdd_crit'):
        pmin[p] = 100.0
        pmax[p] = 700.0
    elif (p == 'fpg'):
        pmin[p] = 0.7 #0.50
        pmax[p] = 0.95 #1.00
    else:
        pmin[p] = pdict[p]*0.5
        pmax[p] = pdict[p]*1.5

ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)

dataset = Dataset(ncfile)
qoi='gpp'
twelve=12

lat_id=34
lon_id=23


monthnames=['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
print("Dimensions #######################")
for ikey in dataset.dimensions.keys(): # time(360),lon(5),lat(7)
    print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size)) # 7

print("Variables #######################")
pnames=[]
ptrain_list=[]
prange_list=[]
for ivar in dataset.variables.keys():
    print(ivar+str(dataset.variables[ivar].shape))
    for attr in dataset.variables[ivar].ncattrs():
        print attr , '=', getattr(dataset.variables[ivar], attr)
    if ivar=='lon':
        lons=dataset.variables[ivar][:]-360
    elif ivar=='lat':
        lats=dataset.variables[ivar][:]
    elif ivar=='time':
        times=dataset.variables[ivar][:]
    elif ivar=='pft_frac':
        pfts=dataset.variables[ivar][:]
    elif ivar=='gpp':
        pass
        #gpps=dataset.variables[ivar][:]
    elif ivar=='lai':
        pass
        #lais=dataset.variables[ivar][:]
    elif ivar=='nue_grass' or ivar=='nue_tree' or ivar=='slatop_decid' or ivar=='slatop_everg':
        pass
    else:
        pnames.append(ivar)
        ptrain_list.append(dataset.variables[ivar][:])
        prange_list.append([pmin[ivar],pmax[ivar]])
        pass

ptrain=np.array(ptrain_list).T
prange=np.array(prange_list)
print pnames
print prange.shape
print ptrain.shape
np.savetxt('ptrain_all.dat',ptrain)
np.savetxt('prange.dat',prange)

sys.exit()
lons=dataset.variables['lon'][:]-360. #.shape
lats=dataset.variables['lat'][:] #.shape


nens=dataset.variables[qoi].shape[0]
print nens

fig=plt.figure(figsize=(12,6))
ytrain=np.empty((nens,twelve))
for ens_id in range(nens):
    thisvar=np.average(dataset.variables[qoi][ens_id,0,:,lat_id,lon_id].reshape(30,twelve),axis=0)
    ytrain[ens_id,:]=thisvar
    #plot(thisvar,linewidth=0.1)

np.savetxt('ytrain_all.dat',ytrain)
#show()
sys.exit()
