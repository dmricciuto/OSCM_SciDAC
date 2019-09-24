import numpy
from netCDF4 import Dataset
#import matplotlibpyplot as plt
from math import sin, cos, sqrt, atan2, radians
import time, os
import utils

if os.environ['USER']=='ksargsy':
  print('Hello Khachik')
  oscm_dir=os.environ['HOME']+'/research/OSCM_SciDAC/'
elif os.environ['USER']=='csafta':
  print('Hello Cosmin')
  oscm_dir=os.environ['HOME']+'/Projects/OSCM_SciDAC.dmr/'
else:
  oscm_dir='../../'



class MyModel(object):

    def __init__(self):
        self.name = 'sELM'
        self.parms = {'gdd_crit': 500.0, 'crit_dayl': 39300., 'ndays_on':30, 'ndays_off': 15,          \
                      'nue_tree': 15.0, 'nue_grass':9.0, 'slatop_everg':0.01, 'slatop_decid':0.03,     \
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
        #Line 1:      Phenology (4 parameters)
        #Line 2:      Photosyntesis (2 parameters)
        #Line 3:      Nitrogen content (3 parmameters)
        #Lines 4-5:   Allocation (5 parameters)
        #Line 6:      Autotrophic respiration (5 parameters)
        #Line 7:      Plant turnover (4 parameters)
        #Lines 8-12:  Soil decomposition and turnover (16 parameters)
        #Line 13:     Litter chemistry(5 parameters)
        #Line 14:     Nutrient limitation (2 parameters)
        self.site='none'
        self.pdefault=self.parms
        self.pmin = {}
        self.pmax = {}
        self.nparms = 0
        for p in self.parms:
            if (p == 'crit_dayl'):
                self.pmin[p] = 36000.
                self.pmax[p] = 43000.
            elif (p == 'gdd_crit'):
                self.pmin[p] = 100.0
                self.pmax[p] = 700.0
            elif (p == 'fpg'):
                self.pmin[p] = 0.70
                self.pmax[p] = 0.95
            else:
                self.pmin[p] = self.parms[p]*0.50
                self.pmax[p] = self.parms[p]*1.50

            self.nparms = self.nparms+1
        self.issynthetic = False
        self.ne = 1

        #Model outputs
        self.outvars = ['gpp','npp','gr', 'mr','hr','nee','lai','leafc','leafc_stor','frootc','frootc_stor','livestemc', \
                        'deadstemc','livecrootc','deadcrootc','ctcpools','totecosysc','totsomc','totlitc', 'cstor']


    def selm_instance(self, parms, spinup_cycles=0, pft=0):

        #--------------- Initialize ------------------------
        #Flux variables
        gpp = self.output['gpp']
        npp = self.output['npp']
        gr  = self.output['gr']
        mr  = self.output['mr']
        hr  = self.output['hr']
        nee = self.output['nee']
        #State variables
        lai         = self.output['lai']
        leafc       = self.output['leafc']
        leafc_stor  = self.output['leafc_stor']
        frootc      = self.output['frootc']
        frootc_stor = self.output['frootc_stor']
        livestemc   = self.output['livestemc']
        deadstemc   = self.output['deadstemc']
        livecrootc  = self.output['livecrootc']
        deadcrootc  = self.output['deadcrootc']
        totecosysc  = self.output['totecosysc']
        ctcpools    = self.output['ctcpools']
        totsomc     = self.output['totsomc']
        totlitc     = self.output['totlitc']
        cstor       = self.output['cstor']

        #Set initial States
        leafc_stor[0] = 0.0
        leafc[0]      = 0.0
        if (pft == 0 or pft == 2):
            leafc_stor[0] = 10.0
        else:
            leafc[0] = 10.0
        frootc_stor[0] = 0.0
        frootc[0]      = 0.0
        livestemc[0]   = 0.0
        deadstemc[0]   = 0.0
        livecrootc[0]  = 0.0
        deadcrootc[0]  = 0.0
        ctcpools[:,0]  = 0.0
        totecosysc[0]  = 0.0
        totsomc[0]     = 0.0
        totlitc[0]     = 0.0
        cstor[0]       = 0.0
        #Set initial soil carbon for long-lived pool
        ctcpools[6] = parms['soil4ci']

        #Forcings
        tmax = self.forcings['tmax']
        tmin = self.forcings['tmin']
        rad  = self.forcings['rad']
        doy  = self.forcings['doy']
        cair = self.forcings['cair']
        dayl = self.forcings['dayl']
        btran = self.forcings['btran']
        #Coefficents for ACM (GPP submodel)
        if (pft == 2):
          nue = parms['nue_grass']
        else:
          nue = parms['nue_tree']
        a = [nue, 0.0156935, 4.22273, 208.868, 0.0453194, 0.37836, 7.19298, 0.011136, \
             2.1001, 0.789798]

        #Turnover times for CTC model
        k_ctc = [parms['k_l1'],parms['k_l2'],parms['k_l3'],parms['k_s1'], \
                 parms['k_s2'],parms['k_s3'],parms['k_s4'],parms['k_frag']]
        #Respiration fractions for CTC model pools
        rf_ctc = [parms['rf_l1s1'],parms['rf_l2s2'],parms['rf_l3s3'] , \
                  parms['rf_s1s2'],parms['rf_s2s3'],parms['rf_s3s4'], 1.0, 0.0]
        #transfer matrix for CTC model
        tr_ctc = numpy.zeros([8,8],numpy.float)
        tr_ctc[0,3] = 1.0 - parms['rf_l1s1']
        tr_ctc[1,4] = 1.0 - parms['rf_l2s2']
        tr_ctc[2,5] = 1.0 - parms['rf_l3s3']
        tr_ctc[3,4] = 1.0 - parms['rf_s1s2']
        tr_ctc[4,5] = 1.0 - parms['rf_s2s3']
        tr_ctc[5,6] = 1.0 - parms['rf_s3s4']
        tr_ctc[7,1] = parms['cwd_flig']
        tr_ctc[7,2] = 1.0 - parms['cwd_flig']


        #Initialize local variables
        gdd     = 0
        leafon  = 0
        leafoff = 0
        annsum_npp = 0
        annsum_npp_temp = 0

        #Run the model
        for s in range(0,spinup_cycles+1):
          totecosysc_last = totecosysc[0]
          if (s > 0):
            leafc_stor[0]  = leafc_stor[self.nobs-1]
            leafc[0]       = leafc[self.nobs-1]
            frootc_stor[0] = frootc_stor[self.nobs-1]
            frootc[0]      = frootc[self.nobs-1]
            livestemc[0]   = livestemc[self.nobs-1]
            deadstemc[0]   = deadstemc[self.nobs-1]
            livecrootc[0]  = livecrootc[self.nobs-1]
            deadcrootc[0]  = deadcrootc[self.nobs-1]
            ctcpools[:,0]  = ctcpools[:,self.nobs-1]
            totecosysc[0]  = totecosysc[self.nobs-1]
            totsomc[0]     = totsomc[self.nobs-1]
            totlitc[0]     = totlitc[self.nobs-1]
            cstor[0]       = cstor[self.nobs-1]

          for v in range(0,self.nobs):
            # --------------------1.  Phenology -------------------------
            #Calculate leaf on
            leafc_trans = 0.0
            frootc_trans = 0.0
            if (pft == 0 or pft == 2):     #Decidous phenology
              gdd_last = gdd
              dayl_last = dayl[v-1]
              gdd_base = 0.0
              gdd = (doy[v] > 1) * (gdd + max(0.5*(tmax[v]+tmin[v])-gdd_base, 0.0))
              if (gdd >= parms['gdd_crit'] and gdd_last < parms['gdd_crit']):
                  leafon = parms['ndays_on']
                  leafc_trans_tot  = leafc_stor[v]*parms['fstor2tran']
                  frootc_trans_tot = frootc_stor[v]*parms['fstor2tran']
              if (leafon > 0):
                  leafc_trans  = leafc_trans_tot  / parms['ndays_on']
                  frootc_trans = frootc_trans_tot / parms['ndays_on']
                  leafon = leafon - 1
              #Calculate leaf off
              if (dayl_last >= parms['crit_dayl']/3600. and dayl[v] < parms['crit_dayl']/3600.):
                   leafoff = parms['ndays_off']
                   leafc_litter_tot  = leafc[v]
                   frootc_litter_tot = frootc[v]
              if (leafoff > 0):
                   leafc_litter  = min(leafc_litter_tot  / parms['ndays_off'], leafc[v])
                   frootc_litter = min(frootc_litter_tot / parms['ndays_off'], frootc[v])
                   leafoff = leafoff - 1
              else:
                   leafc_litter  = 0.0
                   frootc_litter = 0.0
            else:               #Evergreen phenology
                leafc_litter  = leafc[v]  * 1.0 / (parms['leaf_long']*365. )
                frootc_litter = frootc[v] * 1.0 / (parms['froot_long']*365.)

            if (pft == 1):
              slatop = parms['slatop_everg']
            else:
              slatop = parms['slatop_decid']
            lai[v+1] = leafc[v] * slatop

            #---------------------2. GPP -------------------------------------
            #Calculate GPP flux using the ACM model (Williams et al., 1997)
            if (lai[v] > 1e-3):
                rtot = 1.0
                psid = -2.0
                leafn = 1.0/(parms['leafcn'] * slatop)
                gs = abs(psid)**a[9]/((a[5]*rtot+(tmax[v]-tmin[v])))
                pp = max(lai[v],0.5)*leafn/gs*a[0]*numpy.exp(a[7]*tmax[v])
                qq = a[2]-a[3]
                #internal co2 concentration
                ci = 0.5*(cair[v]+qq-pp+((cair[v]+qq-pp)**2-4.*(cair[v]*qq-pp*a[2]))**0.5)
                e0 = a[6]*max(lai[v],0.5)**2/(max(lai[v],0.5)**2+a[8])
                cps   = e0*rad[v]*gs*(cair[v]-ci)/(e0*rad[v]+gs*(cair[v]-ci))
                gpp[v+1] = cps*(a[1]*dayl[v]+a[4])
                #ACM is not valid for LAI < 0.5, so reduce GPP linearly for low LAI
                if (lai[v] < 0.5):
                  gpp[v+1] = gpp[v+1]*lai[v]/0.5
            else:
                gpp[v+1] = 0.0
            gpp[v+1] = gpp[v+1]*btran[v]
            if ((tmax[v]+tmin[v])/2 < 0.0):
              gpp[v+1] = 0.0   #Zero GPP if average temperature below freezing
            #--------------------3.  Maintenace respiration ------------------------
            #Maintenance respiration
            trate = parms['q10_mr']**((0.5*(tmax[v]+tmin[v])-25.0)/25.0)
            mr[v+1] = (leafc[v]/parms['leafcn'] + frootc[v]/parms['frootcn'] + \
                       (livecrootc[v]+livestemc[v])/parms['livewdcn'])* \
                       (parms['br_mr']*24*3600)*trate
            #Nutrient limitation
            availc      = max(gpp[v+1]-mr[v+1],0.0)
            availc      = availc * parms['fpg']
            cstor_alloc = availc * (1.0 - parms['fpg'])
            xsmr        = max(mr[v+1]-gpp[v+1],0.0)

            #---------------4.  Allocation and growth respiration -------------------
            frg  = parms['rg_frac']
            flw  = parms['f_livewd']
            f1   = parms['froot_leaf']
            if (pft < 2):
              f2   = max(parms['stem_leaf']/(1.0+numpy.exp(-0.004*(annsum_npp - \
                             300.0))) - 0.4, 0.1)
              f3   = parms['croot_stem']
            else:
              f2 = 0
              f3 = 0
            fall = (1.0+frg)*(1.0 + f1 + f2*(1+f3))
            if (pft == 0 or pft ==2):
              leafc_alloc      = 0.
              frootc_alloc     = 0.
              leafcstor_alloc  = availc * 1.0/fall
              frootcstor_alloc = availc * f1/fall
            else:
              leafcstor_alloc  = 0.
              frootcstor_alloc = 0.
              leafc_alloc      = availc * 1.0/fall
              frootc_alloc     = availc * f1/fall
            livestemc_alloc  = availc * flw*f2/fall
            deadstemc_alloc  = availc * (1.0-flw) * f2/fall
            livecrootc_alloc = availc * flw*(f2*f3)/fall
            deadcrootc_alloc = availc * (1.0-flw) * f2*f3/fall

            livestemc_litter =  parms['r_mort'] / 365.0 * livestemc[v]
            deadstemc_litter =  parms['r_mort'] / 365.0 * deadstemc[v]
            livecrootc_litter = parms['r_mort'] / 365.0 * livecrootc[v]
            deadcrootc_litter = parms['r_mort'] / 365.0 * deadcrootc[v]
            #Excess maintance respration taken from wood pools
            if (f2+f3 < 0):
              xsmr_deadstemc    = f2/(f2+f3)*xsmr
              xsmr_deadcrootc   = f3/(f2+f3)*xsmr
              xsmr_frootc       = 0.0
            else:
              xsmr_deadstemc    = 0.0
              xsmr_deadcrootc   = 0.0
              xsmr_frootc       = xsmr

            #Cacluate live wood turnover
            livestemc_turnover  = parms['lwtop_ann'] / 365. * livestemc[v]
            livecrootc_turnover = parms['lwtop_ann'] / 365. * livecrootc[v]
            cstor_turnover      = 1.0 / (parms['cstor_tau'] * 365) * cstor[v] * trate

            #increment plant pools
            leafc[v+1]       = leafc[v]       + leafc_alloc + leafc_trans - leafc_litter
            leafc_stor[v+1]  = leafc_stor[v]  + leafcstor_alloc - leafc_trans
            frootc[v+1]      = frootc[v]      + frootc_alloc + frootc_trans - frootc_litter - xsmr_frootc
            frootc_stor[v+1] = frootc_stor[v] + frootcstor_alloc - frootc_trans
            livestemc[v+1]   = livestemc[v]   + livestemc_alloc - livestemc_litter \
                                              - livestemc_turnover
            deadstemc[v+1]   = deadstemc[v]   + deadstemc_alloc - deadstemc_litter \
                                              + livestemc_turnover - xsmr_deadstemc
            livecrootc[v+1]  = livecrootc[v]  + livecrootc_alloc - livecrootc_litter \
                                              - livecrootc_turnover
            deadcrootc[v+1]  = deadcrootc[v]  + deadcrootc_alloc - deadcrootc_litter \
                                              + livecrootc_turnover - xsmr_deadcrootc
            cstor[v+1]       = cstor[v]       + cstor_alloc - cstor_turnover
            #Calculate growth respiration and NPP
            gr[v+1]          = availc * frg * (1.0 + f1+f2*(1+f3))/fall
            npp[v+1] = gpp[v+1] - mr[v+1] - gr[v+1]
            if (doy[v] == 1):
                annsum_npp = annsum_npp_temp
                annsum_npp_temp = 0
            annsum_npp_temp = annsum_npp_temp+npp[v]

            # ----------------- Litter and SOM decomposition model (CTC) --------------------
            trate = parms['q10_hr']**((0.5*(tmax[v]+tmin[v])-10)/10.0)

            ctc_input    = numpy.zeros([8],numpy.float)  #inputs to pool
            ctc_output   = numpy.zeros([8],numpy.float)  #Outputs from pool
            ctc_resp     = numpy.zeros([8],numpy.float)  #Respiration from pool
            #Litter inputs to the system
            ctc_input[0] = leafc_litter*parms['lf_flab'] + frootc_litter* \
                                parms['fr_flab']
            ctc_input[1] = leafc_litter*parms['lf_flig'] + frootc_litter* \
                                parms['fr_flig']+(livestemc_litter+deadstemc_litter \
                                +livecrootc_litter+deadcrootc_litter)*parms['cwd_flig']
            ctc_input[2] = leafc_litter*(1.0 - parms['lf_flab'] - parms['lf_flig']) + frootc_litter* \
                           (1.0-parms['fr_flab']-parms['fr_flig'])+(livestemc_litter+deadstemc_litter \
                           +livecrootc_litter+deadcrootc_litter)*(1.0 - parms['cwd_flig'])
            for p1 in range(0,8):
                if (p1 < 3):
                   ctc_output[p1] = k_ctc[p1]*ctcpools[p1,v]*trate*parms['fpi']
                else:
                   ctc_output[p1] = k_ctc[p1]*ctcpools[p1,v]*trate  #Decomposition (output)
                ctc_resp[p1]   = ctc_output[p1]*rf_ctc[p1]       #HR from this pool
                for p2 in range(0,8):
                    #Transfer (input) and respiration fluxes
                    ctc_input[p2] = ctc_input[p2] + ctc_output[p1]*tr_ctc[p1,p2]
            hr[v+1]=0
            for p in range(0,8):
                ctcpools[p,v+1] = ctcpools[p,v] + ctc_input[p] - ctc_output[p]
                hr[v+1] = hr[v+1] + ctc_resp[p]

            #Calculate NEE
            nee[v+1]    = hr[v+1] - npp[v+1]
            #Total system carbon
            totecosysc[v+1] = leafc[v+1]+leafc_stor[v+1]+frootc[v+1]+frootc_stor[v+1]+ \
                              livestemc[v+1]+deadstemc[v+1]+livecrootc[v+1]+deadcrootc[v+1] + \
                              sum(ctcpools[:,v+1])
            totlitc[v+1] = sum(ctcpools[0:3,v+1])
            totsomc[v+1] = sum(ctcpools[3:7,v+1])


    def run_selm(self, spinup_cycles=0, lat_bounds=[-999,-999], lon_bounds=[-999,-999], \
                     do_monthly_output=False, do_output_forcings=False, pft=-1,          \
                     prefix='model', ensemble=False, myoutvars=[], use_MPI=False):

        ens_torun=[]
        indx_torun=[]
        indy_torun=[]
        pfts_torun=[]
        pftfracs_torun=[]
        n_active=0
        if (self.site == 'none'):
         if (use_MPI):
           from mpi4py import MPI
           comm=MPI.COMM_WORLD
           rank=comm.Get_rank()
           size=comm.Get_size()
         else:
           rank = 0
           size = 0
         if (rank == 0):
          mydomain = Dataset(oscm_dir+"/models/pftdata/domain.360x720_ORCHIDEE0to360.100409.nc4",'r')
          landmask = mydomain.variables['mask']
          myinput = Dataset(oscm_dir+"/models/pftdata/surfdata_360x720_DALEC.nc4")
          pct_pft    = myinput.variables['PCT_NAT_PFT']
          pct_natveg = myinput.variables['PCT_NATVEG']
          self.hdlatgrid = myinput.variables['LATIXY']
          self.hdlongrid = myinput.variables['LONGXY']
          self.x1 = int(round((lon_bounds[0]-0.25)*2))
          if (self.x1 < 0):
             self.x1 = self.x1+720
          self.x2 = int(round((lon_bounds[1]-0.25)*2))
          if (self.x2 < 0):
             self.x2 = self.x2+720
          self.nx = self.x2-self.x1+1
          self.y1 = int(round((lat_bounds[0]+89.75)*2))
          self.y2 = int(round((lat_bounds[1]+89.75)*2))
          self.ny = self.y2-self.y1+1
          lats_torun=[]
          lons_torun=[]
          vegfrac_torun=[]
          if (self.ne > 1 and size > 1 and size < self.nx*self.ny):
            all_ensembles_onejob = True
            k_max = 1
          else:
            all_ensembles_onejob = False
            k_max = self.ne
          for i in range(0,self.nx):
              for j in range(0,self.ny):
                vegfrac    = pct_natveg[self.y1+j,self.x1+i]
                bareground = pct_pft[0,self.y1+j,self.x1+i]
                if (vegfrac > 0.1 and landmask[self.y1+j,self.x1+i] > 0):
                  if (bareground < 95.0):
                    mypfts=[]
                    mypftfracs=[]
                    if (pft < 0):
                      mypftfracs.append(sum(pct_pft[6:9,self.y1+j,self.x1+i])+pct_pft[3,self.y1+j,self.x1+i])
                      mypftfracs.append(sum(pct_pft[1:3,self.y1+j,self.x1+i])+pct_pft[4,self.y1+j,self.x1+i]+sum(pct_pft[9:12,self.y1+j,self.x1+i]))
                      mypftfracs.append(sum(pct_pft[12:,self.y1+j,self.x1+i]))
                    else:
                      mypftfracs=[0.0,0.0,0.0]
                      mypftfracs[pft] = 100.0
                    if (mypftfracs[0] > 5.0):
                      mypfts.append(0)
                    if (mypftfracs[1] > 5.0):
                      mypfts.append(1)
                    if (mypftfracs[2] > 5.0):
                      mypfts.append(2)
                    for k in range(0,k_max):
                      for p in mypfts:
                        lons_torun.append(self.hdlongrid[self.y1+j,self.x1+i])
                        lats_torun.append(self.hdlatgrid[self.y1+j,self.x1+i])
                        pfts_torun.append(p)
                        pftfracs_torun.append(mypftfracs[p])
                        indx_torun.append(i)
                        indy_torun.append(j)
                        ens_torun.append(k)
                        vegfrac_torun.append((100.0-bareground)/100.0)
                        n_active = n_active+1
          #Load all forcing data into memory
          self.get_regional_forcings()
          #get forcings for one point to get relevant info
          self.load_forcings(lon=lons_torun[0], lat=lats_torun[0])
        else:
          #site forcing has already been loaded
          all_ensembles_onejob = False
          rank = 0
          size = 0
          n_active = self.ne
          if (n_active > 1 and use_MPI):
             from mpi4py import MPI
             comm=MPI.COMM_WORLD
             rank=comm.Get_rank()
             size=comm.Get_size()
          if (rank == 0):
            for k in range(0,self.ne):
               pfts_torun.append(pft)
               pftfracs_torun.append(100.0)
               indx_torun.append(0)
               indy_torun.append(0)
               ens_torun.append(k)
            self.nx = 1
            self.ny = 1

        if (rank == 0):
          print('%d simulation units to run'%(n_active))
          n_done=0
          if (do_monthly_output):
             self.nt = (self.end_year-self.start_year+1)*12
             istart=0
          else:
             self.nt = int(self.end_year-self.start_year+1)*365
             istart=1

          model_output={}
          if (len(myoutvars) == 0):
            for v in self.outvars:
              if (v != 'ctcpools'):
                  myoutvars.append(v)
                  model_output[v] = numpy.zeros([self.ne,3,self.nt,self.ny,self.nx], numpy.float)
            for v in self.forcvars:
                if (v != 'time' and do_output_forcings):
                  myoutvars.append(v)
                  model_output[v] = numpy.zeros([3,self.nt,self.ny,self.nx], numpy.float)
          else:
             for v in myoutvars:
                 model_output[v] = numpy.zeros([self.ne,3,self.nt,self.ny,self.nx], numpy.float)
                 if (v in self.forcvars):
                     do_output_forcings=True
          self.pftfrac = numpy.zeros([self.ny,self.nx,3], numpy.float)

          if (self.site == 'none'):
            self.load_forcings(lon=lons_torun[0], lat=lats_torun[0])

          if ((n_active == 1 and self.ne == 1) or size == 0):
            #No MPI
            for i in range(0,n_active):
                if (self.site == 'none'):
                  self.load_forcings(lon=lons_torun[i], lat=lats_torun[i])
                if (self.ne > 1):
                  for p in range(0,len(self.ensemble_pnames)):
                    self.parms[self.ensemble_pnames[p]] = self.parm_ensemble[i,p]
                self.selm_instance(self.parms, spinup_cycles=spinup_cycles, pft=pfts_torun[i])
                self.pftfrac[indy_torun[i],indx_torun[i],pfts_torun[i]] = pftfracs_torun[i]
                for v in myoutvars:
                  if (v in self.outvars):
                    if (do_monthly_output):
                      model_output[v][ens_torun[i],pfts_torun[i],:,indy_torun[i],indx_torun[i]] = \
                         utils.daily_to_monthly(self.output[v][1:])
                    else:
                      model_output[v][ens_torun[i],pfts_torun[i],:,indy_torun[i],indx_torun[i]] = \
                         self.output[v][1:]
                  elif (v in self.forcvars):
                    if (do_monthly_output):
                       model_output[v][ens_torun,pfts_torun[i],:,indy_torun[i],indx_torun[i]] = \
                            utils.daily_to_monthly(self.forcings[v])
                    else:
                       model_output[v][ens_torun,pfts_torun[i],:,indy_torun[i],indx_torun[i]] = \
                            self.forcings[v][:]
            self.write_nc_output(model_output, do_monthly_output=do_monthly_output, prefix=prefix)
          else:
           #send first np-1 jobs where np is number of processes
           for n_job in range(1,size):
            comm.send(n_job, dest=n_job, tag=1)
            comm.send(0,     dest=n_job, tag=2)
            if (self.site == 'none'):
              self.load_forcings(lon=lons_torun[n_job-1], lat=lats_torun[n_job-1])
            parms = self.parms
            if (not all_ensembles_onejob and self.ne > 1):
              for p in range(0,len(self.ensemble_pnames)):
                parms[self.ensemble_pnames[p]] = self.parm_ensemble[ens_torun[n_job-1],p]

            comm.send(all_ensembles_onejob, dest=n_job, tag=300)
            comm.send(do_output_forcings,   dest=n_job, tag=400)
            comm.send(self.forcings,   dest = n_job, tag=6)
            comm.send(self.start_year, dest = n_job, tag=7)
            comm.send(self.end_year,   dest = n_job, tag=8)
            comm.send(self.nobs,       dest = n_job, tag=9)
            comm.send(self.lat,        dest = n_job, tag=10)
            comm.send(self.forcvars,   dest = n_job, tag=11)
            if (all_ensembles_onejob):
              comm.send(self.parm_ensemble,   dest=n_job, tag=100)
              comm.send(self.ensemble_pnames, dest=n_job, tag=101)
            else:
              comm.send(parms,           dest = n_job, tag=100)
            comm.send(myoutvars,         dest = n_job, tag=200)
            comm.send(pfts_torun[n_job-1], dest=n_job, tag=500)

           #Assign rest of jobs on demand
           for n_job in range(size,n_active+1):
            process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
            thisjob = comm.recv(source=process, tag=4)
            myoutput = comm.recv(source=process, tag=5)
            print('Received %d'%(thisjob))
            n_done = n_done+1
            comm.send(n_job, dest=process, tag=1)
            comm.send(0,     dest=process, tag=2)
            if (self.site == 'none'):
              self.load_forcings(lon=lons_torun[n_job-1], lat=lats_torun[n_job-1])
            if (not all_ensembles_onejob and self.ne > 1):
              for p in range(0,len(self.ensemble_pnames)):
                parms[self.ensemble_pnames[p]] = self.parm_ensemble[ens_torun[n_job-1],p]

            comm.send(all_ensembles_onejob, dest=process, tag=300)
            comm.send(do_output_forcings,   dest=process, tag=400)
            comm.send(self.forcings,   dest = process, tag=6)
            comm.send(self.start_year, dest = process, tag=7)
            comm.send(self.end_year,   dest = process, tag=8)
            comm.send(self.nobs,       dest = process, tag=9)
            comm.send(self.lat,        dest = process, tag=10)
            comm.send(self.forcvars,   dest = process, tag=11)
            if (all_ensembles_onejob):
              comm.send(self.parm_ensemble,   dest=process, tag=100)
              comm.send(self.ensemble_pnames, dest=process, tag=101)
            else:
              comm.send(parms,           dest = process, tag=100)
            comm.send(myoutvars,         dest = process, tag=200)
            comm.send(pfts_torun[n_job-1], dest=process, tag=500)
            #write output
            for v in myoutvars:
              if (all_ensembles_onejob):
                for k in range(0,self.ne):
                  model_output[v][k,pfts_torun[thisjob-1],:,indy_torun[thisjob-1],indx_torun[thisjob-1]] = \
                          myoutput[v][k,:]
              else:
                model_output[v][ens_torun[thisjob-1],pfts_torun[thisjob-1],:,indy_torun[thisjob-1],indx_torun[thisjob-1]] \
                                  = myoutput[v][0,:]
            self.pftfrac[indy_torun[thisjob-1],indx_torun[thisjob-1],pfts_torun[thisjob-1]] = pftfracs_torun[thisjob-1]

           #receive remaining messages and finalize
           while (n_done < n_active):
            process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
            thisjob = comm.recv(source=process, tag=4)
            myoutput = comm.recv(source=process, tag=5)
            vnum = 0
            print('Received %d'%(thisjob))
            n_done = n_done+1
            comm.send(-1, dest=process, tag=1)
            comm.send(-1, dest=process, tag=2)
            #write output
            for v in myoutvars:
              if (all_ensembles_onejob):
                for k in range(0,self.ne):
                  model_output[v][k,pfts_torun[thisjob-1],:,indy_torun[thisjob-1],indx_torun[thisjob-1]] = \
                        myoutput[v][k,:]
              else:
                model_output[v][ens_torun[thisjob-1],pfts_torun[thisjob-1],:,indy_torun[thisjob-1],indx_torun[thisjob-1]] \
                        = myoutput[v][0,:]
            self.pftfrac[indy_torun[thisjob-1],indx_torun[thisjob-1],pfts_torun[thisjob-1]] = pftfracs_torun[thisjob-1]
           self.write_nc_output(model_output, do_monthly_output=do_monthly_output, prefix=prefix)
           MPI.Finalize()
        #Slave
        else:
          status=0
          while status == 0:
            myjob = comm.recv(source=0, tag=1)
            status = comm.recv(source=0, tag=2)
            if (status == 0):
              all_ensembles_onejob = comm.recv(source=0, tag=300)
              do_output_forcings   = comm.recv(source=0, tag=400)
              self.forcings = comm.recv(source = 0, tag = 6)
              self.start_year = comm.recv(source=0, tag=7)
              self.end_year   = comm.recv(source=0, tag=8)
              self.nobs       = comm.recv(source=0, tag=9)
              self.lat        = comm.recv(source=0, tag=10)
              self.forcvars   = comm.recv(source=0, tag=11)
              if (all_ensembles_onejob):
                self.parm_ensemble = comm.recv(source=0, tag=100)
                self.ensemble_pnames = comm.recv(source=0, tag=101)
              else:
                myparms         = comm.recv(source=0, tag=100)
              myoutvars = comm.recv(source=0, tag=200)
              mypft     = comm.recv(source=0, tag=500)
              #Initialize output arrays
              self.output = {}
              self.output_ens = {}
              if (all_ensembles_onejob):
                 k_max = self.ne
              else:
                 k_max = 1
              for var in self.outvars:
                if (var == 'ctcpools'):
                   self.output[var] = numpy.zeros([8,self.nobs+1], numpy.float)
                else:
                   self.output[var] = numpy.zeros([self.nobs+1], numpy.float)

              thisoutput = {}
              thisoutput_ens = {}
              for k in range(0,k_max):
                if (all_ensembles_onejob):
                  myparms = self.pdefault
                  for p in range(0,len(self.ensemble_pnames)):
                    myparms[self.ensemble_pnames[p]] = self.parm_ensemble[k,p]
                self.selm_instance(myparms, spinup_cycles=spinup_cycles, pft=mypft )
                for v in myoutvars:
                   if (v in self.outvars):
                     if (do_monthly_output):
                         thisoutput[v] = utils.daily_to_monthly(self.output[v][1:])
                     else:
                         thisoutput[v] = self.output[v][1:]
                   elif (v in self.forcvars):
                     if (do_monthly_output):
                         thisoutput[v] = utils.daily_to_monthly(self.forcings[v])
                     else:
                         thisoutput[v] = self.forcings[v]
                   if (k == 0):
                      thisoutput_ens[v] = numpy.zeros([k_max, len(thisoutput[v])], numpy.float)
                   thisoutput_ens[v][k,:] = thisoutput[v]

              comm.send(rank,  dest=0, tag=3)
              comm.send(myjob, dest=0, tag=4)
              comm.send(thisoutput_ens, dest=0, tag=5)
          print('%d complete'%(rank))
          MPI.Finalize()

    def write_nc_output(self, output, do_monthly_output=False, prefix='model'):
         #set up output file
         output_nc = Dataset(prefix+'_output.nc', 'w', format='NETCDF4')
         output_nc.createDimension('pft',3)
         output_nc.createDimension('lon',self.nx)
         output_nc.createDimension('lat',self.ny)
         if (self.ne > 1):
           output_nc.createDimension('ensemble',self.ne)
           #ens_out = output_nc.createVariable('ensemble','i4',('ensemble',))
           #ens_out.axis="E"
           #ens_out.CoordinateAxisType = "Ensemble"
           pnum=0
           for p in self.ensemble_pnames:
             #write parameter values to file
             pvars={}
             pvars[p] = output_nc.createVariable(p, 'f8', ('ensemble',))
             pvars[p][:] = self.parm_ensemble[:,pnum]
             pnum=pnum+1
         if (self.site == 'none'):
           lat_out = output_nc.createVariable('lat','f8',('lat',))
           lon_out = output_nc.createVariable('lon','f8',('lon',))
           lat_out[:] = self.hdlatgrid[self.y1:self.y2+1,self.x1]
           lon_out[:] = self.hdlongrid[self.y1,self.x1:self.x2+1]
         else:
           lat_out = self.latdeg
           lon_out = self.londeg
         pft_out = output_nc.createVariable('pft_frac','f4',('lat','lon','pft'))
         for n in range(0,self.ne):
           pft_out[:,:,:] = self.pftfrac[:,:,:]
         if (do_monthly_output):
            output_nc.createDimension('time',(self.end_year-self.start_year+1)*12)
            time = output_nc.createVariable('time','f8',('time',))
            dpm = [31,28,31,30,31,30,31,31,30,31,30,31]
            time[0] = 15.5
            mlast = 0
            for i in range(1,(self.end_year-self.start_year+1)*12):
              time[i] = time[i-1]+(dpm[mlast]+dpm[i % 12])/2.0
              mlast = i % 12
            istart=0
         else:
            output_nc.createDimension('time',(self.end_year-self.start_year+1)*365)
            time = output_nc.createVariable('time','f8',('time',))
            for i in range(0,self.nobs):
              time[i] = i+0.5
            istart=1
         time.units='Days since '+str(self.start_year)+'-01-01 00:00'

         ncvars={}
         for v in output:
            if (self.ne > 1):
              ncvars[v] = output_nc.createVariable(v, 'f4',('ensemble','pft','time','lat','lon'))
              #for n in range(0,self.ne):
              #  for t in range(0,self.nt):

              ncvars[v][:,:,:,:,:] = output[v][:,:,:,:,:]
                  #ncvars[v][n,t,:,:] = (output[v][t,:,:,0,n]*self.pftfrac[:,:,0]/100.0 + \
                  #                      output[v][t,:,:,1,n]*self.pftfrac[:,:,1]/100.0 + \
                  #                      output[v][t,:,:,2,n]*self.pftfrac[:,:,2]/100.0).squeeze()
            else:
              ncvars[v] = output_nc.createVariable(v, 'f4',('pft','time','lat','lon',))
              for t in range(0,self.nt):
                 ncvars[v] = output[v][0,:,:,:,:]
#                ncvars[v][t,:,:] = (output[v][t,:,:,0]).squeeze()*pft_out[:,:,0]/100.0 + \
#                                   (output[v][t,:,:,1]).squeeze()*pft_out[:,:,1]/100.0 + \
#                                   (output[v][t,:,:,2]).squeeze()*pft_out[:,:,2]/100.0
         output_nc.close()

         # #output for eden vis system - customize as needed
         # eden_out = numpy.zeros([self.ne,pnum+1],numpy.float)
         # for n in range(0,self.ne):
         #   eden_out[n,0:pnum]      = self.parm_ensemble[n,:]
         #   eden_out[n,pnum:pnum+1] = numpy.mean(output['gpp'][n,0,0:60,0,0])*365.
         # numpy.savetxt("foreden.csv",eden_out,delimiter=",")

    def generate_synthetic_obs(self, parms, err):
        #generate synthetic observations from model with Gaussian error
        self.obs = numpy.zeros([self.nobs], numpy.float)
        self.obs_err = numpy.zeros([self.nobs], numpy.float)+err
        self.selm_instance(parms)
        for v in range(0,self.nobs):
            self.obs[v] = self.output[v]+numpy.random.normal(0,err,1)
        self.issynthetic = True
        self.actual_parms = parms

    def get_regional_forcings(self):
        print('Loading regional forcings')
        self.regional_forc={}
        fnames = ['TMAX','TMIN','FSDS','BTRAN']
        vnames = ['TSA', 'TSA', 'FSDS','BTRAN']
        fnum=0
        for f in fnames:
          os.system('mkdir -p '+oscm_dir+'/models/elm_drivers')
          driver_path = os.path.abspath(oscm_dir+'/models/elm_drivers')
          myfile = "GSWP3_fromELMSP_"+f+"_1980-2009.nc4"
          if (not os.path.exists(driver_path+'/'+myfile)):
            print('%s not found.  Downloading.'%(myfile))
            os.system('wget --no-check-certificate https://acme-webserver.ornl.gov/dmricciuto/elm_drivers/'+myfile)
            os.rename(myfile, driver_path+'/'+myfile)
          print('%s/%s'%(driver_path,myfile))
          myinput = Dataset(driver_path+'/'+myfile,'r')
          self.regional_forc[f.lower()] = myinput.variables[vnames[fnum]][:,:,:]
          myinput.close()
          fnum = fnum+1
        print ('Loading complete')

    def load_forcings(self, site='none', lat=-999, lon=-999):
         #Get single point data from E3SM style cpl_bypass input files
         self.forcvars = ['tmax','tmin','rad','cair','doy','dayl','btran','time']
         self.forcings = {}
         self.site=site
         if (self.site != 'none'):
             #Get data for requested site
             myinput = Dataset(oscm_dir+"/models/site_drivers/"+self.site+'_forcing.nc4','r',format='NETCDF4')
             npts = myinput.variables['TBOT'].size              #number of half hours or hours
             tair = myinput.variables['TBOT'][0,:]              #Air temperature (K)
             fsds  = myinput.variables['FSDS'][0,:]             #Solar radiation (W/m2)
             btran = fsds * 0.0 + 1.0
             self.latdeg = myinput.variables['LATIXY'][0]            #site latitude
             self.londeg = myinput.variables['LONGXY'][0]            #site longitude
             self.start_year = int(myinput.variables['start_year'][:])    #starting year of data
             self.end_year   = int(myinput.variables['end_year'][:])   #ending year of data
             npd = int(npts/(self.end_year - self.start_year + 1)/365)   #number of obs per day
             self.nobs = int((self.end_year - self.start_year + 1)*365)  #number of days
             myinput.close()
             for fv in self.forcvars:
               self.forcings[fv] = []
             self.lat = self.latdeg*numpy.pi/180.
             #populate daily forcings
             for d in range(0,self.nobs):
               self.forcings['tmax'].append(max(tair[d*npd:(d+1)*npd])-273.15)
               self.forcings['tmin'].append(min(tair[d*npd:(d+1)*npd])-273.15)
               self.forcings['rad'].append(sum(fsds[d*npd:(d+1)*npd]*(86400/npd)/1e6))
               self.forcings['btran'].append(1.0)
               self.forcings['cair'].append(360)
               self.forcings['doy'].append((float(d % 365)+1))
               self.forcings['time'].append(self.start_year+d/365.0)
               #Calculate day length
               dec  = -23.4*numpy.cos((360.*(self.forcings['doy'][d]+10.)/365.)*numpy.pi/180.)*numpy.pi/180.
               mult = numpy.tan(self.lat)*numpy.tan(dec)
               if (mult >= 1.):
                 self.forcings['dayl'].append(24.0)
               elif (mult <= -1.):
                  self.forcings['dayl'].append(0.)
               else:
                 self.forcings['dayl'].append(24.*numpy.arccos(-mult)/numpy.pi)

         elif (lat >= -90 and lon >= -180):
             #Get closest gridcell from reanalysis data
             self.latdeg=lat
             if (lon > 180):
                 lon=lon-360.
             if (lat > 9.5 and lat < 79.5 and lon > -170.5 and lon < -45.5):
                xg = int(round((lon + 170.25)*2))
                yg = int(round((lat - 9.75)*2))
                tmax = self.regional_forc['tmax'][:,yg,xg]
                tmin = self.regional_forc['tmin'][:,yg,xg]
                btran = self.regional_forc['btran'][:,yg,xg]
                fsds = self.regional_forc['fsds'][:,yg,xg]
             else:
                print('regions outside North America not currently supported')
                sys.exit(1)
             self.start_year = 1980
             self.end_year   = 2009
             npd = 1
             self.nobs = (self.end_year - self.start_year + 1)*365
             self.lat = self.latdeg*numpy.pi/180.

             #populate daily forcings
             self.forcings['tmax']  = tmax-273.15
             self.forcings['tmin']  = tmin-273.15
             self.forcings['btran'] = btran
             self.forcings['rad']   = fsds*86400/1e6
             self.forcings['cair']  = numpy.zeros([self.nobs], numpy.float) + 360.0
             self.forcings['doy']   = (numpy.cumsum(numpy.ones([self.nobs], numpy.float)) - 1) % 365 + 1
             self.forcings['time']  = self.start_year + (numpy.cumsum(numpy.ones([self.nobs], numpy.float)-1))/365.0
             self.forcings['dayl']  = numpy.zeros([self.nobs], numpy.float)
             for d in range(0,self.nobs):
               #Calculate day length
               dec  = -23.4*numpy.cos((360.*(self.forcings['doy'][d]+10.)/365.)*numpy.pi/180.)*numpy.pi/180.
               mult = numpy.tan(self.lat)*numpy.tan(dec)
               if (mult >= 1.):
                 self.forcings['dayl'][d] = 24.0
               elif (mult <= -1.):
                 self.forcings['dayl'][d] = 0.
               else:
                 self.forcings['dayl'][d] = 24.*numpy.arccos(-mult)/numpy.pi

         #define the x axis for plotting output (time, or one of the inputs)
         self.xlabel = 'Time (years)'
         #Initialize output arrays
         self.output = {}
         for var in self.outvars:
             if (var == 'ctcpools'):
                 self.output[var] = numpy.zeros([8,self.nobs+1], numpy.float)
             else:
                 self.output[var] = numpy.zeros([self.nobs+1], numpy.float)

    #Load actual observations and uncertainties
    def load_obs(self, site):
        obsvars = ['gpp','nee']
        myobs = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
        site_name = myobs.variables['site_name']
        lnum=0
        firstind = (self.start_year-1991)*365
        lastind  = (self.end_year-1991+1)*365
        self.obs={}
        nsm_0=site_name[:,:].shape[0]
        nsm_1=site_name[:,:].shape[1]
        siteName_str='US-UMB' #numpy.array([[site_name[i,j].decode("utf-8") for j in range(nsm_1)] for i in range(nsm_0)])
        #for s in site_name:
        for s in siteName_str:
          if site in ''.join(s):
             self.obs['gpp'] = myobs.variables['GPP'][lnum,firstind:lastind]*24*3600*1000
             self.obs['nee'] = myobs.variables['NEE'][lnum,firstind:lastind]*24*3600*1000
          lnum+=1

    #Plot model outputs, with obersvations for NEE and GPP if requested
    def plot_output(self, var='all', startyear=-1,endyear=-1, obs=True, figname_postfix=''):
      var_list = []
      if (var == 'all'):
         for key in self.output:
            var_list.append(key)
      else:
         var_list.append(var)
      for var in var_list:
          fig = plt.figure()
          ax = fig.add_subplot(111)
          if (var != 'ctcpools'):
              ax.plot(self.forcings['time'], self.output[var][1:])
              if (os.path.exists('./plots') == False):
                  os.mkdir('./plots')
              if obs and (var == 'gpp' or var == 'nee'):
                  ax.plot(self.forcings['time'], self.obs[var])
                  ax.legend(['Model','Observations'])
              if (startyear == -1):
                  startyear = self.start_year[0]
              if (endyear == -1):
                  endyear = self.end_year[0]
              plt.xlim(startyear, endyear+1)
              plt.xlabel('Year')
              plt.ylabel(var)
              plt.savefig('./plots/'+var+figname_postfix+'.pdf')

    def generate_ensemble(self, n_ensemble, pnames, fname='', normalized=False):
      self.parm_ensemble = numpy.zeros([n_ensemble,len(pnames)])
      self.ensemble_pnames = pnames
      self.ne = n_ensemble
      if (fname != ''):
        print('Generating parameter ensemble from %s'%(fname))
        pvals = numpy.genfromtxt(fname)
        if (normalized):
          for n in range(0,n_ensemble):
            for p in range(0,len(pnames)):
                self.parm_ensemble[n,p] = self.pmin[pnames[p]]+0.5* \
                     (pvals[n,p]+1.0)*(self.pmax[pnames[p]]-self.pmin[pnames[p]])
        else:
          for n in range(0,n_ensemble):
            for p in range(0,len(pnames)):
                self.parm_ensemble[n,p] = pvals[n,p]
        # inparms = open(fname,'r')
        # lnum = 0
        # for s in inparms:
        #   pvals = s.split()
        #   if (lnum < n_ensemble):
        #     for p in range(0,len(pnames)):
        #       if (normalized):
        #         self.parm_ensemble[lnum,p] = self.pmin[pnames[p]]+0.5* \
        #              (float(pvals[p]))*(self.pmax[pnames[p]]-self.pmin[pnames[p]])
        #       else:
        #         self.parm_ensemble[lnum,p] = float(pvals[p])
        #   lnum=lnum+1
        # inparms.close()
      else:
        for n in range(0,n_ensemble):
          # HACK: to have identical ensemble members uncomment line below
          # numpy.random.seed(2018)
          for p in range(0,len(pnames)):
            #Sample uniformly from the parameter space
            self.parm_ensemble[n,p] = numpy.random.uniform(low=self.pmin[pnames[p]], \
                  high=self.pmax[pnames[p]])
        #numpy.savetxt('inputs.txt', self.parm_ensemble)

