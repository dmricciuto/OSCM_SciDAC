import numpy
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from math import sin, cos, sqrt, atan2, radians
import os

class MyModel(object):
    
    def __init__(self):
        self.name = 'DALEC'
        self.parms = {'gdd_crit': 50.0, 'crit_dayl': 10.917, 'ndays_on':30, 'ndays_off': 15,           \
                      'nue': 15.0, 'slatop':0.03,                                                      \
                      'livewdcn': 50, 'leafcn': 25, 'frootcn': 42,                                     \
                      'fstor2tran': 0.5, 'stem_leaf': 2.7, 'croot_stem': 0.3, 'f_livewd':0.1,          \
                      'froot_leaf': 0.5,                                                               \
                      'rg_frac': 0.3, 'br_mr': 0.216, 'q10_mr': 1.5, 'cstor_tau':3.0,                  \
                      'r_mort': 0.02, 'lwtop_ann': 0.7, 'leaf_long': 1.5, 'froot_long': 1.5,           \
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
        self.pdefault=self.parms
        self.pmin = {}
        self.pmax = {}
        self.nparms = 0
        for p in self.parms:
            self.pmin[p] = self.parms[p]*0.75
            self.pmax[p] = self.parms[p]*1.25
            self.nparms = self.nparms+1
        self.issynthetic = False

        #Model outputs
        self.outvars = ['gpp','npp','gr', 'mr','hr','nee','lai','leafc','leafc_stor','frootc','frootc_stor','livestemc', \
                        'deadstemc','livecrootc','deadcrootc','ctcpools','totecosysc','totsomc','totlitc', 'cstor']


    def run(self, spinup_cycles=0, deciduous=True):

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

        #Set initial States (if non-zero)
        if (deciduous):
            leafc_stor[0] = 10.0
        else:
            leafc[0] = 10.0
        #Set initial soil carbon for long-lived pool
        ctcpools[6] = self.parms['soil4ci']

        #Forcings
        tmax = self.forcings['tmax']
        tmin = self.forcings['tmin']
        rad  = self.forcings['rad']
        doy  = self.forcings['doy']
        cair = self.forcings['cair']
        dayl = self.forcings['dayl']

        #Coefficents for ACM (GPP submodel)
        a = [self.parms['nue'], 0.0156935, 4.22273, 208.868, 0.0453194, 0.37836, 7.19298, 0.011136, \
             2.1001, 0.789798]
             
        #Turnover times for CTC model
        k_ctc = [self.parms['k_l1'],self.parms['k_l2'],self.parms['k_l3'],self.parms['k_s1'], \
                 self.parms['k_s2'],self.parms['k_s3'],self.parms['k_s4'],self.parms['k_frag']]
        #Respiration fractions for CTC model pools
        rf_ctc = [self.parms['rf_l1s1'], self.parms['rf_l2s2'], self.parms['rf_l3s3'] , \
          self.parms['rf_s1s2'], self.parms['rf_s2s3'], self.parms['rf_s3s4'], 1.0, 0.0]
        #transfer matrix for CTC model
        tr_ctc = numpy.zeros([8,8],numpy.float)
        tr_ctc[0,3] = 1.0 - self.parms['rf_l1s1']
        tr_ctc[1,4] = 1.0 - self.parms['rf_l2s2']
        tr_ctc[2,5] = 1.0 - self.parms['rf_l3s3']
        tr_ctc[3,4] = 1.0 - self.parms['rf_s1s2']
        tr_ctc[4,5] = 1.0 - self.parms['rf_s2s3']
        tr_ctc[5,6] = 1.0 - self.parms['rf_s3s4']
        tr_ctc[7,1] = self.parms['cwd_flig']
        tr_ctc[7,2] = 1.0 - self.parms['cwd_flig']


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
                                  
                                  #print s, (totecosysc[0]-totecosysc_last)/(self.nobs/365)
          for v in range(0,self.nobs):
            # --------------------1.  Phenology -------------------------
            #Calculate leaf on
            leafc_trans = 0.0
            frootc_trans = 0.0
            if (deciduous):     #Decidous phenology
              gdd_last = gdd
              dayl_last = dayl[v-1]
              gdd = (doy[v] > 1) * (gdd + max(0.5*(tmax[v]+tmin[v])-10.0, 0.0))
              if (gdd >= self.parms['gdd_crit'] and gdd_last < self.parms['gdd_crit']):
                  leafon = self.parms['ndays_on']
                  leafc_trans_tot  = leafc_stor[v]*self.parms['fstor2tran']
                  frootc_trans_tot = frootc_stor[v]*self.parms['fstor2tran']
              if (leafon > 0):
                  leafc_trans  = leafc_trans_tot  /self.parms['ndays_on']
                  frootc_trans = frootc_trans_tot /self.parms['ndays_on']
                  leafon = leafon - 1
              #Calculate leaf off
              if (dayl_last >= self.parms['crit_dayl'] and dayl[v] < self.parms['crit_dayl']):
                   leafoff       = self.parms['ndays_off']
                   leafc_litter_tot  = leafc[v]
                   frootc_litter_tot = frootc[v]
              if (leafoff > 0):
                   leafc_litter  = leafc_litter_tot /self.parms['ndays_off']
                   frootc_litter = frootc_litter_tot/self.parms['ndays_off']
                   leafoff = leafoff - 1
              else:
                   leafc_litter  = 0.0
                   frootc_litter = 0.0
            else:               #Evergreen phenology
                leafc_litter  = leafc[v]  * self.parms['leaf_long']/365.
                frootc_litter = frootc[v] * self.parms['froot_long']/365.

            lai[v+1] = leafc[v] * self.parms['slatop']

            #---------------------2. GPP -------------------------------------
            #Calculate GPP flux using the ACM model (Williams et al., 1997)
            if (lai[v] > 1e-3):
                rtot = 1.0
                psid = -2.0
                leafn = 1.0/(self.parms['leafcn'] * self.parms['slatop'])
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

            #--------------------3.  Maintenace respiration ------------------------
            #Maintenance respiration
            trate = self.parms['q10_mr']**((0.5*(tmax[v]+tmin[v])-25.0)/25.0)
            mr[v+1] = (leafc[v]/self.parms['leafcn'] + frootc[v]/self.parms['frootcn'] + \
                       (livecrootc[v]+livestemc[v])/self.parms['livewdcn'])* \
                       self.parms['br_mr']*trate
            #Nutrient limitation
            availc      = max(gpp[v+1]-mr[v+1],0.0)
            availc      = availc * self.parms['fpg']
            cstor_alloc = availc * (1.0 - self.parms['fpg'])
            xsmr        = max(mr[v+1]-gpp[v+1],0.0)

            #---------------4.  Allocation and growth respiration -------------------
            frg  = self.parms['rg_frac']
            flw  = self.parms['f_livewd']
            f1   = self.parms['froot_leaf']
            f2   = max(self.parms['stem_leaf']/(1.0+numpy.exp(-0.004*(annsum_npp - \
                           300.0))) - 0.4, 0.1)
            f3   = self.parms['croot_stem']
            fall = (1.0+frg)*(1.0 + f1 + f2*(1+f3))
            if (deciduous):
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

            livestemc_litter = self.parms['r_mort'] / 365.0 * livestemc[v]
            deadstemc_litter = self.parms['r_mort'] / 365.0 * deadstemc[v]
            livecrootc_litter =self.parms['r_mort'] / 365.0 * livecrootc[v]
            deadcrootc_litter =self.parms['r_mort'] / 365.0 * deadcrootc[v]
            #Excess maintance respration taken from wood pools
            xsmr_deadstemc    = f2/(f2+f3)*xsmr
            xsmr_deadcrootc   = f3/(f2+f3)*xsmr
            
            #Cacluate live wood turnover
            livestemc_turnover  = self.parms['lwtop_ann'] / 365. * livestemc[v]
            livecrootc_turnover = self.parms['lwtop_ann'] / 365. * livecrootc[v]
            cstor_turnover      = 1.0 / (self.parms['cstor_tau'] * 365) * cstor[v] * trate
 
            #increment plant pools
            leafc[v+1]       = leafc[v]       + leafc_alloc + leafc_trans - leafc_litter
            leafc_stor[v+1]  = leafc_stor[v]  + leafcstor_alloc - leafc_trans
            frootc[v+1]      = frootc[v]      + frootc_alloc + frootc_trans - frootc_litter
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
            trate = self.parms['q10_hr']**((0.5*(tmax[v]+tmin[v])-10)/10.0)

            ctc_input    = numpy.zeros([8],numpy.float)  #inputs to pool
            ctc_output   = numpy.zeros([8],numpy.float)  #Outputs from pool
            ctc_resp     = numpy.zeros([8],numpy.float)  #Respiration from pool
            #Litter inputs to the system
            ctc_input[0] = leafc_litter*self.parms['lf_flab'] + frootc_litter* \
                                self.parms['fr_flab']
            ctc_input[1] = leafc_litter*self.parms['lf_flig'] + frootc_litter* \
                                self.parms['fr_flig']+(livestemc_litter+deadstemc_litter \
                                +livecrootc_litter+deadcrootc_litter)*self.parms['cwd_flig']
            ctc_input[2] = leafc_litter*(1.0 - self.parms['lf_flab'] - self.parms['lf_flig']) + frootc_litter* \
                           (1.0-self.parms['fr_flab']-self.parms['fr_flig'])+(livestemc_litter+deadstemc_litter \
                           +livecrootc_litter+deadcrootc_litter)*(1.0 - self.parms['cwd_flig'])
            for p1 in range(0,8):
                if (p1 < 3):
                   ctc_output[p1] = k_ctc[p1]*ctcpools[p1,v]*trate*self.parms['fpi']
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
 
    def generate_synthetic_obs(self, parms, err):
        #generate synthetic observations from model with Gaussian error
        self.obs = numpy.zeros([self.nobs], numpy.float)
        self.obs_err = numpy.zeros([self.nobs], numpy.float)+err
        self.run(parms)
        for v in range(0,self.nobs):
            self.obs[v] = self.output[v]+numpy.random.normal(0,err,1)
        self.issynthetic = True
        self.actual_parms = parms

    def load_forcings(self, site='none', lat=-999, lon=-999):
         #Get data from E3SM style cpl_bypass input files
         if (site != 'none'):
             #Get data for requested site
             myinput = Dataset("../site_drivers/"+site+'_forcing.nc4','r',format='NETCDF4')
             npts = myinput.variables['TBOT'].size              #number of half hours or hours
             tair = myinput.variables['TBOT'][0,:]              #Air temperature (K)
             fsds  = myinput.variables['FSDS'][0,:]             #Solar radiation (W/m2)
             latdeg = myinput.variables['LATIXY'][0]            #site latitude
             londeg = myinput.variables['LONGXY'][0]            #site longitude
             self.start_year = myinput.variables['start_year'][:]    #starting year of data
             self.end_year   = myinput.variables['end_year'][:]      #ending year of data
             npd = npts/(self.end_year[0] - self.start_year[0] + 1)/365   #number of obs per day
             self.nobs = (self.end_year[0] - self.start_year[0] + 1)*365  #number of days
             myinput.close()
         elif (lat >= -90 and lon >= -180):
             #Get closest gridcell from reanalysis data
             latdeg=lat
             zone_map = open("../global_drivers/zone_mappings.txt",'r')
             mindist = 99999
             myzone  = 0
             mygcell = 0
             for s in zone_map:
                data=s.split()
                mylon = float(data[0]) 
                if mylon > 180:
                  mylon-mylon-360.0
                mylat = float(data[1])
                zone  = int(data[2])
                gcell = int(data[3])
                a = (sin(radians(lat)-radians(mylat))/2.0)**2 + cos(radians(lat)) * cos(radians(mylat)) \
                    * sin((radians(lon)-radians(mylon))/2.0)**2
                c = 2 * atan2(sqrt(a), sqrt(1-a))
                dist = 6373.0 * c               
                if dist < mindist:
                  mindist=dist
                  myzone = str(100+zone)
                  mygcell=gcell
             zone_map.close()
             mytinput = Dataset("../global_drivers/GSWP3.v1_TBOT_1901-2010_z" +myzone[1:]+ '.nc','r')  
             tair = mytinput.variables['TBOT'][mygcell,:]
             myrinput = Dataset("../global_drivers/GSWP3.v1_FSDS_1901-2010_z" +myzone[1:]+ '.nc','r')
             fsds = myrinput.variables['FSDS'][mygcell,:]
             mytinput.close()
             myrinput.close()
             self.start_year = 1901
             self.end_year   = 2010
             npd = 8
             self.nobs = (self.end_year - self.start_year + 1)*365

         forcvars = ['tmax','tmin','rad','cair','doy','dayl','time']
         self.forcings = {}
         for fv in forcvars:
             self.forcings[fv] = []
         self.lat = latdeg*numpy.pi/180.
         #populate daily forcings
         for d in range(0,self.nobs):
             self.forcings['tmax'].append(max(tair[d*npd:(d+1)*npd])-273.15)
             self.forcings['tmin'].append(min(tair[d*npd:(d+1)*npd])-273.15)
             self.forcings['rad'].append(sum(fsds[d*npd:(d+1)*npd]*(86400/npd)/1e6))
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
        myobs = Dataset("../site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
        site_name = myobs.variables['site_name']
        lnum=0
        firstind = (self.start_year-1991)*365
        lastind  = (self.end_year-1991+1)*365
        self.obs={}
        for s in site_name:
          if site in ''.join(s):
             self.obs['gpp'] = myobs.variables['GPP'][lnum,firstind:lastind]*24*3600*1000
             self.obs['nee'] = myobs.variables['NEE'][lnum,firstind:lastind]*24*3600*1000
          lnum+=1




    #Plot model outputs, with obersvations for NEE and GPP if requested
    def plot_output(self, var='all', startyear=-1,endyear=-1, obs=True):
      var_list = []
      if (var == 'all'):
         for key in self.output:
            var_list.append(key)
      else:  
         var_list.append(var)
      print var_list
      for var in var_list:
          fig = plt.figure()
          ax = fig.add_subplot(111)
          if (var != 'ctcpools'):
              ax.plot(self.forcings['time'], self.output[var][1:])
              if (os.path.exists('./plots') == False):
                  os.mkdir('./plots')
              if (var == 'gpp' or var == 'nee' and obs):
                  ax.plot(self.forcings['time'], self.obs[var])
                  ax.legend(['Model','Observations'])
              if (startyear == -1):
                  startyear = self.start_year[0]
              if (endyear == -1):
                  endyear = self.end_year[0]
              plt.xlim(startyear, endyear+1)
              plt.xlabel('Year')
              plt.ylabel(var)
              plt.savefig('./plots/'+var+'.pdf')
