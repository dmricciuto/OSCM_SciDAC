import numpy

class MyModel(object):
    
    def __init__(self):
        self.name = 'DALEC'
        self.nparms = 22
        self.parm_names = ['gdd_min', 'gdd_max', 'laimax', 'tsmin', 'leaffall', 'nue', 'rg_frac', \
                           'br_mr', 'q10_mr', 'astem', 'troot', 'tstem', 'q10_hr', 'br_lit', \
                           'br_som', 'dr', 'lma', 'leafcn', 'stemc_init', 'rootc_init', 'litc_init', \
                           'somc_init']  #Parameter names
        self.pdef = [100.0, 200.0, 4.0,  5.0, 0.1,  7.0, 0.2, 0.00010, 2.0, 0.70,  5.0, 50.0, \
                     2.0, 3.0,   30.0, 0.00100,  80.0, 25.0, 5000.0,  500.0,   600.0,  7000.0]
        self.pmin = [  0.0, 150.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.00003, 1.0, 0.20,  2.0, 20.0, \
                     1.0, 1.0,   10.0, 0.00025,  20.0, 10.0, 1000.0,  100.0,    50.0,  1000.0]
        self.pmax = [150.0, 500.0, 7.0, 20.0, 1.0, 20.0, 0.6, 0.00030, 4.0, 0.95, 10.0, 100.0, \
                     4.0, 10.0, 100.0, 0.00400, 150.0, 50.0, 15000.0, 3000.0, 1000.0, 25000.0]
        self.pstd = numpy.zeros(self.nparms, numpy.float)
        self.issynthetic = False
        
        #units of the model output
        self.ylabel = 'NEE (umol m-2 s-1)'
    

    def run(self, parms):
        #initial pool sizes
        self.stemc[0] = parms[18]
        self.rootc[0] = parms[19]
        self.litrc[0] = parms[20]
        self.somc[0]  = parms[21]
        gdd = 0
        #Model parameters
        gdd_min = parms[0]
        gdd_max = parms[1]
        laimax  = parms[2]
        tsmin   = parms[3]
        leaffall= parms[4]
        nue     = parms[5]
        rg_frac = parms[6]
        br_mr   = parms[7]
        q10_mr  = parms[8]
        astem   = parms[9]
        troot   = 1/(parms[10]*365.0)
        tstem   = 1/(parms[11]*365.0)
        q10_hr  = parms[12]
        br_lit  = 1/(parms[13]*365.0)
        br_som  = 1/(parms[14]*365.0)
        dr      = parms[15]
        lma     = parms[16]
        leafcn  = parms[17]
        #Coeffs for Photosynthesis model
        a = [nue, 0.0156935, 4.22273, 208.868, 0.0453194, 0.37836, 7.19298, 0.011136, \
             2.1001, 0.789798]

        #Run the model
        for v in range(0,self.nobs):
            
            #Phenology
            gdd = (self.doy[v] > 1) * (gdd + max(0.5*(self.tmax[v]+self.tmin[v])-10.0, 0.0))
            if (self.doy[v] < 200):
                self.lai[v+1] = laimax * min((max(gdd,gdd_min) - gdd_min)/(gdd_max-gdd_min), 1.0)
            elif (self.tmin[v] < tsmin and self.lai[v] > 0):
                self.lai[v+1] = max(self.lai[v] - leaffall*laimax, 0.0)
            else:
                self.lai[v+1] = self.lai[v]
            dleaf = (self.lai[v+1]-self.lai[v])*lma
            
            #Calculate GPP flux
            if (self.lai[v] > 1e-3):
                rtot = 1.0
                psid = -2.0
                leafn = lma/leafcn
                gs = abs(psid)**a[9]/((a[5]*rtot+(self.tmax[v]-self.tmin[v])))
                pp = self.lai[v]*leafn/gs*a[0]*numpy.exp(a[7]*self.tmax[v])
                qq = a[2]-a[3]
                #internal co2 concentration
                ci = 0.5*(self.cair[v]+qq-pp+((self.cair[v]+qq-pp)**2-4.*(self.cair[v]*qq-pp*a[2]))**0.5)
                e0 = a[6]*self.lai[v]**2/(self.lai[v]**2+a[8])
                cps   = e0*self.rad[v]*gs*(self.cair[v]-ci)/(e0*self.rad[v]+gs*(self.cair[v]-ci))
                self.gpp[v+1] = cps*(a[1]*self.dayl[v]+a[4])
            else:
                self.gpp[v+1] = 0.0

            #Autotrophic repiration fluxes
            trate = q10_mr**((0.5*(self.tmax[v]+self.tmin[v])-10)/10.0)
            self.mr[v+1] = (self.leafc[v]+self.rootc[v])*br_mr*trate
            self.gr[v+1] = max(rg_frac*(self.gpp[v+1] - self.mr[v+1]), 0.0)
            if (dleaf > self.gpp[v+1] - self.mr[v+1] and dleaf > 0):
                self.gr[v+1] = self.gr[v+1] + (dleaf - max(self.gpp[v+1]-self.mr[v+1], 0.0))*rg_frac
            
            #Update vegetation and litter pools (allocation, litterfall and decomp)
            trate = q10_hr**((0.5*(self.tmax[v]+self.tmin[v])-10)/10.0)
            self.npp[v] = self.gpp[v] - self.mr[v] - self.gr[v]
            leaf_litter = -1.0 * min(dleaf, 0.0)
            stem_litter = tstem * self.stemc[v]
            root_litter = troot * self.rootc[v]
            self.leafc[v+1] = self.leafc[v] + dleaf
            self.stemc[v+1] = self.stemc[v] - stem_litter + (self.npp[v]-max(dleaf,0.0))*astem
            self.rootc[v+1] = self.rootc[v] - root_litter + (self.npp[v]-max(dleaf,0.0))*(1.0-astem)
            self.litrc[v+1] = self.litrc[v] + leaf_litter + stem_litter + root_litter - \
                                dr*self.litrc[v] - br_lit*self.litrc[v]*trate
            self.somc[v+1] = self.somc[v] + dr*self.litrc[v] - br_som*self.somc[v]*trate
            self.nee[v+1] = -1.0*(self.leafc[v+1] + self.stemc[v+1] + self.rootc[v+1] + self.litrc[v+1] + \
                           self.somc[v+1]) + (self.leafc[v] + self.stemc[v] + self.rootc[v] + \
                           self.litrc[v] + self.somc[v])
                                
        self.output = self.nee[1:]
                           
    def generate_synthetic_obs(self, parms, err):
        #generate synthetic observations from model with Gaussian error
        self.obs = numpy.zeros([self.nobs], numpy.float)
        self.obs_err = numpy.zeros([self.nobs], numpy.float)+err
        self.run(parms)
        for v in range(0,self.nobs):
            self.obs[v] = self.output[v]+numpy.random.normal(0,err,1)
        self.issynthetic = True
        self.actual_parms = parms

    def load_forcings(self, site):
         myinput = open('./site_drivers/'+site+'_drivers.csv','r')
         self.tmin = []
         self.tmax  = []
         self.rad = []
         self.cair = []
         self.doy = []
         self.dayl = []
         self.x = []
         self.nobs = 0
         for s in myinput:
             if (self.nobs < 365):
                 self.tmin.append(float(s.split(',')[1])-273.15)
                 self.tmax.append(float(s.split(',')[2])-273.15)
                 self.rad.append(float(s.split(',')[3]))
                 self.cair.append(float(s.split(',')[4]))
                 self.doy.append((float(self.nobs % 365)+1))
                 self.x.append((self.nobs)/365.0)
                 self.nobs = self.nobs+1
         myinput.close()
         self.lat = 45.0
         self.lat = self.lat*numpy.pi/180.
         for d in range(0,self.nobs):
             #Calculate day length
             dec  = -23.4*numpy.cos((360.*(self.doy[d]+10.)/365.)*numpy.pi/180.)*numpy.pi/180.
             mult = numpy.tan(self.lat)*numpy.tan(dec)
             if (mult >= 1.):
                 self.dayl.append(24.0)
             elif (mult <= -1.):
                 self.dayl.append(0.)
             else:
                 self.dayl.append(24.*numpy.arccos(-mult)/numpy.pi)

         #define the x axis for plotting output (time, or one of the inputs)
         self.xlabel = 'Time (years)'
         #initialize the output arrays (speedier than doing here this every model run...)
         self.output = numpy.zeros([self.nobs], numpy.float)
         self.lai = numpy.zeros([self.nobs+1], numpy.float)
         self.gpp = numpy.zeros([self.nobs+1], numpy.float)
         self.npp = numpy.zeros([self.nobs+1], numpy.float)
         self.gr = numpy.zeros([self.nobs+1], numpy.float)
         self.mr = numpy.zeros([self.nobs+1], numpy.float)
         self.leafc = numpy.zeros([self.nobs+1], numpy.float)
         self.stemc = numpy.zeros([self.nobs+1], numpy.float)
         self.rootc = numpy.zeros([self.nobs+1], numpy.float)
         self.litrc = numpy.zeros([self.nobs+1], numpy.float)
         self.somc = numpy.zeros([self.nobs+1], numpy.float)
         self.nee = numpy.zeros([self.nobs+1], numpy.float)

    #Load actual observations and uncertainties here (or in the init function, if that is easier)
    def load_obs(self):
        self.obs = [0.2,0.3,0.4,0.5,0.8,1.1,1.3,1.7,2.1,2.0,2.5,3.1,3.8,4.2,4.4]
        self.obs_err = numpy.zeros([self.nobs], numpy.float) + 0.3 
        self.issynthetic = False
            
