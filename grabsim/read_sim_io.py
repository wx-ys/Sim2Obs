from abc import ABC, abstractmethod

import h5py
import numpy as np
from scipy.spatial import KDTree
from pynbody.array import SimArray

#import configparser

'''
source:
normal-star:  age > 30 Myr, BruzualCharlotSEDFamily
star-forming: age 0 ~ 30 Myr: ToddlersSEDFamily, age 0 ~ 10 Myr: MappingsSEDFamily, 

medium:
dust:         temp < 1e4 K,
electrons:    temp >1e4 K, 
gas:          
'''

class Read_sim(ABC):
    @abstractmethod
    def get_star(self,key):
        # key: 
        # x,y,z,        coordinate,         (kpc)
        # h,            smoothing length    (kpc)
        # vx,vy,vz,     bulk velocity       (km/s)
        # vdisp         SubfindVelDisp      (km/s)
        
        # mass          the mass            (Msun)
        
        # mbirth    the initial mass    (Msun)
        # metals        metallicity         (1)
        # age           age                 (Gyr)
        
        # BpassSEDFamily,           mbirth, metals, age
        # BruzualCharlotSEDFamily,  mbirth, metals, age     **
        # CastelliKuruczSEDFamily,  stellar_radius, effective stellar surface temperature, metals, gravity at the stellar surface
        # FileSSPSEDFamily,         mbirth, metals, age, [U(1) hasIonizationParameter]
        # FSPSSEDFamily,            mbirth, metals, age
        # MarastonSEDFamily         mbirth, metals, age
        # Starburst99SEDFamily      mbirth, metals, age
        
        
        
        # ToddlersSEDFamily         age, metals, sfe (1), ncloud, mbirth
        # represents the TODDLERS (Time evolution of Dust Diagnostics and Line Emission from Regions containing young Stars) 
        # family of star-forming region SEDs
            # ages = 0.1 ~ 30Myr\
            # metals = 0.001 ~ 0.04
            # sfe = 0.01 ~ 0.15
            # ncloud, cloud number densities, = 10 ~ 2560 cm-3
            # mbirth
            
            
        # MappingsSEDFamily         sfr, metals, compactness, pressure, PDR covering factor,
        # represents the family of MAPPINGS III star-forming region template SEDs
            # sfr, mbirth/ 10Myr, (Msun/yr)
            # metals, 
            # logC, compactness,
            # logCompactnessMean = 5, logCompactnessStd = 0.4, # from Kapoor et al. 2021
                    #logC = np.random.normal(loc=np.float32(self.config['logCompactnessMean']), 
                    #  scale=np.float32(self.config['logCompactnessStd']), size=size) # from Kapoor et al. 2021
            # pressure, logPressure = 5, log10[(PRead_AnastrisTNGressure/k_B)/cm^-3 K] = logPressure
            # pressure = (10**np.float32(self.config['logPressure']) * const.k_B) * (u.J * u.cm**-3).to(u.J * u.m**-3) # J * m**3 == Pa
            # coveringFactor = 0.2, from Groves et al. 2008
        pass
    
    @abstractmethod
    def get_gas(self,key):
        # key: 
        # x,y,z,        coordinate,         (kpc)
        # h,            smoothing length    (kpc)
        # mass,electrons   the mass            (Msun or 1)
        # metals        metallicity         (1)
        # temp          Temperature         (K)
        # vx,vy,vz,     bulk velocity       (km/s)
        # Bx,By,Bz,     magnetic field      (G)
        pass
    
    @abstractmethod
    def load_init(self):

        pass

    

class Read_AnastrisTNG(Read_sim):
    try:
        import AnastrisTNG
    except ImportError:
        print('no AnastrisTNG')
    _gas_transkey = {'Bx': 'MagneticField_x',
                           'By': 'MagneticField_y',
                           'Bz': 'MagneticField_z',
                           'electrons': 'ElectronAbundance',
            }
    _star_transkey = {'Bx': 'MagneticField_x',
                           'By': 'MagneticField_y',
                           'Bz': 'MagneticField_z',
                           'mbirth': 'GFM_InitialMass',
            }
    def __init__(self,config,): #config.get['SimIO','file_operator']
        
        self._config = config
        self.filepath = self._config.get('SimIO','SIM_FILE_PATH')
        
        self.star_nth = self._config.getint('SimIO','SIM_STAR_SMOOTH_NTH',fallback = 32)
        self.star_smooth_max = self._config.getfloat('SimIO','SIM_STAR_SMOOTH_MAX',fallback = 0.8)
        self.gas_smooth = self._config.getfloat('SimIO','SIM_GAS_SMOOTH', fallback = 1)
        self.gas_smooth_max = self._config.getfloat('SimIO','SIM_GAS_SMOOTH_MAX', fallback = 3)
        self.snap = self._config.getint('SimIO','SIM_SNAP', fallback = 99)
        import AnastrisTNG
        self.snapshot = AnastrisTNG.TNGsimulation.Snapshot(BasePath=self.filepath,Snap=self.snap)
        self.snapshot.load_particle_para['star_fields']=self._config.get('SimIO','SIM_STAR_FIELDS').split(',')
        self.snapshot.load_particle_para['gas_fields']=self._config.get('SimIO','SIM_GAS_FIELDS').split(',')
        self._data={}
        self._data_kdtree={}
        #self.load_init()
        
    def load_init(self):
        if 'SIM_SUBHALOID' in self._config['SimIO']:
            if self._config['SimIO']['SIM_SUBHALOID']:
                self.load_particle(self._config.getint('SimIO','SIM_SUBHALOID'),groupType='Subhalo')
                self.groupType='Subhalo'
                self.ID = self._config.getint('SimIO','SIM_SUBHALOID')
        if 'SIM_HALOID' in self._config['SimIO']:
            if self._config['SimIO']['SIM_HALOID']:
                self.load_particle(self._config.getint('SimIO','SIM_HALOID'),groupType='Halo')
                self.groupType='Halo'
                self.ID = self._config.getint('SimIO','SIM_HALOID')
                
    def _clear(self):
        self._data={}
        self._data_kdtree={}
        
    def load_particle(self,ID,groupType='Subhalo'):
        if f'{groupType}_{ID}' in self._data:
            pass
        else:
            print(f'Loading {groupType}_{ID} from {self.snapshot} \n')
            self._data[f'{groupType}_{ID}'] = self.snapshot.load_particle(ID,groupType)
            self._data[f'{groupType}_{ID}'].physical_units()
            self._data[f'{groupType}_{ID}'].face_on()
            self._print_properties(groupType,ID)
            
    def _print_properties(self,groupType,ID):
        print(f'{groupType}_{ID} properties:')
        Ms = self._data[f'{groupType}_{ID}'].s['mass'].sum()
        print(f'--- log(Mstar): {np.log10(Ms+1):.2f} {Ms.units}')
        Mg = self._data[f'{groupType}_{ID}'].g['mass'].sum()
        print(f'--- log(Mgas): {np.log10(Mg+1):.2f} {Mg.units}')
        krot = self._data[f'{groupType}_{ID}'].krot()
        print(f'--- krot: {krot:.2f}')
        re = self._data[f'{groupType}_{ID}'].re
        print(f'--- re: {re:.2f} {re.units} \n')
        
        
    def get_star(self,key,**kwargs):
        
        self.load_init()
        
        ID = kwargs.get('ID',self.ID)
        groupType = kwargs.get('groupType', self.groupType)
        
        particle = kwargs.get('particle',None)
        
        self.load_particle(ID,groupType)
        NEED_PA = {}
        NEED_PA['mask'] = np.ones(len(self._data[f'{groupType}_{ID}'].s),dtype=np.bool_)
        if particle:
            cons = self._config[particle]
            for i in cons:
                isp = i.split('_')
                if isp[0] =='s':
                    if len(isp) == 2:
                        NEED_PA[isp[-1]] = float(cons[i])
                    if len(isp) == 3:
                        if isp[-1].lower() == 'max':
                            NEED_PA['mask'] = (NEED_PA['mask']) & (self._data[f'{groupType}_{ID}'].s[isp[1]] < float(cons[i]))
                        if isp[-1].lower() == 'min':
                            NEED_PA['mask'] = (NEED_PA['mask']) & (self._data[f'{groupType}_{ID}'].s[isp[1]] > float(cons[i]))
                            
        NEED_PA['mask'] = (NEED_PA['mask']) & (self._data[f'{groupType}_{ID}'].s['age'] > 0)
        if key in NEED_PA:
            return NEED_PA[key]*np.ones(len(self._data[f'{groupType}_{ID}'].s))[NEED_PA['mask']]
        if key == 'sfr':
            return (self._data[f'{groupType}_{ID}'].s['mass']/float(cons['s_age_max'])/1e9)[NEED_PA['mask']] # Msun/yr
        
        
        
        if key == 'h':
            if f'{groupType}_{ID}_star_h' in self._data_kdtree:
                return self._data_kdtree[f'{groupType}_{ID}_star_h'][NEED_PA['mask']]
            kdtree_data = KDTree(self._data[f'{groupType}_{ID}'].s['pos'])
            self._data_kdtree[f'{groupType}_{ID}_star'] = kdtree_data
            distances, _ = kdtree_data.query(self._data[f'{groupType}_{ID}'].s['pos'], k=self.star_nth) # smooth length
            h = distances[:, -1]
            h[h>self.star_smooth_max] = self.star_smooth_max
            h = SimArray(h,units=self._data[f'{groupType}_{ID}'].s['pos'].units)
            self._data_kdtree[f'{groupType}_{ID}_star_h'] = h
            return h[NEED_PA['mask']]
        
        
        if key in self._star_transkey:
            return self._data[f'{groupType}_{ID}'].s[self._star_transkey[key]][NEED_PA['mask']]
        
        return self._data[f'{groupType}_{ID}'].s[key][NEED_PA['mask']]
    
    def get_gas(self,key,**kwargs):
        
        self.load_init()
        
        ID = kwargs.get('ID',self.ID)
        groupType = kwargs.get('groupType', self.groupType)
        self.load_particle(ID,groupType)
        
        particle = kwargs.get('particle',None)

        
        NEED_PA = {}
        NEED_PA['mask'] = np.ones(len(self._data[f'{groupType}_{ID}'].g),dtype=np.bool_)
        if particle:
            cons = self._config[particle]
            for i in cons:
                isp = i.split('_')
                if isp[0] =='m':
                    if len(isp) == 2:
                        NEED_PA[isp[-1]] = float(cons[i])
                    if len(isp) == 3:
                        if isp[-1].lower() == 'max':
                            NEED_PA['mask'] = (NEED_PA['mask']) & (self._data[f'{groupType}_{ID}'].g[isp[1]] < float(cons[i]))
                        if isp[-1].lower() == 'min':
                            NEED_PA['mask'] = (NEED_PA['mask'] )&( self._data[f'{groupType}_{ID}'].g[isp[1]] > float(cons[i]))
        if key in NEED_PA:
            return NEED_PA[key]*np.ones(len(self._data[f'{groupType}_{ID}'].g))[NEED_PA['mask']]
        
        if key == 'h':
            h = np.power(self._data[f'{groupType}_{ID}'].g['mass']/self._data[f'{groupType}_{ID}'].g['rho']/(4*np.pi),1/3)
            h = h*self.gas_smooth
            h.units = self._data[f'{groupType}_{ID}'].g['pos'].units
            h[h>self.gas_smooth_max] = self.gas_smooth_max
            return h[NEED_PA['mask']]
        if key in self._gas_transkey:
            return self._data[f'{groupType}_{ID}'].g[self._gas_transkey[key]][NEED_PA['mask']]
        
        return self._data[f'{groupType}_{ID}'].g[key][NEED_PA['mask']]
    
    