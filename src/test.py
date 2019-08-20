import numpy as np
import maces_utilities as utils
import random
import importlib
from datetime import date
from TAIHydroMOD import tai_hydro_mod as taihydro

# https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2012JF002363
# Test case: Venice Lagoon

def icdf_wind_gen(x):
    if x>=0 and x<0.57071:
        return 2.5 / 0.57071 * x
    elif x>=0.57071 and x<0.83141:
        return 2.5 + 2.5 / 0.2607 * (x - 0.57071)
    elif x>=0.83141 and x<0.93381:
        return 5.0 + 2.5 / 0.1024 * (x - 0.83141)
    elif x>=0.93381 and x<0.97232:
        return 7.5 + 2.5 / 0.03851 * (x - 0.93381)
    elif x>=0.97232 and x<0.99062:
        return 10.0 + 2.5 / 0.0183 * (x - 0.97232)
    elif x>=0.99062 and x<0.99724:
        return 12.5 + 2.5 / 0.00662 * (x - 0.99062)
    elif x>=0.99724 and x<0.99915:
        return 15.0 + 2.5 / 0.00191 * (x - 0.99724)
    elif x>=0.99915 and x<0.99974:
        return 17.5 + 2.5 / 0.00059 * (x - 0.99915)
    elif x>=0.99974 and x<0.99996:
        return 20.0 + 2.5 / 0.00022 * (x - 0.99974)
    elif x>=0.99996 and x<=1:
        return 22.5 + 2.5 / 0.00004 * (x - 0.99996)

# construct model platform
venice_segments = {}
venice_segments['length'] = np.array([6.0, 5.0, 4.0, 3.4, 1.6, 5.2, 10.4], order='F')
venice_segments['zhs'] = np.array([-4.0, -3.0, -2.0, -1.0, 0.0, 0.25, 0.5, 1.5], order='F')
venice_segments['pop'] = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 319.2, 319.2], order='F')
site_x, site_zh, site_pop = utils.construct_tai_platform(venice_segments)
xref = utils.get_refshore_coordinate(site_x, site_zh)
pft_grids = {}
pft_grids['x'] = np.array([0, 0.8, 2.4, 6.8, 17.2], order='F') + xref
pft_grids['pft'] = np.array([2, 2, 0, 0, 0], dtype=np.int32, order='F')
site_pft = utils.construct_platform_pft(pft_grids, site_x)
nx = len(site_x)

site_rslr = 0.0     # relative sea level rise (mm/yr)
site_rslr = site_rslr * 1e-3 / 8.64e4 / 365     # m/s
rhoSed = 2650.0     # sediment density (kg/m3)
rhoOM = 1200.0      # OM density (kg/m3)
porSed = 0.4        # porosity
d50 = 25e-6         # sediment median diameter (m)
Twav = 2.0          # wave period (s)
MHT = 0.75          # mean high water level (msl)
Ttide = 12.0        # tide period (hr)

npft = 9
Cz0 = 65.0
Kdf = 100.0
cbc = 0.015
fr = 0.78
alphaA = np.zeros(npft, dtype=np.float64, order='F')
betaA = np.zeros(npft, dtype=np.float64, order='F')
alphaD = np.zeros(npft, dtype=np.float64, order='F')
betaD = np.zeros(npft, dtype=np.float64, order='F')
cD0 = np.zeros(npft, dtype=np.float64, order='F')
ScD = np.zeros(npft, dtype=np.float64, order='F')
alphaA[2:6] = 8.0
betaA[2:6] = 0.5
alphaD[2:6] = 0.005
betaD[2:6] = 0.3
cD0[2:6] = 1.1
ScD[2:6] = -0.3

# driving data for eco-geomorphology model
# In an old version of ALBM, I have a method to calculate Tsoi from Tair

# temporal variables
npool = 2
nvar = 5
site_OM = np.zeros((nx,npool), dtype=np.float64, order='F')
site_Esed = np.zeros(nx, dtype=np.float64, order='F')
site_Dsed = np.zeros(nx, dtype=np.float64, order='F')
site_Bag = np.zeros(nx, dtype=np.float64, order='F')
site_Bbg = np.zeros(nx, dtype=np.float64, order='F')

odir = '/Users/tanz151/Downloads/'
site_id = 466

# instantiate hydrodynamics model
taihydro.inithydromod(site_x, site_zh, nvar, npft)
taihydro.setmodelparams(d50, Cz0, Kdf, cbc, fr, alphaA, betaA, 
                        alphaD, betaD, cD0, ScD)

# instantiate ecogeomorphology models
mac_params = {'d50': d50, 'rhoSed': rhoSed, 'porSed': porSed}
mac_module = importlib.import_module('minac_mod')
mac_class = getattr(mac_module, 'F06MOD')
mac_mod = mac_class(mac_params)

omac_params = {'rhoOM': rhoOM, 'phi': 2.2, 'aa': 15.5, 'bb': -18.55, 
               'cc': -1.364}
omac_module = importlib.import_module('omac_mod')
omac_class = getattr(omac_module, 'NULLMOD')
omac_mod = omac_class(omac_params)

models = {'taihydro': taihydro, 'mac_mod': mac_mod, 'omac_mod': omac_mod}

uhydro_tol = np.array([1e-3,1e-4,1e-6,1e-3,1e-3], dtype=np.float64, order='F')  # tolerance
dyncheck = np.array([1,0,1,1,1], dtype=np.int32, order='F') # check negative value

# first spinup for one year
try:
    date0 = date(2009,1,1)
    date1 = date(2009,2,1)
    # construct hydrodynamic model boundary conditions
    nday = (date1-date0).days
    nhour = 24 * nday
    h0 = MHT * np.sin(2*np.pi*np.arange(nhour)/Ttide)
    dh0 = np.gradient(h0)
    U0 = dh0 / np.max(np.abs(dh0)) * 1.0
    U10 = np.array([icdf_wind_gen(random.uniform(0,1)) for ii in range(nday)], 
                    dtype=np.float64)
    Css0 = 8.0e-3 * np.ones(nday, dtype=np.float64)
    Cj0 = 28.0 * np.ones(nday, dtype=np.float64)
    input_data = {'rk4_mode': 101, 'date0': date0, 'date1': date1, 'tstep': 'day',
                  'x': site_x, 'pft': site_pft, 'zh': site_zh, 'Bag': site_Bag, 
                  'Bbg': site_Bbg, 'SOM': site_OM, 'Esed': site_Esed, 
                  'Dsed': site_Dsed, 'MHT': MHT, 'U10': U10, 'Twav': Twav, 
                  'h0': h0, 'U0': U0, 'Css0': Css0, 'Cj0': Cj0, 
                  'rslr': site_rslr, 'tol': uhydro_tol, 'dyncheck': dyncheck}
    tai_state, uhydro_out, ecogeom_out = utils.run_tai_maces(input_data, \
        models, True, True)
except AssertionError as errstr:
    # print error message
    print("Model stops due to that", errstr)
finally:
    site_zh = tai_state['zh']
    site_pft = tai_state['pft']
    site_Bag = tai_state['Bag']
    site_Bbg = tai_state['Bbg']
    site_OM = tai_state['SOM']
    site_Esed = tai_state['Esed']
    site_Dsed = tai_state['Dsed']

# regular simulation
try:
    date0 = date(2009,1,1)
    date1 = date(2010,1,1)
    # construct hydrodynamic model boundary conditions
    nday = (date1-date0).days
    nhour = 24 * nday
    h0 = MHT * np.sin(2*np.pi*np.arange(nhour)/Ttide)
    dh0 = np.gradient(h0)
    U0 = dh0 / np.max(np.abs(dh0)) * 1.0
    U10 = np.array([icdf_wind_gen(random.uniform(0,1)) for ii in range(nday)], 
                    dtype=np.float64)
    Css0 = 8.0e-3 * np.ones(nday, dtype=np.float64)
    Cj0 = 28.0 * np.ones(nday, dtype=np.float64)
    input_data = {'rk4_mode': 101, 'date0': date0, 'date1': date1, 'tstep': 'day',
                  'x': site_x, 'pft': site_pft, 'zh': site_zh, 'Bag': site_Bag, 
                  'Bbg': site_Bbg, 'SOM': site_OM, 'Esed': site_Esed, 
                  'Dsed': site_Dsed, 'MHT': MHT, 'U10': U10, 'Twav': Twav, 
                  'h0': h0, 'U0': U0, 'Css0': Css0, 'Cj0': Cj0, 
                  'rslr': site_rslr, 'tol': uhydro_tol, 'dyncheck': dyncheck}
    tai_state, uhydro_out, ecogeom_out = utils.run_tai_maces(input_data, \
        models, False, True)
except AssertionError as errstr:
    # print error message
    print("Model stops due to that", errstr)
finally:
    # write simulation outputs
    utils.write_hydro_outputs(odir, site_id, uhydro_out)
    utils.write_ecogeom_outputs(odir, site_id, ecogeom_out)
    # deallocate
    taihydro.finalizehydromod()