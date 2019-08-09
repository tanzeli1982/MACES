import numpy as np
import maces_utilities as utils
import random
import minac_mod as mac
import omac_mod as omac
from TAIHydroMOD import taihydro
from TAIHydroMOD import rungekutta4 as rk4

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
venice_segments['length'] = np.array([6.0, 5.0, 4.0, 3.4, 1.6, 5.2, 10.4])
venice_segments['zhs'] = np.array([-4.0, -3.0, -2.0, -1.0, 0.0, 0.25, 0.5, 1.5])
venice_segments['pop'] = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 319.2, 319.2])
site_x, site_zh, site_pop = utils.construct_tai_platform(venice_segments)
xref = utils.get_refshore_coordinate(venice_segments)
pft_grids = {}
pft_grids['x'] = np.array([0, 0.8, 2.4, 6.8, 17.2]) + xref
pft_grids['pft'] = np.array([2, 2, 0, 0, 0], dtype=np.int32)
site_pft = utils.construct_platform_pft(pft_grids, site_x)
Nx = len(site_x)

# initialize soil organic matter pools
site_OM = np.zeros((2,Nx), dtype=np.float64)

site_rslr = 0.0     # relative sea level rise (mm/yr)
site_rslr = 1e-3 / 8.64e4 / 365 * site_rslr     # m/s
rhoSed = 2650.0     # sediment density (kg/m3)
rhoOM = 1200.0      # OM density (kg/m3)
porSed = 0.4        # porosity
hydro_params = {}
hydro_params['d50'] = 1e-4
hydro_params['Cz0'] = 65.0
hydro_params['Kdf'] = 100.0
hydro_params['cbc'] = 0.015
hydro_params['fr'] = 0.78
hydro_params['alphaA'] = 8.0 * np.ones(8, dtype=np.float64)
hydro_params['betaA'] = 0.5 * np.ones(8, dtype=np.float64)
hydro_params['alphaD'] = 0.005 * np.ones(8, dtype=np.float64)
hydro_params['betaD'] = 0.3 * np.ones(8, dtype=np.float64)
hydro_params['cD0'] = 1.1 * np.ones(8, dtype=np.float64)
hydro_params['ScD'] = -0.3 * np.ones(8, dtype=np.float64)

# eco-geomorphology models
mac_params = {}
mac_mod = mac.M12MOD(mac_params)

omac_params = {}
omac_mod = omac.M12MOD(omac_params)

# hydrodynamic model initialization and set up
taihydro.Initialize(site_x, site_zh)

taihydro.SetModelParameters(hydro_params['d50'], hydro_params['Cz0'], 
                            hydro_params['Kdf'], hydro_params['cbc'], 
                            hydro_params['fr'], hydro_params['alphaA'], 
                            hydro_params['betaA'], hydro_params['alphaD'], 
                            hydro_params['betaD'], hydro_params['cD0'], 
                            hydro_params['ScD'])

# construct hydrodynamic model boundary conditions
nyear = 10
nday = nyear * 365
nhour = 24 * nday
Ttide = 12.0
h0 = 0.75 * np.sin(2*np.pi*np.arange(nhour)/Ttide)
dh0 = np.gradient(h0)
U0 = dh0 / np.max(np.abs(dh0)) * 1.0
U10 = np.array([icdf_wind_gen(random.uniform(0,1)) for ii in range(nhour)], 
                dtype=np.float64)
Twav = 2.0
Hwav_ocean = 0.27 * U10^2 / utils.G
Hwav0 = np.sinh(utils.Karman*h0) / np.sinh(utils.Karman*30.0) * Hwav_ocean
Css0 = 8.0 * np.ones(nhour, dtype=np.float64)
Cj0 = 28.0 * np.ones(nhour, dtype=np.float64)

# driving data for eco-geomorphology model
# in an old version of ALBM, there is a method to link Tair and Tsoi

# temporal variables
site_Esed = np.zeros(Nx, dtype=np.float64)
site_Dsed = np.zeros(Nx, dtype=np.float64)
site_Bag = np.zeros(Nx, dtype=np.float64)
site_DepOM = np.zeros(Nx, dtype=np.float64)
tmp_uhydro = np.zeros((5,Nx), dtype=np.float64)
uhydro_out = {}
uhydro_out['h'] = 1e20 * np.ones((nhour,Nx), dtype=np.float32)
uhydro_out['U'] = 1e20 * np.ones((nhour,Nx), dtype=np.float32)
uhydro_out['Hwav'] = 1e20 * np.ones((nhour,Nx), dtype=np.float32)
uhydro_out['tau'] = 1e20 * np.ones((nhour,Nx), dtype=np.float32)
uhydro_out['Css'] = 1e20 * np.ones((nhour,Nx), dtype=np.float32)
uhydro_out['Cj'] = 1e20 * np.ones((nhour,Nx), dtype=np.float32)
ecogeom_out = {}
ecogeom_out['zh'] = 1e20 * np.ones((nday,Nx), dtype=np.float32)
ecogeom_out['Esed'] = 1e20 * np.ones((nday,Nx), dtype=np.float32)
ecogeom_out['Dsed'] = 1e20 * np.ones((nday,Nx), dtype=np.float32)
ecogeom_out['DepOM'] = 1e20 * np.ones((nday,Nx), dtype=np.float32)
ecogeom_out['Bag'] = 1e20 * np.ones((nday,Nx), dtype=np.float32)
ecogeom_out['Bbg'] = 1e20 * np.ones((nday,Nx), dtype=np.float32)
ecogeom_out['SOM'] = 1e20 * np.ones((nday,2,Nx), dtype=np.float32)

odir = '/qfs/projects/taim/TAIMOD/test/'
site_id = 466

# run simulation
MAX_OF_STEP = 1800  # maximum simulation time step (s)
rk4_mode = 101      # adaptive mode
uhydro_tol = np.array([1e-6,1e-6,1e-6,1e-6,1e-6,1e-6], dtype=np.float64)  # toleratance
t = 0
tf = 8.64e4 * nday
hindx = -1
dindx = -1
ncount = 0
isHourNode = False
isDayNode = False
curstep = np.array([50], dtype=np.float64)
nextstep = np.array([MAX_OF_STEP], dtype=np.float64)
try:
    while t < tf:
        if t>=3.6e3*hindx and hindx<nhour:
            print('time step', int(hindx))
            isHourNode = True
            hindx = hindx + 1
            if np.mod(hindx,24)==0:
                isDayNode = True
                dindx = dindx + 1
        # simulate hydrodynamics
        taihydro.ModelSetup(site_zh, site_pft, site_Bag, site_Esed, site_Dsed, 
                            Twav, U10[hindx], h0[hindx], U0[hindx], 
                            Hwav0[hindx], Css0[hindx], Cj0[hindx])
        error = np.array([0], dtype=np.int32)
        rk4.RK4Fehlberg(taihydro.TAIHydroEquations, taihydro.m_uhydro, 
                        tmp_uhydro, rk4_mode, uhydro_tol, curstep, nextstep, 
                        error)
        assert error[0]==0, "Runge-Kutta iteration is more than MAXITER"
        taihydro.ModelCallback(tmp_uhydro)
        # simulate eco-geomorphology
        mac_inputs = {}
        site_Esed = mac_mod.mineral_suspend(mac_inputs)
        site_Dsed = mac_mod.mineral_deposition(mac_inputs)
        site_Lbed = mac_mod.bed_loading(mac_inputs)
        omac_inputs = {}
        site_Bag = omac_mod.aboveground_biomass(omac_inputs)
        omac_inputs['Bag'] = site_Bag
        site_Bbg = omac_mod.belowground_biomass(omac_inputs)
        site_DepOM = omac_mod.organic_deposition(omac_inputs)
        site_DecayOM = omac_mod.soilcarbon_decay(omac_inputs)
        # update soil OM pool
        DepOM_pools = np.zeros((2,Nx), dtype=np.float64)
        DepOM_pools[0] = 0.158 * site_DepOM
        DepOM_pools[1] = 0.842 * site_DepOM
        site_OM = site_OM + (DepOM_pools - site_DecayOM) * curstep
        # update platform elevation
        site_zh = site_zh + ((site_Dsed/rhoSed + site_Lbed/rhoSed + \
            site_DepOM/rhoOM - site_Esed/rhoSed)/(1.0-porSed) - \
            site_rslr) * curstep
        # archive hydrodynamic state variables
        if isHourNode:
            uhydro_out['h'][hindx] = taihydro.m_uhydro[0,:]
            uhydro_out['U'][hindx] = taihydro.m_U
            uhydro_out['Hwav'][hindx] = taihydro.m_uhydro[2,:]
            uhydro_out['tau'][hindx] = taihydro.m_tau
            uhydro_out['css'][hindx] = taihydro.m_uhydro[3,:]
            uhydro_out['sal'][hindx] = taihydro.m_uhydro[4,:]
        if isDayNode:
            # archive daily mean eco-geomorphology variables
            ecogeom_out['zh'] = site_zh
            ecogeom_out['Esed'][dindx] = site_Esed
            ecogeom_out['Dsed'][dindx] = site_Dsed
            ecogeom_out['DepOM'][dindx] = site_DepOM
            ecogeom_out['Bag'][dindx] = site_Bag
            ecogeom_out['Bbg'][dindx] = site_Bbg
            ecogeom_out['SOM'][dindx] = site_OM
        # check small time step
        isHourNode = False
        isDayNode = False
        if curstep[0]<0.1:
            ncount = ncount + 1
            err_msg = 'Run diverge at step' + '{:d}'.format(hindx)
            assert ncount<=100, err_msg
            nextstep = 50.0
        else:
            ncount = 0
        t = t + curstep[0]
        curstep = nextstep
        nextstep = MAX_OF_STEP
        
except AssertionError as errstr:
    # print error message
    print("Model stops due to that", errstr)
finally:
    # write simulation outputs
    utils.write_outputs(odir, site_id, uhydro_out, ecogeom_out)
    # deallocate
    taihydro.Destruct()
