import numpy as np
import maces_utilities as utils
import random
import minac_mod as mac
import omac_mod as omac
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
site_rslr = 1e-3 / 8.64e4 / 365 * site_rslr     # m/s
rhoSed = 2650.0     # sediment density (kg/m3)
rhoOM = 1200.0      # OM density (kg/m3)
porSed = 0.4        # porosity
npft = 9
hydro_params = {}
hydro_params['d50'] = 25e-6
hydro_params['Cz0'] = 65.0
hydro_params['Kdf'] = 100.0
hydro_params['cbc'] = 0.015
hydro_params['fr'] = 0.78
hydro_params['alphaA'] = np.zeros(npft, dtype=np.float64, order='F')
hydro_params['betaA'] = np.zeros(npft, dtype=np.float64, order='F')
hydro_params['alphaD'] = np.zeros(npft, dtype=np.float64, order='F')
hydro_params['betaD'] = np.zeros(npft, dtype=np.float64, order='F')
hydro_params['cD0'] = np.zeros(npft, dtype=np.float64, order='F')
hydro_params['ScD'] = np.zeros(npft, dtype=np.float64, order='F')
hydro_params['alphaA'][2:6] = 8.0
hydro_params['betaA'][2:6] = 0.5
hydro_params['alphaD'][2:6] = 0.005
hydro_params['betaD'][2:6] = 0.3
hydro_params['cD0'][2:6] = 1.1
hydro_params['ScD'][2:6] = -0.3

# eco-geomorphology models
npool = 2
site_OM = np.zeros((nx,npool), dtype=np.float64, order='F')

mac_params = {}
mac_mod = mac.F06MOD(mac_params)

omac_params = {'phi': 2.2, 'aa': 15.5, 'bb': -18.55, 'cc': -1.364}
omac_mod = omac.NULLMOD(omac_params)

# hydrodynamic model initialization and set up
nvar = 5
taihydro.inithydromod(site_x, site_zh, nvar, npft)
taihydro.setmodelparams(hydro_params['d50'], hydro_params['Cz0'], \
    hydro_params['Kdf'], hydro_params['cbc'], hydro_params['fr'], \
    hydro_params['alphaA'], hydro_params['betaA'], hydro_params['alphaD'], \
    hydro_params['betaD'], hydro_params['cD0'], hydro_params['ScD'])

# construct hydrodynamic model boundary conditions
nyear = 10
nday = nyear * 365
nhour = 24 * nday
Ttide = 12.0
h0 = 0.75 * np.sin(2*np.pi*np.arange(nhour)/Ttide)
dh0 = np.gradient(h0)
U0 = dh0 / np.max(np.abs(dh0)) * 1.0
U10 = np.array([icdf_wind_gen(random.uniform(0,1)) for ii in range(nday)], 
                dtype=np.float64)
Twav = 2.0
Css0 = 8.0 * np.ones(nday, dtype=np.float64)
Cj0 = 28.0 * np.ones(nday, dtype=np.float64)

# driving data for eco-geomorphology model
# In an old version of ALBM, I have a method to calculate Tsoi from Tair

# temporal variables
site_Esed = np.zeros(nx, dtype=np.float64, order='F')
site_Dsed = np.zeros(nx, dtype=np.float64, order='F')
site_Bag = np.zeros(nx, dtype=np.float64, order='F')
site_DepOM = np.zeros(nx, dtype=np.float64, order='F')
sim_h = np.zeros(nx, dtype=np.float64, order='F')
sim_U = np.zeros(nx, dtype=np.float64, order='F')
sim_Hwav = np.zeros(nx, dtype=np.float64, order='F')
sim_tau = np.zeros(nx, dtype=np.float64, order='F')
sim_Css = np.zeros(nx, dtype=np.float64, order='F')
sim_Cj = np.zeros(nx, dtype=np.float64, order='F')
uhydro_out = {}
uhydro_out['h'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
uhydro_out['U'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
uhydro_out['Hwav'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
uhydro_out['tau'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
uhydro_out['Css'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
uhydro_out['Cj'] = 1e20 * np.ones((nhour,nx), dtype=np.float32)
ecogeom_out = {}
ecogeom_out['pft'] = -1 * np.ones((nday,nx), dtype=np.int8)
ecogeom_out['zh'] = 1e20 * np.ones((nday,nx), dtype=np.float32)
ecogeom_out['Esed'] = 1e20 * np.ones((nday,nx), dtype=np.float32)
ecogeom_out['Dsed'] = 1e20 * np.ones((nday,nx), dtype=np.float32)
ecogeom_out['DepOM'] = 1e20 * np.ones((nday,nx), dtype=np.float32)
ecogeom_out['Bag'] = 1e20 * np.ones((nday,nx), dtype=np.float32)
ecogeom_out['Bbg'] = 1e20 * np.ones((nday,nx), dtype=np.float32)
ecogeom_out['SOM'] = 1e20 * np.ones((nday,nx,npool), dtype=np.float32)

print('simulation starts')

odir = '/Users/tanz151/Downloads/'
site_id = 466

# run simulation
MAX_OF_STEP = 1800  # maximum simulation time step (s)
rk4_mode = 101      # adaptive mode
uhydro_tol = np.array([1e-3,1e-4,1e-6,1e-3,1e-3], dtype=np.float64, order='F')  # tolerance
dyncheck = np.array([1,0,1,1,1], dtype=np.int32, order='F') # check negative value
assert len(uhydro_tol)==nvar, "size of uhydro_tol is not equal to # of variables"
t = 0
tf = 8.64e4 * nday
hindx = -1
dindx = -1
ncount = 0
isHourNode = False
isDayNode = False
curstep = 50.0
nextstep = MAX_OF_STEP
try:
    while t < tf:
        if t>=3.6e3*(hindx+1) and hindx+1<nhour:
            isHourNode = True
            hindx = hindx + 1
            print('time step', int(hindx))
            if np.mod(hindx,24)==0:
                isDayNode = True
                dindx = dindx + 1
        # simulate hydrodynamics
        if isHourNode:
            h0_abs = -site_zh[0] + h0[hindx]
            U10_hr = U10[dindx] + (hindx/24-dindx)*(U10[dindx+1]-U10[dindx])
            Hwav0 = utils.estimate_Hwav_seaward(U10_hr, h0_abs)
            #np.set_printoptions(precision=3, suppress=True)
            np.set_printoptions(precision=3, suppress=False)
            #print(np.array([U10[dindx],h0_abs,U0[hindx],Hwav0,Css0[dindx],Cj0[dindx]]))
            #print(np.array([h0_abs,U0[hindx],Hwav0]))
        taihydro.modelsetup(site_zh, site_pft, site_Bag, site_Esed, site_Dsed, 
                            Twav, U10_hr, h0_abs, U0[hindx], Hwav0, 
                            Css0[dindx], Cj0[dindx])
        curstep, nextstep, error = taihydro.modelrun(rk4_mode, uhydro_tol, 
                                                     dyncheck, curstep)
        assert error==0, "runge-Kutta iteration is more than MAXITER"
        taihydro.modelcallback()
        sim_h, sim_U, sim_Hwav, sim_tau, sim_Css, sim_Cj = taihydro.getmodelsims(nx)
        assert np.all(np.isfinite(sim_h)), "NaN h found"
        assert np.all(np.isfinite(sim_U)), "NaN U found"
        assert np.all(np.isfinite(sim_Hwav)), "NaN Hwav found"
        assert np.all(np.isfinite(sim_tau)), "NaN tau found"
        assert np.all(np.isfinite(sim_Css)), "NaN Css found"
        assert np.all(np.isfinite(sim_Cj)), "NaN Cj found"
        # get wet area length
        if isHourNode:
            wetL = site_x[sim_h>1e-6][-1]
            tmp_zh = site_zh - h0[hindx]
            wetL_potential = site_x[tmp_zh<0][-1]
            #print(np.array([wetL_potential, wetL, Hwav0, np.max(sim_Hwav)]))
            print(np.array([h0_abs, U0[hindx], sim_h[1], sim_U[1]]))
        # simulate eco-geomorphology
        mac_inputs = {'x': site_x, 'Css': sim_Css, 'tau': sim_tau, 
                      'd50': hydro_params['d50'], 'Rous': rhoSed}
        site_Esed = mac_mod.mineral_suspension(mac_inputs)
        site_Dsed = mac_mod.mineral_deposition(mac_inputs)
        site_Lbed = mac_mod.bed_loading(mac_inputs)
        omac_inputs = {'x': site_x, 'zh': site_zh, 'MHT': 0.75, 'pft': site_pft, 
                       'SOM': site_OM}
        site_Bag = omac_mod.aboveground_biomass(omac_inputs)
        omac_inputs['Bag'] = site_Bag
        site_Bbg = omac_mod.belowground_biomass(omac_inputs)
        site_DepOM = omac_mod.organic_deposition(omac_inputs)
        site_DecayOM = omac_mod.soilcarbon_decay(omac_inputs)
        # update soil OM pool
        DepOM_pools = np.zeros((nx,npool), dtype=np.float64, order='F')
        DepOM_pools[:,0] = 0.158 * site_DepOM
        DepOM_pools[:,1] = 0.842 * site_DepOM
        site_OM = site_OM + (DepOM_pools - site_DecayOM) * curstep
        # update platform elevation
        site_Lbed[:] = 0.0
        site_DepOM[:] = 0.0
        site_zh = site_zh + ((site_Dsed/rhoSed + site_Lbed/rhoSed + \
            site_DepOM/rhoOM - site_Esed/rhoSed)/(1.0-porSed) - \
            site_rslr) * curstep
        # archive hydrodynamic state variables
        if isHourNode:
            uhydro_out['h'][hindx] = sim_h
            uhydro_out['U'][hindx] = sim_U
            uhydro_out['Hwav'][hindx] = sim_Hwav
            uhydro_out['tau'][hindx] = sim_tau
            uhydro_out['Css'][hindx] = sim_Css
            uhydro_out['Cj'][hindx] = sim_Cj
        if isDayNode:
            # archive daily mean eco-geomorphology variables
            ecogeom_out['zh'][dindx] = site_zh
            ecogeom_out['Esed'][dindx] = site_Esed
            ecogeom_out['Dsed'][dindx] = site_Dsed
            ecogeom_out['DepOM'][dindx] = site_DepOM
            ecogeom_out['Bag'][dindx] = site_Bag
            ecogeom_out['Bbg'][dindx] = site_Bbg
            ecogeom_out['SOM'][dindx] = site_OM
            ecogeom_out['pft'][dindx] = site_pft
        # check small time step
        isHourNode = False
        isDayNode = False
        if curstep<0.1:
            ncount = ncount + 1
            err_msg = 'run diverge at step ' + '{:d}'.format(hindx)
            assert ncount<=100, err_msg
            nextstep = 50.0
        else:
            ncount = 0
        t = t + curstep
        curstep = nextstep
        nextstep = MAX_OF_STEP
        
except AssertionError as errstr:
    # print error message
    print("Model stops due to that", errstr)
finally:
    # write simulation outputs
    #utils.write_outputs(odir, site_id, uhydro_out, ecogeom_out)
    # deallocate
    taihydro.finalizehydromod()
    print('simulation ends')
