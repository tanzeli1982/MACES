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

# eco-geomorphology models
mac_params = {}
mac_mod = mac.M12MOD(mac_params)

omac_params = {}
omac_mod = omac.M12MOD(omac_params)

# hydrodynamic model initialization and set up
taihydro.Initialize(site_x, site_zh)

taihydro.SetModelParameters(d50, )

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

# create numpy array as Fortran-contiguous ordered
rk4.k1 = np.zeros((4,20), dtype=np.float64, order='F')
rk4.k2 = np.zeros((4,20), dtype=np.float64, order='F')
rk4.k3 = np.zeros((4,20), dtype=np.float64, order='F')
rk4.k4 = np.zeros((4,20), dtype=np.float64, order='F')
rk4.k5 = np.zeros((4,20), dtype=np.float64, order='F')
rk4.k6 = np.zeros((4,20), dtype=np.float64, order='F')
rk4.nxt4th = np.zeros((4,20), dtype=np.float64, order='F')
rk4.nxt5th = np.zeros((4,20), dtype=np.float64, order='F')
rk4.interim = np.zeros((4,20), dtype=np.float64, order='F')
rk4.rerr = np.zeros((4,20), dtype=np.float64, order='F')
invars = np.zeros((4,20), dtype=np.float64, order='F')
outvars = np.zeros((4,20), dtype=np.float64, order='F')
invars[0,:] = 1
invars[1,:] = 1
invars[2,:] = 1e-3
invars[3,:] = 1e-1
curstep = np.array([1.0], dtype=np.float64, order='F')
nextstep = np.array([1.0], dtype=np.float64, order='F')
mode = 102
outerr = np.array([0], dtype=np.int32, order='F')
tol = np.array([1e-2, 1e-2, 1e-5, 1e-4], dtype=np.float64, order='F')

rk4.rk4fehlberg(odeFunc, invars, outvars, mode, tol, curstep, nextstep, outerr)

# release allocated memory
taihydro.Destruct()
