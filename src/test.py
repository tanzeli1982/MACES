import numpy as np
import maces_utilities as util
from TAIHydroMOD import taihydro
from TAIHydroMOD import rungekutta4 as rk4

# https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2012JF002363
# Test case: Venice Lagoon

venice_segments = []
util.construct_tai_platform(venice_segments)

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
print(outvars)
print(outerr)
print(nextstep)

rk4.k1 = None
rk4.k2 = None
rk4.k3 = None
rk4.k4 = None 
rk4.k5 = None 
rk4.k6 = None 
rk4.nxt4th = None 
rk4.nxt5th = None 
rk4.interim = None 
