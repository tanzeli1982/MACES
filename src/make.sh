f2py --fcompiler=gfortran --opt='-O3 -Wno-unused-variable' -L/usr/lib -I. -lblas -llapack -m TAIHydroMOD -c rungekutta4_mod.f90 tai_hydro_mod.f90 
