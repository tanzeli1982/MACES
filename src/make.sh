gfortran -O3 -mmacosx-version-min=10.9 -c -fPIC data_buffer_mod.f90 hydro_utilities_mod.f90
f2py -c --no-lower --fcompiler=gnu95 --opt='-O3' -I. data_buffer_mod.o hydro_utilities_mod.o -L/usr/lib -lblas -llapack -m TAIHydroMOD tai_hydro_mod.f90
