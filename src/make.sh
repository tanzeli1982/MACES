#!/bin/bash

arg=$( echo $1 | tr '[:upper:]' '[:lower:]' )
if [ -z "$arg" ]; then
   gfortran -O3 -mmacosx-version-min=10.9 -c -fPIC data_buffer_mod.f90 hydro_utilities_mod.f90
   f2py -c --quiet --fcompiler=gnu95 --opt='-O3' -I. data_buffer_mod.o hydro_utilities_mod.o -L/usr/lib -lblas -llapack -m TAIHydroMOD tai_hydro_mod.f90
elif [ $arg = 'clean' ]; then
   sources=$( ls *.f90 )
   objects=$( find . -type f \( -name \$sources -o -name \*.o \) )
   modules=$( find . -type f \( -name \$sources -o -name \*.mod \) )
   rm -f TAIHydroMOD.*.so $objects $modules
elif [ $arg = 'debug' ]; then
   gfortran -g -fbacktrace -fcheck=bounds -mmacosx-version-min=10.9 -c -fPIC data_buffer_mod.f90 hydro_utilities_mod.f90 
   #f2py --debug-capi -c --quiet --fcompiler=gnu95 --opt='-g -fbacktrace -fcheck=bounds' -I. data_buffer_mod.o hydro_utilities_mod.o -L/usr/lib -lblas -llapack -m TAIHydroMOD tai_hydro_mod.f90
   f2py -c --quiet --fcompiler=gnu95 --opt='-g -fbacktrace -fcheck=bounds' -I. data_buffer_mod.o hydro_utilities_mod.o -L/usr/lib -lblas -llapack -m TAIHydroMOD tai_hydro_mod.f90
else
   echo "Wrong Argument: $1!!!"
fi
