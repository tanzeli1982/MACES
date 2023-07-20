#!/bin/bash

module purge
module load intel netcdf/4.7.4 mvapich2

export PATH=/share/apps/python/anaconda3-2020.02/bin:$PATH
source activate /qfs/people/tanz151/.conda/envs/work_env

arg=$( echo $1 | tr '[:upper:]' '[:lower:]' )
if [ -z "$arg" ]; then
   #ifort -O3 -fast -c -fPIC data_buffer_mod.f90 hydro_utilities_mod.f90
   #f2py -c --quiet --fcompiler=intelem --opt='-O3 -fast' -I. data_buffer_mod.o hydro_utilities_mod.o -m TAIHydroMOD tai_hydro_mod.f90
   ifort -O3 -c -fPIC data_buffer_mod.f90 hydro_utilities_mod.f90 
   f2py -c --quiet --fcompiler=intelem --opt='-O3' -I. data_buffer_mod.o hydro_utilities_mod.o -m TAIHydroMOD tai_hydro_mod.f90 
elif [ $arg = 'clean' ]; then
   sources=$( ls *.f90 )
   objects=$( find . -type f \( -name \$sources -o -name \*.o \) )
   modules=$( find . -type f \( -name \$sources -o -name \*.mod \) )
   rm -f TAIHydroMOD.*.so $objects $modules
elif [ $arg = 'debug' ]; then
   ifort -CB -g -traceback -fpe0 -c -fPIC data_buffer_mod.f90 hydro_utilities_mod.f90 
   #f2py --debug-capi -c --quiet --fcompiler=intelem --opt='-CB -g -traceback -fpe0' -I. data_buffer_mod.o hydro_utilities_mod.o -m TAIHydroMOD tai_hydro_mod.f90
   f2py -c --quiet --fcompiler=intelem --opt='-CB -g -traceback -fpe0' -I. data_buffer_mod.o hydro_utilities_mod.o -m TAIHydroMOD tai_hydro_mod.f90
else
   echo "Wrong Argument: $1!!!"
fi

conda deactivate
