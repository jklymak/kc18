## Complie directory

To compile execute
```
module load cray-netcdf-hdf5parallel

./MITgcm66h/tools/genmake2 -optfile=../build_options/conrad -mods=../code/ -rootdir=../MITgcm66h -mpi

make depend
make
```
