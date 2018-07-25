import xarray as xr
import subprocess
# 'Channel5k1000_01',
for td in [ 'KC18r02']:
    todo = f'../results/{td}/input/levels.nc'

    with xr.open_dataset(todo) as ds:
        ds = ds.isel(record=range(-20, -1))
        for v in ['hFacC', 'hFacS', 'hFacW']:
            del ds[v]
        ds.to_netcdf(f'../reduceddata/LevelLast{td}.nc')
subprocess.run(['rsync', '-av', '../reduceddata/','valdez.seos.uvic.ca:kc18/reduceddata'])
