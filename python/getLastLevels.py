import xarray as xr
# 'Channel5k1000_01',
for td in [ 'KC18r01']:
    todo = f'../results/{td}/input/levels.nc'

    with xr.open_dataset(todo) as ds:
        print(ds)
        ds = ds.isel(record=-1)
        print(ds)
        for v in ['hFacC', 'hFacS', 'hFacW']:
            del ds[v]
        ds.to_netcdf(f'../reduceddata/LevelLast{td}.nc')
