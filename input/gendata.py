import numpy as np
import matplotlib
matplotlib.use('Agg')
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from shutil import copy
from os import mkdir
import otis_tide_pred as otp
import os
import logging
import shutil,os,glob



logging.basicConfig(level=logging.DEBUG)

_log = logging.getLogger(__name__)


def lininc(n,Dx,dx0):
  a=(Dx-n*dx0)*2./n/(n+1)
  dx = dx0+arange(1.,n+1.,1.)*a
  return dx



def copy_dirs(outdir):
    #### Set up the output directory
    backupmodel=1
    if backupmodel:
        try:
            mkdir(outdir0)
        except:
            import datetime
            import time
            ts = time.time()
            st=datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d%H%M%S')
            shutil.move(outdir0[:-1],outdir0[:-1]+'.bak'+st)
            mkdir(outdir0)

            _log.info(outdir0+' Exists')
        outdir=outdir0
        try:
            mkdir(outdir)
        except:
            _log.info(outdir+' Exists')
        outdir=outdir+'input/'
        mkdir(outdir)
        mkdir(outdir+'/figs/')
        copy('gendata.py', outdir)
        for tocp in ['data', 'data.pkg', 'data.obcs', 'data.kl10', 'README',
                    'eedata', 'otis_tide_pred.py', 'runModel.sh',
                    'data.diagnostics']:
            copy(tocp, outdir)

    ## Copy some other files
    _log.info( "Copying files")

    try:
      shutil.rmtree(outdir+'/../code/')
    except:
      _log.info("code is not there anyhow")
    shutil.copytree('../code', outdir+'/../code/')
    shutil.copytree('../python', outdir+'/../python/')

    try:
      shutil.rmtree(outdir+'/../build/')
    except:
      _log.info("build is not there anyhow")
    _log.info(outdir+'/../build/')
    mkdir(outdir+'/../build/')
    builddir = outdir+'/../build/'
    copy('../build/Makefile', builddir)
    copy('../build/mitgcmuv', builddir)

    # copy any data that is in the local indata
    shutil.copytree('../indata/', outdir+'/../indata/')


Hmax = 5200
runname = 'KC18r02'
comments = 'Kauaii Channel HighResRun'
outdir0 = '../results/' + runname + '/'
outdir = outdir0 + '/input/'
indir = outdir0 + '/indata/'
figdir = outdir0 + 'input/figs'

copy_dirs(outdir0)

# These must match ../code/SIZE.h
nx = 16 * 50  # 8*64
ny = 16 * 50 # 8*64
nz = 250
innerw = int(440 / 2)

#########
# dz

Hmax = 5200
dz = np.ones(250) * 10
for i in range(100, 250):
    dz[i] = dz[i-1] * 1.012
z = np.cumsum(dz)
Hmax = z[-1]

with open(outdir+"/delZvar.bin", "wb") as f:
	dz.tofile(f)
f.close()

if 1:
    fig, ax = plt.subplots()
    ax.plot(dz, z)
    fig.savefig(figdir + '/dz.pdf')

# FLIP!
lon0 = -158-37.77/60;
lat0 =21+40.780/60;
mpernm = 1852

#### Horizontal Grid

def ll2xy(lon, lat):
    x = (lon - lon0) * mpernm * 60 * np.cos(lat0 * np.pi / 180)
    y = (lat - lat0) * mpernm * 60
    return x, y


def xy2ll(x, y):
    lon = lon0 + x / (mpernm * 60 * np.cos(lat0 * np.pi / 180))
    lat = lat0 + y / (mpernm * 60)
    return lon, lat

def doubleinterp(x, y, z, xint, yint):
    """
    Interpolate z first in x and then in y.
    """
    m, n = np.shape(z)
    mnew = len(yint)
    nnew = len(xint)

    znew0 = np.zeros((m, nnew))
    for i in range(m):
        ind = np.isfinite(z[i, :])
        znew0[i, :] = np.interp(xint, x[ind], z[i, ind])

    znew = np.zeros((mnew, nnew))
    for j in range(nnew):
        ind = np.isfinite(znew0[:, j])
        znew[:, j] = np.interp(yint, y[ind], znew0[ind, j])

    return znew

# dx
expandrate = 1.025

dx = np.ones(nx)
inner = range(int(nx/2) - innerw, int(nx/2) + innerw)
dx[inner] = 100
for i in range(inner[-1], nx):
    dx[i] = dx[i-1] * expandrate
for i in range(inner[0], -1, -1):
    dx[i] = dx[i+1] * expandrate

x = np.cumsum(dx)
x = x - np.mean(x)

with open(outdir+"/delXvar.bin", "wb") as f:
	dx.tofile(f)
f.close()

# dy

dy = np.ones(ny)
inner = range(int(ny/2) - innerw, int(ny/2) + innerw)
dy[inner] = 100
for i in range(inner[-1], ny):
    dy[i] = dy[i-1] * expandrate
for i in range(inner[0], -1, -1):
    dy[i] = dy[i+1] * expandrate

y = np.cumsum(dy)
y = y - np.mean(y)

with open(outdir+"/delYvar.bin", "wb") as f:
	dy.tofile(f)
f.close()

# plot
if 1:
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(x/1000., dx)
    ax[1].plot(y/1000., dy)
    fig.savefig(figdir+'/dxdy.pdf')

############################
# Topography
ds = xr.open_dataset('../indata/OahuFilled.nc')
ds['x'], ds['y'] = ll2xy(ds['lon'], ds['lat'])

ss = xr.open_dataset('../indata/SmithSandwellNearHi.nc')
ss['x'], ss['y'] = ll2xy(ss['longitude'] - 360, ss['latitude'])

topo0 = doubleinterp(ss.x, ss.y, ss.ROSE, x, y)
iny = (ds.y.values[0] <= y) & ( y <= ds.y.values[-1])
inx = (ds.x.values[0] <= x) & (x <= ds.x.values[-1])
print(len(np.where(inx)[0]))
topo = topo0.copy()
toponew = doubleinterp(ds.x, ds.y, ds.depths, x[inx], y[iny])
for nn, i in enumerate(np.where(inx)[0]):
    topo[iny, i] = toponew[:, nn]

topo[topo < -Hmax] = -Hmax


with open(outdir+"/topo.bin", "wb") as f:
	topo.tofile(f)
f.close()

# plot
if 1:
    fig, ax = plt.subplots()
    ax.pcolormesh(x/1.e3, y/1e3, topo, rasterized=True, vmin=-6000, vmax=0)
    ind = np.where(np.diff(x) == 100)[0]
    dx = x[ind[-1]] - x[ind[0]]
    dy = y[ind[-1]] - y[ind[0]]
    rec = mpatches.Rectangle((x[ind[0]] / 1e3, y[ind[0]] / 1e3),
            dx / 1e3, dy / 1e3,
            facecolor='none', edgecolor='r')
    ax.add_artist(rec)
    ind = np.where(np.diff(x) <= 1000)[0]
    dx = x[ind[-1]] - x[ind[0]]
    dy = y[ind[-1]] - y[ind[0]]
    rec = mpatches.Rectangle((x[ind[0]] / 1e3, y[ind[0]] / 1e3),
            dx / 1e3, dy / 1e3,
            facecolor='none', edgecolor='y')
    ax.add_artist(rec)
    indx = np.arange(14, len(x)-14)
    indy = np.arange(14, len(y)-14)
    dx = x[indx[-1]] - x[indx[0]]
    dy = y[indy[-1]] - y[indy[0]]
    rec = mpatches.Rectangle((x[indx[0]] / 1e3, y[indy[0]] / 1e3),
            dx / 1e3, dy / 1e3,
            facecolor='none', edgecolor='c')
    ax.add_artist(rec)
    ax.set_aspect(1)
    fig.savefig(figdir+'/topo.pdf')


########
# Temperature profile...
tp = xr.open_dataset('../indata/TempProfile.nc')
z = np.cumsum(dz)
T0 = np.interp(z, tp.z.values, tp['T'].values)

with open(outdir+"/TRef.bin", "wb") as f:
	T0.tofile(f)
f.close()

# plot:
if 1:
    fig, ax = plt.subplots()
    ax.plot(T0,z)
    fig.savefig(figdir+'/TO.pdf')

#############
# U, V on boundaries from tidal predictions...
#

# hourly updates (3600 s)
dates = np.arange(np.datetime64('2002-08-28'), np.datetime64('2002-09-26'),
                  dtype='datetime64[h]')
print(f'{len(dates)} dates')

# South
lons, lats = xy2ll(x, y[0] + 0 *x)
h, u, v, depths = otp.tide_pred('../indata/OtisHOME/Model_haw', lons, lats,
                        dates)
bad = depths < 50; u[:, bad] = 0; v[:, bad] = 0

# must be time, z, x/y
u = u[:, np.newaxis, :] + 0 * z[np.newaxis, :, np.newaxis]
v = v[:, np.newaxis, :] + 0 * z[np.newaxis, :, np.newaxis]

with open(outdir+"/Us.bin","wb") as f:
  u.tofile(f)
with open(outdir+"/Vs.bin","wb") as f:
  v.tofile(f)

# North
lons, lats = xy2ll(x, y[-1] + 0 *x)
h, u, v, depths = otp.tide_pred('../indata/OtisHOME/Model_haw', lons, lats, dates)
bad = depths < 50; u[:, bad] = 0; v[:, bad] = 0

u = u[:, np.newaxis, :] + 0 * z[np.newaxis, :, np.newaxis]
v = v[:, np.newaxis, :] + 0 * z[np.newaxis, :, np.newaxis]

with open(outdir+"/Un.bin","wb") as f:
  u.tofile(f)
with open(outdir+"/Vn.bin","wb") as f:
  v.tofile(f)

# West
lons, lats = xy2ll(x[0] + 0 * y, y)
h, u, v, depths = otp.tide_pred('../indata/OtisHOME/Model_haw', lons, lats, dates)
bad = depths < 50; u[:, bad] = 0; v[:, bad] = 0

u = u[:, np.newaxis, :] + 0 * z[np.newaxis, :, np.newaxis]
v = v[:, np.newaxis, :] + 0 * z[np.newaxis, :, np.newaxis]

with open(outdir+"/Uw.bin","wb") as f:
  u.tofile(f)
with open(outdir+"/Vw.bin","wb") as f:
  v.tofile(f)

# East
lons, lats = xy2ll(x[-1] + 0 * y, y)
h, u, v, depths = otp.tide_pred('../indata/OtisHOME/Model_haw', lons, lats, dates)
bad = depths < 50; u[:, bad] = 0; v[:, bad] = 0


u = u[:, np.newaxis, :] + 0 * z[np.newaxis, :, np.newaxis]
v = v[:, np.newaxis, :] + 0 * z[np.newaxis, :, np.newaxis]

print(np.shape(u))

with open(outdir+"/Ue.bin","wb") as f:
  u.tofile(f)
with open(outdir+"/Ve.bin","wb") as f:
  v.tofile(f)



_log.info('Writing info to README')
############ Save to README
with open('README','r') as f:
  data=f.read()
with open('README','w') as f:
  import datetime
  import time
  ts = time.time()
  st=datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
  f.write( st+'\n')
  f.write( outdir+'\n')
  f.write(comments+'\n\n')
  f.write(data)

_log.info('All Done!')

_log.info('Archiving to home directory')

try:
    shutil.rmtree('../archive/'+runname)
except:
    pass

shutil.copytree(outdir0+'/input/', '../archive/'+runname+'/input')
shutil.copytree(outdir0+'/python/', '../archive/'+runname+'/python')
shutil.copytree(outdir0+'/code', '../archive/'+runname+'/code')
