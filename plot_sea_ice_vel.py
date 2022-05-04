import netCDF4
import h5py
from os import listdir
from os.path import isfile, join
import numpy as np
from numpy import cos, dtype, meshgrid, sin
import cartopy.crs as ccrs
import seaborn as sns
import cartopy
import matplotlib.pyplot as plt


OCEAN_VELOCITY_FOLDER = "data/ocean_velo_test/"

SEA_ICE_CONCENTRATION_FOLDER = "data/sea_ice_concentration"


SEA_VELOCITY_TEST = "data/ocean_velo_test/Full_DOT_data_Arco.nc"
test_v = netCDF4.Dataset(SEA_VELOCITY_TEST,'r')

SEA_ICE_VELOCITY = "data/sea_ice_velocity/icemotion_weekly_nh_25km_20200101_20201231_v4.1.nc"
sea_ice_v = netCDF4.Dataset(SEA_ICE_VELOCITY,'r')

print(test_v)

def easegrid(iopt, alat, alon):
    """EASE grid transformation
    (thelon thelat)=easegrid(iopt,lat,lon,ascale)

    computes the forward "ease" grid transform

    given a lat,lon (alat,alon) and the scale (ascale) the image
    transformation coordinates (thelon,thelat) are comuted
    using the "ease grid" (version 1.0) transformation given in fortran
    source code supplied by nsidc.

    the radius of the earth used in this projection is imbedded into
    ascale while the pixel dimension in km is imbedded in bscale
    the base values are: radius earth= 6371.228 km
                 pixel dimen =25.067525 km
    then, bscale = base_pixel_dimen
          ascale = radius_earth/base_pixel_dimen

    iopt is ease type: iopt=11=north, iopt=12=south, iopt=13=cylindrical
    """
    # ported from easegrid.m by JPB 21 Sept 2011
    R=6371.228
    C=25.06752
    ascale = 2*R/C
    pi2 = np.pi / 2.0
    dtr = pi2 / 90.0

    if iopt == 11:  # ease grid north
        thelon = ascale * sin(alon * dtr) * sin(dtr * (45.0 - 0.5 * alat))
        thelat = ascale * cos(alon * dtr) * sin(dtr * (45.0 - 0.5 * alat))
    elif iopt == 12:  # ease grid south
        thelon = ascale * sin(alon * dtr) * cos(dtr * (45.0 - 0.5 * alat))
        thelat = ascale * cos(alon * dtr) * cos(dtr * (45.0 - 0.5 * alat))
    elif iopt == 13:  # ease cylindrical
        thelon = ascale * pi2 * alon * cos(30.0 * dtr) / 90.0
        thelat = ascale * sin(alat * dtr) / cos(30.0 * dtr)

    cols = thelon * 25159
    rows = -thelat * 25159
    
    return cols, rows

x = np.array(sea_ice_v.variables["x"])
y = np.array(sea_ice_v.variables["y"])
xx, yy = np.meshgrid(x, y)

def ease_ssh(lat, lon, ssha):
    cols = np.zeros(lat.shape)
    rows = np.zeros(lon.shape)
    for j in range(lat.shape[1]):
        cols[:,j], rows[:,j] = easegrid(11, lat[:,j], lon[:,j])
    col = cols.ravel()
    row = rows.ravel()
    ssh_ = ssha.reshape(-1)
    indices = np.argwhere(ssh_>1e3)
    ssh_[indices] = 0
    #return LinearNDInterpolator((col, row), ssh_)
    return griddata((col, row), ssh_, (xx, yy), method="linear")

def get_u_v(ssh_):
    f = 1.4e-4
    g = 9.8
    dssh_dy, dssh_dx = np.gradient(ssh_)
    U = (-g*dssh_dy) / f
    V = (g*dssh_dx) / f
    return U, V

def plot(U, V, N=4, scale=30000):
    fig = plt.figure(figsize=[10, 5])
    ax = plt.axes(projection=ccrs.NorthPolarStereo())

    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='gray', facecolor=sns.xkcd_rgb['tan'])
    extent = 2500000
    #plt.contour(x, y, mask_land)

    _x = x[::N]
    _y = y[::N]
    _U = U[::N, ::N]
    _V = V[::N, ::N]

    ax.set_extent((-extent,extent,-extent,extent), crs=ccrs.NorthPolarStereo())
    ax.quiver(_x, _y, _U, _V, transform=ccrs.NorthPolarStereo(), scale=scale)
    plt.show()

lat_grid = test_v.variables["lats"][:]
lon_grid = test_v.variables["lons"][:]

dot_sum = np.mean(test_v.variables["DOT_smoothed"][60:100,:,:], axis=0)
dot_win = np.mean(test_v.variables["DOT_smoothed"][:60,:,:], axis=0)

ssh_grid_sum = ease_ssh(lat_grid, lon_grid, dot_sum)
ssh_grid_win = ease_ssh(lat_grid, lon_grid, dot_win)

u_sum, v_sum = get_u_v(ssh_grid_sum)
u_win, v_win = get_u_v(ssh_grid_win)

plot(u_sum, v_sum)
plot(u_win, v_win)