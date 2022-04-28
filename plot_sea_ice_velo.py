import netCDF4
from os import listdir
from os.path import isfile, join
import numpy as np
import cartopy.crs as ccrs
import seaborn as sns
import cartopy
import matplotlib.pyplot as plt
SEA_ICE_VELOCITY = "data/sea_ice_velocity/icemotion_weekly_nh_25km_20200101_20201231_v4.1.nc"
sea_ice_v = netCDF4.Dataset(SEA_ICE_VELOCITY,'r')

def plot(x, y, U, V, N=4):
    fig = plt.figure(figsize=[10, 5])
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='gray', facecolor=sns.xkcd_rgb['tan'])
    extent = 2500000
    x = x[::N]
    y = y[::N]
    U = U[::N,::N]
    V = V[::N,::N]
    ax.set_extent((-extent,extent,-extent,extent), crs=ccrs.NorthPolarStereo())
    ax.quiver(x, y, U, V, transform=ccrs.NorthPolarStereo(), scale=100)
    plt.show()


def easegrid(iopt, alat, alon):
    R=6371.228
    C=25.06752
    ascale = 2* R/C
    pi2 = np.pi / 2.0
    dtr = pi2 / 90.0
    if iopt == 11:
        thelon = ascale * np.sin(alon * dtr) *  np.sin(dtr * (45.0 - 0.5 * alat))
        thelat = ascale * np.cos(alon * dtr) *  np.sin(dtr * (45.0 - 0.5 * alat))
    elif iopt == 12:  
        thelon = ascale * np.sin(alon * dtr) * np.cos(dtr * (45.0 - 0.5 * alat))
        thelat = ascale * np.cos(alon * dtr) * np.cos(dtr * (45.0 - 0.5 * alat))
    elif iopt == 13:  
        thelon = ascale * pi2 * alon * np.cos(30.0 * dtr) / 90.0
        thelat = ascale * np.sin(alat * dtr) / np.cos(30.0 * dtr)
    return thelat, thelon

u = sea_ice_v.variables["u"][1,:,:]
v = sea_ice_v.variables["v"][1,:,:]
x = np.array(sea_ice_v.variables["x"])
y = np.array(sea_ice_v.variables["y"])



plot(x, y, u, v)