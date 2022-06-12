
import os
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import netCDF4
from os import listdir
from os.path import isfile, join
import numpy as np
import seaborn as sns
import cartopy
from numpy import cos, dtype, meshgrid, sin
from scipy.interpolate.ndgriddata import griddata
import imagesc as imagesc
from sympy import *
import h5py
from numpy import cos, dtype, meshgrid, sin
import cartopy.crs as ccrs
import seaborn as sns
import matplotlib.pyplot as plt

easegrid = Dataset('easegrid_xy.nc')

SEA_ICE_CONCENTRATION = '/Users/yuyaoning/Desktop/research/2012 ice conc/month'

SEA_VELOCITY_TEST = "/Users/yuyaoning/Desktop/research/Full_DOT_data_Arco.nc"
test_v = netCDF4.Dataset(SEA_VELOCITY_TEST, 'r')

SEA_ICE_VELOCITY = "ice velocity/icemotion_weekly_nh_25km_20030101_20031231_v4.1.nc"
sea_ice_v = netCDF4.Dataset(SEA_ICE_VELOCITY, 'r')

PATH_U1 = 'anl_isentrop.033_ugrd.reg_tl319.201001_201012.nc'
PATH_V1 = 'anl_isentrop.034_vgrd.reg_tl319.201001_201012.nc'

PATH_U2 = 'anl_isentrop.033_ugrd.reg_tl319.200901_200912.nc'
PATH_V2 = 'anl_isentrop.034_vgrd.reg_tl319.200901_200912.nc'

PATH_U3 = 'anl_isentrop.033_ugrd.reg_tl319.201101_201112.nc'
PATH_V3 = 'anl_isentrop.034_vgrd.reg_tl319.201101_201112.nc'

# transfer to lon/lats

def psn2ll(x, y):
    phi_c = 70
    a = 6378137.0
    e = 0.08181919
    lambda_0 = -45

    phi_c = phi_c*np.pi/180
    lambda_0 = lambda_0*np.pi/180

    t_c = np.tan(np.pi/4 - phi_c/2)/((1-e*np.sin(phi_c))/(1+e*np.sin(phi_c)))**(e/2)
    m_c = np.cos(phi_c)/np.sqrt(1-e**2*(np.sin(phi_c))**2)
    rho = np.sqrt(x**2+y**2)
    t = rho*t_c/(a*m_c)

    chi = np.pi/2 - 2 * np.arctan(t)
    phi = chi+(e**2/2 + 5*e**4/24 + e**6/12 + 13*e**8/360)*np.sin(2*chi) + (7*e**4/48 + 29*e**6/240 + 811*e**8/11520)*np.sin(4*chi) + (7*e**6/120+81*e**8/1120)*np.sin(6*chi) + (4279*e**8/161280)*np.sin(8*chi)
    lambda_1 = lambda_0 + np.arctan2(x,-y)

    lambda_1=np.mod(lambda_1+np.pi,2*np.pi)-np.pi

    lat=phi*180/np.pi
    lon=lambda_1*180/np.pi

    return [lat, lon]


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

def get_u_v(ssh_):
    f = 1.4e-4
    g = 9.8
    dssh_dy, dssh_dx = np.gradient(ssh_)
    U = (-g*dssh_dy) / f
    V = (g*dssh_dx) / f
    return U, V

x_grid = np.array(easegrid.variables['x_grid'])
y_grid = np.array(easegrid.variables['y_grid'])

wind_v1 = netCDF4.Dataset(PATH_V1,'r')
wind_u1 = netCDF4.Dataset(PATH_U1,'r')

wind_v2 = netCDF4.Dataset(PATH_V2,'r')
wind_u2 = netCDF4.Dataset(PATH_U2,'r')

wind_v3 = netCDF4.Dataset(PATH_V3,'r')
wind_u3 = netCDF4.Dataset(PATH_U3,'r')


xx, yy = np.meshgrid(x_grid, y_grid)
def ease_ssh(lat, lon, value):
    cols = np.zeros(lat.shape)
    rows = np.zeros(lon.shape)
    for j in range(0,320):
        cols[j], rows[j] = easegrid(11, lat[j], lon[j])
    col = cols.ravel()
    row = rows.ravel()
    value_ = np.array(value)
    value_ = value_.ravel()
    # value_[indices] = 0
    #return LinearNDInterpolator((col, row), ssh_)
    return griddata((col, row), value_, (xx, yy), method="linear")

def same_dimension(wind):
    new_wind = np.zeros((320,320))
    for i in range(len(wind)):
        new_wind[i] = wind[i][::2]
    return new_wind


def plot(x, y, U, V, N=8):
    fig = plt.figure(figsize=[10, 10])
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='gray', facecolor=sns.xkcd_rgb['tan'])
    extent = 2500000
    x = x[::N]
    y = y[::N]
    U = U[::N,::N]
    V = V[::N,::N]
    ax.set_extent((-extent, extent, -extent, extent), crs=ccrs.NorthPolarStereo())
    ax.quiver(x, y, U, V, transform=ccrs.NorthPolarStereo(), scale=100)
    plt.show()



lats1 = np.array(wind_u1.variables['g4_lat_2'])
lons1 = np.array(wind_u1.variables['g4_lon_3'])

lats2 = np.array(wind_u2.variables['g4_lat_2'])
lons2 = np.array(wind_u2.variables['g4_lon_3'])

lats3 = np.array(wind_u3.variables['g4_lat_2'])
lons3 = np.array(wind_u3.variables['g4_lon_3'])

wind_v3 = np.array(wind_v3.variables['VGRD_GDS4_THEL_S123'][0,0])
wind_u3 = np.array(wind_u3.variables['UGRD_GDS4_THEL_S123'][0,0])

wind_v3 = same_dimension(wind_v3)[0]
wind_u3 = same_dimension(wind_u3)[0]


transferv = ease_ssh(lats3, lons3[::2], wind_v3)
transferu = ease_ssh(lats3, lons3[::2], wind_u3)


x = np.array(sea_ice_v.variables["x"])
y = np.array(sea_ice_v.variables["y"])

plot(x, y, transferu, transferv)