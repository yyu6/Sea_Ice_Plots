from netCDF4 import Dataset
import numpy as np
easegrid = Dataset('easegrid_xy.nc')

x_grid = np.array(easegrid.variables['x_grid'])
y_grid = np.array(easegrid.variables['y_grid'])

PATH_U = 'anl_isentrop.033_ugrd.reg_tl319.201001_201012.nc'
PATH_V = 'anl_isentrop.034_vgrd.reg_tl319.201001_201012.nc'

wind_v = netCDF4.Dataset(PATH_V,'r')
wind_u = netCDF4.Dataset(PATH_U,'r')

u_lats = wind_u.variables['UGRD_GDS4_THEL_S123'][0,0,:,:].data
u_lons = wind_u.variables['UGRD_GDS4_THEL_S123'][0,0,:,:].data


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

