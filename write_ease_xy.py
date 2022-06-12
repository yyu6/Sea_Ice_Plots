from netCDF4 import Dataset    # Note: python is case-sensitive!
import numpy as np

PATH = "ice velocity/icemotion_weekly_nh_25km_20100101_20101231_v4.1.nc"

sea_ice_v = Dataset(PATH,'r')
x_ori = np.array(sea_ice_v.variables["x"])
y_ori = np.array(sea_ice_v.variables["y"])

ncfile = Dataset('easegrid_xy.nc',mode='w',format='NETCDF4_CLASSIC') 
x_dim = ncfile.createDimension('x_grid', 361)     # latitude axis
y_dim = ncfile.createDimension('y_grid', 361)

x = ncfile.createVariable('x_grid', np.float32, ('x_grid',))
x.long_name = 'x_grid'

x[:] = x_ori

y = ncfile.createVariable('y_grid', np.float32, ('y_grid',))
y.long_name = 'y_grid'

y[:] = y_ori