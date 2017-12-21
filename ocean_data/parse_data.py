 
import numpy as np
import xarray as xa

filename = 'Norkyst-800m-25d.nc'
#filename = 'Norkyst-800m-small.nc'
data = xa.open_dataset(filename)
type = np.float32
u = data.u
v = data.v
#u = list(u.drop('depth'))
u = xa.DataArray(u)
#u.where(u.data == 'nan', set = 0)
u = u.fillna(0)
u_array = u.data
u = None
tmin = 30
tmax = 150

#u_array = u_array[tmin:tmax,:,:,:]
shape = u_array.shape
u_array = np.ndarray.astype(u_array, type)
t = str(shape[0])
x = str(shape[3])
y = str(shape[2])
file_out = x + '_' + y + '_' + t + '_Norkyst_u.txt'
u_array.tofile(file_out)
u_array = None
v = xa.DataArray(v)
v = v.fillna(0)
v_array = v.data
#v_array = v_array[tmin:tmax,:,:,:]
shape = v_array.shape
v_array = np.ndarray.astype(v_array, type)
t = str(shape[0])
x = str(shape[3])
y = str(shape[2])
file_out = x + '_' + y + '_' + t + '_Norkyst_v.txt'
v_array.tofile(file_out)
a = 1

X = data.X
x = xa.DataArray(X)
x = x.data
x = np.ndarray.astype(x, type)
Y = data.Y
y = xa.DataArray(Y)
y= y.data
y = np.ndarray.astype(y, type)

x.tofile('x_coords_km.txt')
y.tofile('y_coords_km.txt')
