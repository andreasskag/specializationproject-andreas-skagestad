
# coding: utf-8

# In[ ]:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D as ax




def plot_ftle(nx, ny, delta, dt):
        typ = np.float32
        filename_forw =  'output/'+str(nx) + '_' + str(ny) + '_ftle_heat_t5_delta' + delta + \
            '_dt' + dt + 'FORW.txt'
        filename_back =  'output/'+str(nx) + '_' + str(ny) + '_ftle_heat_t5_delta' + delta + \
            '_dt' + dt + 'BACK.txt'

        filename_particles =  'output/'+str(nx) + '_' + str(ny) + '_particle_ends.txt'

        # shape = np.genfromtxt(filename, skip_header=1, max_rows=1)
        shape = [ny, nx]

        data_forw = np.fromfile(filename_forw, dtype=typ)  # np.genfromtxt(filename)#
        data_forw = np.reshape(data_forw, (int(shape[0]), int(shape[1])))

        data_back = np.fromfile(filename_back, dtype=typ)  # np.genfromtxt(filename)#
        data_back = np.reshape(data_back, (int(shape[0]), int(shape[1])))

        particle_data = np.fromfile(filename_particles, dtype = typ)
        particle_data = np.reshape(particle_data, (2, (ny + 2) * (nx + 2)))

        particle_box = np.zeros((2, int((ny + 2) * (nx + 2) / 100)))

        x_start = int((nx + 2) / 2 - (nx + 2) / 20)
        y_start = int((ny + 2) / 2 - (ny + 2) / 20)
        for y in range(0, int((ny + 2) / 10)):
            for x in range(0,int((nx + 2) / 10)):
                particle_box[0,x + y * int(((nx + 2) / 10))] = \
                    particle_data[0, x_start + x + (y + y_start) * (nx + 2)]
                particle_box[1,x + y * int(((nx + 2) / 10))] = \
                    particle_data[1, x_start + x + (y + y_start) * (nx + 2)]

        data_forw[data_forw < 0] = 0.0
        data_back[data_back < 0] = 0.0
        data = data_forw*1- 1*data_back# 
        vmax = np.amax(np.abs(data))
        vmin = -vmax

        xdata = np.fromfile('output/'+str(nx) + '_' + str(ny) + 'grid_xcords.txt', dtype=typ)
        xdata = xdata[1:len(xdata) - 1]
        ydata = np.fromfile('output/'+str(nx) + '_' + str(ny) + 'grid_ycords.txt', dtype=typ)
        ydata = ydata[1:len(ydata) - 1]

        maskfile = 'output/mask.txt'

       #mask = np.fromfile(maskfile, dtype=int)
        #xmask = np.fromfile('ocean_data/x_coords_km.txt', dtype=type)
        #ymask = np.fromfile('ocean_data/y_coords_km.txt', dtype=type)
        #mask = np.ndarray.reshape(mask, (len(ymask), len(xmask)))
        #X,Y = np.meshgrid(xmask, ymask)

        #plt.contourf(X, Y, mask, cmap='gist_earth')
        fig = plt.figure()
        ax2 = fig.add_subplot(111, aspect='equal')
        mesh = plt.pcolormesh(xdata, ydata, data, cmap = plt.get_cmap('seismic'), vmin = vmin, vmax = vmax)
        #plt.colorbar(mesh)
        plt.autoscale(False)
        plt.scatter(particle_box[0,:], particle_box[1,:],c='k', s=0.4)
        
        bounding_x_start = xdata[x_start - 1]
        bounding_y_start = ydata[y_start - 1]
        bounding_x_len = xdata[x_start + int((nx + 2) / 10) - 2] - bounding_x_start
        bounding_y_len = ydata[y_start + int((ny + 2) / 10) - 2] - bounding_y_start
        
        ax2.add_patch(
            patches.Rectangle(
                (bounding_x_start, bounding_y_start),
                        bounding_x_len,
                         bounding_y_len,
                    fill=False      # remove background
                )
            )
        plt.show()

#plot_ftle(199, 99, '0.010000', '0.100000')
plot_ftle(1678, 1102, '0.200000', '7200.000000')
#plot_ftle(1607, 1076, '0.250000', '3600.000000')
#plot_ftle(1642, 792, '0.400000', '3600.000000he')
#plotTrajectory('0.100000');





def plot_error_ftle(nx, ny, delta1, dt1, dt_an):
        filename = 'output/'+str(nx) + '_' + str(ny) + '_ftle_heat_t5_delta' + delta1 + '_dt' + dt1 + '.txt'
        data = np.fromfile(filename, dtype=np.float32)  # np.genfromtxt(filename)#
        shape = [ny, nx]
        data = np.reshape(data, (int(shape[0]), int(shape[1])))

        filename_a ='output/'+ str(nx) + '_' + str(ny) + '_ftle_heat_t5_delta' + delta1 + '_dt' + dt_an + '.txt'
        data_a = np.fromfile(filename_a, dtype=np.float32)  # np.genfromtxt(filename)#
        shape_a = shape
        data_a = np.reshape(data_a, (int(shape[0]), int(shape[1])))
        difference = abs(data - data_a)
        error = difference / (abs(data_a) + 1 * 1e-15)
        print(np.min(abs(data_a[np.nonzero(data_a)])))
        xc = np.linspace(0, 2, int(shape[1]))
        yc = np.linspace(0, 1, int(shape[0]))
        sum_diff = np.sum(difference)
        sum_data_a = np.sum(abs(data_a))
        print(sum_diff / sum_data_a * 100)
        fig = plt.figure()
        plt.pcolormesh(xc, yc, data, cmap='GnBu', label=dt1)
        plt.title(dt1)
        plt.colorbar()
        fig = plt.figure()
        plt.pcolormesh(xc, yc, data_a, cmap='GnBu', label=dt_an)
        plt.title(dt_an)
        plt.colorbar()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X, Y = np.meshgrid(xc, yc)
        Z = error
        ax.plot_surface(X, Y, Z)
        #ax.set_zlim(bottom=-0, top=0.5, emit=True, auto=False)
        fig = plt.figure()
        square = 1;
        plt.show()
        plt.pcolormesh(xc[square:len(xc) - square], yc[square:len(yc) - square], error[square:len(yc) - square, square:len(xc) - square])
        plt.colorbar()
        plt.show()

##OLD FUNCTIONS
def plot_deform(nx, ny, delta, dt):
    filename = 'output/'+ str(nx) + '_' + str(ny) +  '_deform_heat_t5_delta' + delta + '_dt' + dt + '.txt'
    #shape = np.genfromtxt(filename, skip_header=1, max_rows=1)
    data = np.fromfile(filename, dtype=np.float32)#np.genfromtxt(filename)#
    shape = [ny, nx]
    data = np.reshape(data, (int(shape[0]), int(shape[1])))
    fig = plt.figure()
    xc = np.linspace(0, 2, int(shape[1]))
    yc = np.linspace(0, 1, int(shape[0]))
    plt.pcolormesh(xc, yc, data, cmap='coolwarm')
    plt.colorbar()
    #plt.contourf(data.T)
    plt.show()

def plotField() :
    datax = np.genfromtxt('output/testField_x_dx_0.010000_dt_0.100000.txt');
    datay = np.genfromtxt('output/testField_y_dx_0.010000_dt_0.100000.txt');
    dim = datax.shape;
    x = np.linspace(0,dim[0])
    y = np.linspace(0, dim[1])
    fig, ax = plt.subplots();
    #ax.streamplot(X, Y, datax, datay, color='r');
    plt.quiver(datax,datay)


    plt.show()

##
def plot_error_deform(nx, ny, delta1, dt1, dt_an):
    filename = 'output/'+str(nx) + '_' + str(ny) +  '_deform_heat_t5_delta' + delta1 + '_dt' + dt1 + '.txt'
    data = np.fromfile(filename, dtype=np.float32)  # np.genfromtxt(filename)#
    shape = [ny, nx]
    data = np.reshape(data, (int(shape[0]), int(shape[1])))


    filename_a = 'output/'+str(nx) + '_' + str(ny) +  '_deform_heat_t5_delta' + delta1 + '_dt' + dt_an + '.txt'
    data_a = np.fromfile(filename_a, dtype=np.float32)  # np.genfromtxt(filename)#
    shape_a = shape
    data_a = np.reshape(data_a, (int(shape[0]), int(shape[1])))
    error = np.divide(data , data_a)

    fig = plt.figure()
    xc = np.linspace(0, 2, int(shape[1]))
    yc = np.linspace(0, 1, int(shape[0]))

    fig = plt.figure()
    plt.pcolormesh(xc, yc, data, cmap='coolwarm')
    plt.colorbar()
    fig = plt.figure()
    plt.pcolormesh(xc, yc, data_a, cmap='coolwarm')
    plt.colorbar()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X,Y = np.meshgrid(xc,yc)
    Z = error
    ax.plot_surface(X,Y,Z)

    fig = plt.figure()
    plt.contourf(error)
    plt.show()

def plotTrajectory(dt):
    filename = 'output/trajectory_dt_' + dt + '_.txt'
    shape = np.genfromtxt(filename, skip_header=1, max_rows=1)
    data = np.genfromtxt(filename, skip_header=2)
    data = np.reshape(data ,(int(shape[1]),int(shape[0])))
    dim = data.shape
    #fig, ax = plt.subplots();
    start = 0 #int(2 * shape[2] / 4)

    county = 0
    countx = 0
   # for i in range(0, 8):
    fig = plt.plot(data[0,:],data[1,:])
        #fig = plt.plot(data[start + countx + county * 20, 0:int(shape[0] - 1), 0], data[start + countx + county * 20, 0:int(shape[0] - 1), 1])
        #countx = countx + 1
        #if (countx + 1 > 4):
        #    county = county + 1
        #    countx = 0

    plt.show()
def test_load_ocean():
    file = 'ocean_load_test.txt'
    x = 1103
    y = 653
    data = np.fromfile(file, dtype=np.int)
    data = data.reshape(y, x)
    x = np.fromfile('ocean_data/x_coords_km.txt')
    y = np.fromfile('ocean_data/y_coords_km.txt')
    fig = plt.figure()
    X,Y = np.meshgrid(x,y)
    plt.contourf(X,Y,data)
    plt.show()