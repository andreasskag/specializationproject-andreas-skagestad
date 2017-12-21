       
    import numpy as np
    from matplotlib import pyplot as plt
    
        ###
        typ = np.float32
        filename_forw =  'output/'+str(nx) + '_' + str(ny) + '_ftle_heat_t5_delta' + delta + \
            '_dt' + dt + 'FORW.txt'
        filename_back =  'output/'+str(nx) + '_' + str(ny) + '_ftle_heat_t5_delta' + delta + \
            '_dt' + dt + 'BACK.txt'

        

        # shape = np.genfromtxt(filename, skip_header=1, max_rows=1)
        shape = [ny, nx]

        data_forw = np.fromfile(filename_forw, dtype=typ)  # np.genfromtxt(filename)#
        data_forw = np.reshape(data_forw, (int(shape[0]), int(shape[1])))

        data_back = np.fromfile(filename_back, dtype=typ)  # np.genfromtxt(filename)#
        data_back = np.reshape(data_back, (int(shape[0]), int(shape[1])))
        
        

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
        
        ###        
        filename_particles =  'output/'+str(nx) + '_' + str(ny) + '_particle_ends.txt'
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