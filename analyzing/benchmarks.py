

import numpy as np
from matplotlib import pyplot as plt


SET_2 = "502x330x43"
SET_1 = '502x330x667'
SET = SET_1, SET_2
DELTAS = 0.8, 0.2, 0.08
PARTICLES = 18480, 51520, 115920, 1854720, 11592000
num_runs = 5
####
GTX960 = np.zeros((num_runs * 2, 2))
GTX960[0, 0] = 3.748
GTX960[0, 1] = 7.911
GTX960[1, 0] = 1.729
GTX960[1, 1] = 6.172
GTX960[2, 0] = 2.704    #translation step
GTX960[2, 1] = 7.115   #whole program
GTX960[3, 0] = 29.784
GTX960[3, 1] = 34.4
GTX960[4, 0] = 180.41
GTX960[4, 1] = 186.71
###shorts
GTX960[5, 0] = 0.765
GTX960[5, 1] = 2.705
GTX960[6, 0] = 0.801 ##
GTX960[6, 1] = 2.74 ###
GTX960[7, 0] = 0.866 ###
GTX960[7, 1] = 2.721###
GTX960[8, 0] = 2.586
GTX960[8, 1] = 4.576
GTX960[9, 0] = 12.261
GTX960[9, 1] = 15.951

####
####
GTX750TI = np.zeros((num_runs * 2, 2))

GTX750TI[0, 0] = 2.02     #translation step
GTX750TI[0, 1] = 2.02+0.9+0.9     #whole program
GTX750TI[1, 0] = 3.04
GTX750TI[1, 1] = 4.84
GTX750TI[2, 0] = 5.13
GTX750TI[2, 1] = 6.94
GTX750TI[3, 0] = 60.41
GTX750TI[3, 1] = 64.74
GTX750TI[4, 0] = 0
GTX750TI[4, 1] = 0
###shorts
GTX750TI[5, 0] = 1.37
GTX750TI[5, 1] = 7.85
GTX750TI[6, 0] = 1.27 ##
GTX750TI[6, 1] = 2.83 ##
GTX750TI[7, 0] = 1.31 ##
GTX750TI[7, 1] = 2.85 ##
GTX750TI[8, 0] = 2.56 ##
GTX750TI[8, 1] = 6.71 ##
GTX750TI[9, 0] = 9.58 ##
GTX750TI[9, 1] = 30.12 ##
####
####
I7_5930K = np.zeros((num_runs * 2, 2))

I7_5930K[0, 0] = 7.118 ##
I7_5930K[0, 1] = 7.665 ##
I7_5930K[1, 0] = 18.22##
I7_5930K[1, 1] = 18.72 ##
I7_5930K[2, 0] = 39.8     #translation step
I7_5930K[2, 1] = 40.337      #whole program
I7_5930K[3, 0] = 637.987
I7_5930K[3, 1] = 639.42
I7_5930K[4, 0] = 3955.13 ##
I7_5930K[4, 1] = 3963.12 ##
##shorss
I7_5930K[5, 0] = 0.576
I7_5930K[5, 1] = 0.62
I7_5930K[6, 0] = 1.402 ##
I7_5930K[6, 1] = 1.467 ##
I7_5930K[7, 0] = 3.121 ##
I7_5930K[7, 1] = 3.219 ##
I7_5930K[8, 0] = 49.5
I7_5930K[8, 1] = 50.454
I7_5930K[9, 0] = 308.442
I7_5930K[9, 1] = 314.158

####
def plot_times():

    fig = plt.figure()
    plt_n = 1
    for i in range(0,1):
        plt.plot(PARTICLES, (I7_5930K[num_runs*i:num_runs*(i + 1), plt_n]), label=SET[i] + ' - I7_5930K', marker='*')
        plt.plot(PARTICLES,  (GTX750TI[num_runs*i:num_runs*(i + 1), plt_n]), label=SET[i] + ' - GTX750TI', marker='o')
        plt.plot(PARTICLES,  (GTX960[num_runs*i:num_runs*(i + 1), plt_n]), label=SET[i] + ' - GTX960', marker='s')
        plt.legend()
    fig2 = plt.figure()

    for i in range(1,2):
        plt.plot(PARTICLES, (I7_5930K[num_runs*i:num_runs*(i + 1), plt_n]), label=SET[i] + ' - I7_5930K', marker='*')
        plt.plot(PARTICLES,  (GTX750TI[num_runs*i:num_runs*(i + 1), plt_n]), label=SET[i] + ' - GTX750TI', marker='o')
        plt.plot(PARTICLES,  (GTX960[num_runs*i:num_runs*(i + 1), plt_n]), label=SET[i] + ' - GTX960', marker='s')
        plt.legend()
    
    plt.show()


plot_times()

    