
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.lines as mlines

t_steps = 4
x_steps = 3
i_steps = 3


def load_errors():
    dataType = np.float32
    data = np.fromfile("error_Npart.txt", dataType)

    size_x = t_steps * x_steps * i_steps + 2 * x_steps * t_steps
    size_y = 5
    data = np.reshape(data, (size_y, size_x))
    return data


def loadFile(filename, shape):
    dataType = np.float32
    data = np.fromfile(filename, dataType)

    data = np.reshape(data, (shape[0], shape[1]))
    return data


def plot_errors():
    data = load_errors()
    x = [0.1, 0.05, 0.01]
    dt = 0.2, 0.1, 0.05, 0.01
    interp = 'Bilinear', 'Bicubic', 'Bicubic spline'

    meanPlot = plt.figure()

    print("Error scales as $e \approx C * n^a$, where $n$ is number of timesteps.  \n")
    for dx in range(0, x_steps):
        plt.plot(dt, data[0, t_steps * dx:t_steps * dx + t_steps],
                 marker='*', label='Bilinear $\Delta x$ = ' + str(x[dx]))
        plt.plot(
            dt, data[0, t_steps * x_steps + t_steps *
                     dx:t_steps * dx + t_steps * x_steps + t_steps],
            marker='s', label='Bicubic $\Delta x$ = ' + str(x[dx])
        )
        plt.plot(
            dt, data[0, 2 * t_steps * x_steps + t_steps *
                     dx:t_steps * dx + t_steps + 2 * t_steps * x_steps],
            marker='o', label='Splines $\Delta x$ = ' + str(x[dx])
        )
        print(
            "splines,$\Delta t$ = 0.1 $\Delta x$= " + repr(x[dx]) +
            ", C = " +
            repr(data[4, 2 * t_steps * x_steps + t_steps * dx + 3])
            + ", a = " +
            repr(data[3, 2 * t_steps * x_steps + t_steps * dx + 3])
            + " \n"
        )
        print(
            "Average error at tmax: " +
            repr(data[2, 2 * t_steps * x_steps + t_steps * dx + 3]) + "\n"
        )

    for dx in range(0, x_steps):
        plt.plot(
            dt, data[0, i_steps * t_steps * x_steps +
                     t_steps * dx:t_steps * dx + t_steps +
                     i_steps * t_steps * x_steps],
            marker='x',  label='RK4 double steps $\Delta x$ = ' + str(x[dx])
        )

        plt.plot(
            dt, data[0, i_steps * t_steps * x_steps +
                     t_steps * x_steps + t_steps * dx:t_steps * dx +
                     t_steps + t_steps * x_steps + i_steps * t_steps * x_steps],
            marker='s', label='Heuns $\Delta x$ = ' + str(x[dx])
        )

        print(
            "RK4 - double step,$\Delta t$ = 0.1 $\Delta x$ = " + repr(x[dx]) +
            ", C = " +
            repr(data[4, i_steps * t_steps * x_steps +
                      t_steps * dx + 1])
            + ", a = " +
            repr(data[3, i_steps * t_steps * x_steps +
                      t_steps * dx + 1])
            + " \n"
        )
        print(
            "Average error at tmax: " + repr(data[2, i_steps * t_steps * x_steps +
                                                  t_steps * dx + 1]) + "\n"
        )
    plt.legend()
    tmax = 20
    rk4Shape1 = int(tmax / 0.05), 1, 1
    error_RK41 = loadFile('RK4_error1.txt', rk4Shape1)

    rk4Shape2 = int(tmax / 0.1), 1, 1
    error_RK42 = loadFile('RK4_error2.txt', rk4Shape2)

    rk4_2xShape = int(tmax / 0.2), 1, 1
    error_RK4_2xdt = loadFile('RK4_2xdt_error.txt', rk4_2xShape)

    heunsShape = int(tmax / 0.05), 1, 1
    error_heuns = loadFile('heuns_error.txt', heunsShape)

    f, errorplots = plt.subplots(1, 2, figsize=(15,8))
    errorplots[0].plot(np.linspace(0, tmax, rk4Shape1[0]),
                       error_RK41, label='RK4 lin interp in t, $\Delta t$ = 0.05')
    errorplots[0].plot(np.linspace(0, tmax, rk4Shape2[0]),
                       error_RK42, label='RK4 lin interp in t, $\Delta t$ = 0.1')
    errorplots[0].plot(np.linspace(0, tmax, rk4_2xShape[0]),
                       error_RK4_2xdt, label='RK4, $\Delta t$ = 0.2')
    errorplots[0].plot(np.linspace(0, tmax, heunsShape[0]),
                       error_heuns, label='heuns, $\Delta t$ = 0.05')
    errorplots[0].set_title('Comparison of different integrators')
    errorplots[0].set_xlabel('Time / s')
    errorplots[0].legend()
    errorplots[1].plot(np.linspace(0,rk4_2xShape[0], rk4_2xShape[0]), error_RK4_2xdt, label='RK4, $\Delta t$ = 0.2')
    approx = data[4, i_steps * t_steps * x_steps +
                  t_steps * 2 +1] * np.exp(np.linspace(0, rk4_2xShape[0], rk4_2xShape[0]) *
                  (data[3, i_steps * t_steps * x_steps +
                      t_steps * 2 + 1])  ) 
    

    errorplots[1].plot(np.linspace(0, rk4_2xShape[0], rk4_2xShape[0]), approx, 'k-.', label=r'$e(n)$ = $Ce^{n\alpha}$')
    errorplots[1].set_title('Error for RK4 with no interpolation in time')
    errorplots[1].legend()
    #errorplots[0].set_yscale('log')
    #errorplots[1].set_xscale('log')
    errorplots[1].set_xlabel('$n$ steps')
    errorplots[0].set_ylabel('Average error %')
    errorplots[1].set_ylabel('Average error %')
    #f.subplots_adjust(hspace=4.6)
    plt.tight_layout(pad=0.2)
    f.savefig('errorPlot_T20.png', dpi = 200)
    #plt.ylim(0, 1)
    plt.show()


plot_errors()
