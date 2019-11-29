import numpy as np
import AnalyticalSolution
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def initial_conditions(x, a):
    return AnalyticalSolution.AnalyticalSolution(x, 0, alpha=a)


def fitzhugh_nagumo_solver(alpha=0.2, gamma=1, delta=1, max_x=1, max_t=10):
    # Initial conditions: u(x,0) = I(x)
    # Boundry conditions u(0, t) = B = u(L, t)
    # Fourier Conditions: if x' = x + 2nL then u(x', t) = u(x, t) for all integers n
    L = max_x
    meshpoints = 160
    x = np.linspace(-L, L, meshpoints)
    y = np.linspace(0, L, round(meshpoints/4))
    numsteps = 10
    delta_t = max_t/numsteps

    initial = initial_conditions(10*x, alpha)
    # nought = np.zeros(2)
    nought = []
    initial = np.hstack((np.flip(initial), nought, initial))
    x = np.linspace(0, L, meshpoints)
    x = np.hstack((np.flip(-x), x))
    u_current = initial

    u_relevant = u_current[-len(y):]

    fig, ax = plt.subplots()
    ln, = plt.plot(y, u_relevant)
    ax.set_ylim(0, 1.1)

    all_u = [initial]

    def delta_diffusion(u):
        fourier = np.fft.fft(u)
        delta_x = meshpoints / max_x
        delta_a = []
        j = np.linspace(-len(x)/2, len(x)/2, len(x))
        k_j = j * np.pi / (len(x)/2) * delta_x
        a = np.exp(-np.power(k_j, 2) * delta_t)
        delta_fourier = a * fourier
        delta_f = np.fft.ifft(delta_fourier)
        print(delta_f)
        return delta_f

    def delta_extra(u):
        #print(u)
        dg_dt = u * (1 - u) * (u - alpha)
        # print(dg_dt)
        calc_1 = -3*np.power(u, 2) + 2*(1+alpha)*u - alpha
        #print(delta_g)
        d2g_dt2 = dg_dt * calc_1
        #print(delta_g)
        calc_2 = -6*u + 2*(1+alpha)
        #print(delta_g)
        d3g_dt3 = d2g_dt2*calc_1 + dg_dt*calc_2
        #print(delta_g)
        delta_g = dg_dt*delta_t + d2g_dt2*np.power(delta_t, 2)/2 + d3g_dt3*np.power(delta_t, 3)/6
        #print(delta_g)
        return delta_g

    def delta_u(u):
        d_u = delta_diffusion(u) + delta_extra(u)
        # d_u = delta_diffusion(u)
        # print(d_u)
        return d_u

    def update(frame):

        for i in range(1):
            u_new = u_current + delta_u(u_current)
            u_current[:] = u_new
            #u_current[0], u_current[-1] = 1, 1
            u_relevant[:] = u_current[-len(y):]
            all_u.append(u_relevant)
            print(u_current)

        ln.set_data(y, u_relevant)
        ax.set_title('t = {}'.format(10 * frame * delta_t))
        # print(frame)
        return ln,

    ani = FuncAnimation(fig, update, frames=max_t, interval=30, blit=False, repeat=False)
    plt.show()

    plt.plot(x, u_current)
    plt.plot(x, initial)
    plt.show()

fitzhugh_nagumo_solver()
#print(initial_conditions(-1, 0.2))