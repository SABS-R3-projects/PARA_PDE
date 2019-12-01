import numpy as np
import AnalyticalSolution
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def initial_conditions(x, a):
    return AnalyticalSolution.AnalyticalSolution(x, 0, alpha=a)

def fitzhugh_nagumo_solver(alpha=0.2, beta=1, gamma=1, max_x=1, max_t=10):
    """A solver for the fitzhugh-nagumo equation using the  accurate space derivative (ASD) method as laid out in
    Jeno Gazdag's and Jose Canos's paper: Numerical Solution of Fisher's Equation

        This is done by splitting our problem into 2 equations: Simple diffusion: df/dt = d^2u/dx^2 and the extra
        terms dg/dt = u*(1 - u)*(u - alpha). Where the change in u for our equation is the sum of the change in f and
        the change in g. The diffusion equation is then solved using fourier transforms and the other equation is
        solved using taylor expansion.

            :param max_x: upper bound of positional dimension x, this must be specified.
            :param max_t: amount of time to run the simulation for.
            :param alpha: A constant of the equation which should obey 0 < alpha < 0.5, the default value is 0.2.
            :param beta: A constant of the equation, the default value is 1.
            :param gamma: A constant of the equation, the default value is 1.
            :return: A matrix containing transmembrane potential values for each position x at a given t (each column is
            a given x and each row is a given t)

    Source of numerical solution method is below:
    @article{10.2307/3212689,
    ISSN = {00219002},
    URL = {http://www.jstor.org/stable/3212689},
    abstract = {The accurate space derivative (ASD) method for the numerical treatment of nonlinear partial differential equations has been applied to the solution of Fisher's equation, a nonlinear diffusion equation describing the rate of advance of a new advantageous gene, and which is also related to certain water waves and plasma shock waves described by the Korteweg-de-Vries-Burgers equation. The numerical experiments performed indicate how from a variety of initial conditions, (including a step function, and a wave with local perturbation) the concentration of advantageous gene evolves into the travelling wave of minimal speed. For an initial superspeed wave this evolution depends on the cutting off of the right-hand tail of the wave, which is physically plausible; this condition is necessary for the convergence of the ASD method. Detailed comparisons with an analytic solution for the travelling waves illustrate the striking accuracy of the ASD method for other than very small values of the concentration.},
    author = {Jenö Gazdag and José Canosa},
    journal = {Journal of Applied Probability},
    number = {3},
    pages = {445--457},
    publisher = {Applied Probability Trust},
    title = {Numerical Solution of Fisher's Equation},
    volume = {11},
    year = {1974}
    }

    """
    # Initial conditions: u(x,0) = I(x)
    #
    # Boundary conditions for the fourier transform: periodic in x, if x' = x + 2nL then u(x', t) = u(x, t) for all
    # integers n. To translate this boundary condition into our problem we need to leave some buffer space around the
    # area we are concerned with and make it symmetrical around the origin
    L = 2*max_x
    meshpoints = 1000
    x = np.linspace(-L, L, meshpoints)
    y = np.linspace(0, max_x, round(meshpoints/2))  # the area of x that we are actually concerned with
    numsteps = 32
    delta_t = max_t/numsteps

    initial = initial_conditions(10*x, alpha)
    initial = np.hstack((np.flip(initial), initial))
    u_current = initial

    u_relevant = u_current[-len(y):]  # the values that we are concerned with

    all_u = [u_relevant]  # for the collection of our results

    def delta_diffusion(u):
        fourier = np.fft.fft(u)
        delta_x = meshpoints / max_x
        j = np.linspace(-meshpoints, meshpoints, 2*meshpoints)
        k_j = j * np.pi / meshpoints * delta_x
        a = np.exp(-np.power(k_j, 2) * delta_t)
        delta_fourier = a * fourier
        delta_f = np.fft.ifft(delta_fourier)
        return beta * delta_f

    def delta_extra(u):
        dg_dt = u * (1 - u) * (u - alpha)
        calc_1 = -3*np.power(u, 2) + 2*(1+alpha)*u - alpha
        d2g_dt2 = dg_dt * calc_1
        calc_2 = -6*u + 2*(1+alpha)
        d3g_dt3 = d2g_dt2*calc_1 + dg_dt*calc_2
        delta_g = dg_dt*delta_t + d2g_dt2*np.power(delta_t, 2)/2 + d3g_dt3*np.power(delta_t, 3)/6
        return gamma * delta_g

    def delta_u(u):
        d_u = delta_diffusion(u) + delta_extra(u)
        return d_u

    for i in range(0, numsteps-1):
        u_new = u_current + delta_u(u_current)
        u_current = u_new
        u_relevant = u_current[-len(y):]
        all_u.append(np.real(u_relevant))


    all_u = np.array(all_u)

    return all_u, max_t, max_x


def animate(alpha=0.2, beta=1, gamma=1, max_x=1, max_t=10):
    """Produces an animation and graph for the solution of the fitzhugh-nagumo equation using the accurate space
    derivative (ASD) method to solve.

    :param max_x: upper bound of positional dimension x, this must be specified.
    :param max_t: amount of time to run the simulation for.
    :param alpha: A constant of the equation which should obey 0 < alpha < 0.5, the default value is 0.2.
    :param beta: A constant of the equation, the default value is 1.
    :param gamma: A constant of the equation, the default value is 1.
    :return: A matrix containing transmembrane potential values for each position x at a given t (each column is
            a given x and each row is a given t)
    """
    all_u, _, _ = fitzhugh_nagumo_solver(alpha, beta, gamma, max_x, max_t)
    y = np.linspace(0, max_x, 500)
    numsteps = 32
    delta_t = max_t/numsteps

    fig, ax = plt.subplots()
    ln, = plt.plot(y, all_u[0])
    ax.set_ylim(0, 1.1)

    def update(frame):
        ln.set_data(y, all_u[frame])
        ax.set_title('t = {}'.format(frame * delta_t))
        return ln,

    ani = FuncAnimation(fig, update, frames=numsteps, interval=30, blit=False, repeat=False)
    plt.show()

    plt.plot(y, all_u[-1], label='End')
    plt.plot(y, all_u[0], label='Start')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    animate()
