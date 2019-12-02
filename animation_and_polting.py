import matplotlib.pyplot as plt
import time
from matplotlib import animation

    
def animate(U, timeSteps: int, postionSteps: int, timeStepSize: float):

    fig= plt.figure()
    ims = []
    for i in range(timeSteps):
        im = plt.plot(postionSteps, U[:,i] , animated = True, color = 'red')
        ims.append(im)
    ani = animation.ArtistAnimation(fig, ims, interval = (10), blit = True, repeat_delay = 500)
    plt.show()

#animation(u, trial.k_N, trial.x_range, trial.k)