import matplotlib.pyplot as plt
import time
from matplotlib import animation

    
def animation(U, timeSteps: int, postionSteps: int, timeStepSize: float):

    fig= plt.figure()
    ims = []
    for i in range(timeSteps):
        im = plt.plot(postionSteps, U[:,i] , animated = True, color = 'red')
        ims.append(im)
    ani = animation.ArtistAnimation(fig, ims, interval = (100*timeStepSize), blit = True)
    plt.show()