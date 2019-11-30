import numpy as np
import matplotlib.pyplot as plt
import treecode as tree
import time

if __name__ == "__main__":
    n = 100
    uniDim = np.array([100, 100, 100])
    velRan = np.array([100, 100, 100])
    arrCent = np.array([0, 0, 0])
    masDim = 1e24
    arrBods = np.array([])

    dt = 1e4
    numSteps = 30

    # Generate n bodies
    np.random.seed(0)
    randMas = np.random.random(n)*masDim
    np.random.seed(1)
    randPos = np.random.random([n, 3])*uniDim
    np.random.seed(2)
    randVel = np.random.random([n, 3])*uniDim
    for i in range(n):
        arrBods = np.append(arrBods, tree.body(randMas[i], randPos[i], randVel[i], [0, 0, 0]))

    # arrBods = [tree.body(100, arrCent, [0, 0, 0], [0, 0, 0]),
    #            tree.body(100, [10, 0, 0], [0, 0, 0], [0, 0, 0])]
    start = time.time()
    b = tree.basicRun(arrBods, arrCent, uniDim*2, numSteps, dt)
    end = time.time()
    print(end - start)
    pos = np.array(b[0].getPos())
    # print(pos[:, 0])
    times = np.array(range(0, int(numSteps*dt), int(dt)))
    plt.plot(np.array(b[0].getPos())[:, 0])
    plt.plot(np.array(b[1].getPos())[:, 0])
    #plt.plot(times, np.array(b[1].getPos())[:, 0])
    plt.show()