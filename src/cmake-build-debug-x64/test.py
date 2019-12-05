import numpy as np
import matplotlib.pyplot as plt
import treecode as tree
import time

if __name__ == "__main__":
    n = 20
    uniDim = np.array([10, 10, 10])
    velRan = np.array([1, 1, 1])
    arrCent = np.array([0, 0, 0])
    masDim = 1e2
    arrBods = np.array([])

    dt = 1e2
    numSteps = 1000

    # Generate n bodies
    np.random.seed(0)
    randMas = np.ones(n)*masDim # np.random.random(n)*masDim
    np.random.seed(1)
    randPos = np.random.random([n, 3])*uniDim/2
    np.random.seed(2)
    randVel = np.random.random([n, 3])*velRan*0
    for i in range(n):
        arrBods = np.append(arrBods, tree.body(10, randPos[i], randVel[i], [0, 0, 0]))

    # arrBods = [tree.body(100, arrCent, [0, 0, 0], [0, 0, 0]),
    #            tree.body(100, [10, 0, 0], [0, 0, 0], [0, 0, 0])]
    start = time.time()
    b = tree.basicRun(arrBods, arrCent, uniDim*2, numSteps, dt)
    end = time.time()
    # print(end - start)
    pos = np.array(b[0].getPos())
    print(pos[:, 0])
    times = np.array(range(0, int(numSteps*dt), int(dt)))
    for i in range(0, n):
        plt.plot(np.array(b[i].getPos())[:, 0], np.array(b[i].getPos())[:, 1])
        plt.scatter(np.array(b[i].getPos())[-1, 0], np.array(b[i].getPos())[-1, 1])
    # plt.plot(np.array(b[1].getPos())[:, 0])
    #plt.plot(times, np.array(b[1].getPos())[:, 0])
    plt.show()