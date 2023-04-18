import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
	

def set_initial_values(N, V0, delta):
    """Calculates the center of each section, the voltage of each section and if it is a plate or not"""
    center = []
    voltageM = []
    isPlate = []

    for i in range(1, N+1):
        for j in range(1, N+1):
            if (i > (N+1)/2) and (j > (N+1)/2):
                center.append([(j - 0.5)*delta, (i - 0.5)*delta])
                voltageM.append(0)
                isPlate.append(False)
            else:
                center.append([(j - 0.5)*delta, (i - 0.5)*delta])
                isPlate.append(True)
                voltageM.append(V0)

    return center, voltageM, isPlate

def impedance_matrix(N, centers, V, delta, isPlate):
    impedanceM = []
    e0 = 8.85e-12

    for i in range((N**2)):
        impedanceM.append([])
        for j in range((N**2)):
            if i == j:
                impedanceM[i].append((delta/(math.pi * e0))
                                    * (math.log(1 + math.sqrt(2))))
            elif not isPlate[i] or not isPlate[j]:
                impedanceM[i].append(0)
            else:
                impedanceM[i].append(1/(4 * math.pi * e0) * pow(delta, 2) 
                                    / math.sqrt(pow(centers[i][0] - centers[j][0], 2) + pow(centers[i][1] - centers[j][1], 2)))

    return impedanceM


def get_pulses(Z, V):
    return np.linalg.solve(Z, V)


def set_blank(isPlate, a, r):
    """Sets the blank spaces in the plot"""
    for i, b in enumerate(isPlate):
        if not b:
            a[i] = math.nan
            r[i] = [math.nan, math.nan]
    return np.array(a), np.array(r)


def plot_results(N, r, a, V, isPlate):
    a, r = set_blank(isPlate, a, r)

    fig = plt.figure(figsize=(10, 5))
    
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1.plot_surface(np.reshape(r[:, 0], (N, N)), np.reshape(r[:, 1], (N, N)), np.reshape(
        a, (N, N)), cmap=cm.coolwarm, linewidth=0, antialiased=True)
    ax1.set_title('Gráfico 3D')
    
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(a)
    ax2.set_xlabel('N divisões')
    ax2.set_ylabel('Pulsos')
    
    plt.show()


if __name__ == '__main__':
    L = 0.1
    N = 20
    V0 = 10
    delta = L/N

    r, V, iP = set_initial_values(N, V0, delta)
    Z = impedance_matrix(N, r, V, delta, iP)
    a = get_pulses(Z, V)
    print(sum(a))

    plot_results(N, r, a, V, iP)