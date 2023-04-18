import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import cm


def set_initial_values(N, V0, delta):
    """Calculates the center of each section, the voltage of each section and if it is a plate or not"""
    center = []
    voltageM = []
    isPlate = []

    for k in range(2):
        for i in range(1, N+1):
            for j in range(1, N+1):
                if (i > (N+1)/2) and (j > (N+1)/2):
                    center.append([(j - 0.5)*delta, (i - 0.5)*delta, k*d])
                    voltageM.append(0)
                    isPlate.append(False)
                else:
                    center.append([(j - 0.5)*delta, (i - 0.5)*delta, k*d])
                    isPlate.append(True)
                    if k == 1:
                        voltageM.append(V0)
                    else:
                        voltageM.append(0)

    return center, voltageM, isPlate


def impedance_matrix(N, centers, delta, isPlate):
    """Calculates the impedance matrix of the system"""
    e0 = 8.85e-12
    impedanceM = []

    for i in range(2*(N**2)):
        impedanceM.append([])
        for j in range(2*(N**2)):
            if i == j:
                impedanceM[i].append((delta/(math.pi * e0))
                                    * (math.log(1 + math.sqrt(2))))
            elif not isPlate[i] or not isPlate[j]:
                impedanceM[i].append(0)
            else:
                impedanceM[i].append(1/(4 * math.pi * e0) * pow(delta, 2) 
                                    / math.sqrt(pow(centers[i][0] - centers[j][0], 2) + pow(centers[i][1] - centers[j][1], 2) + pow(centers[i][2] - centers[j][2], 2)))

    return impedanceM


def get_pulses(impedanceM, voltageM):
    """Calculates the pulses of each section"""
    return np.linalg.solve(impedanceM, voltageM)


def set_blank(isPlate, a, r):
    """Sets the blank spaces in the plot"""
    for i, b in enumerate(isPlate):
        if not b:
            a[i] = math.nan
            r[i] = [math.nan, math.nan, math.nan]
    return np.array(a), np.array(r)


def plot_results(N, r, a, V, isPlate):
    """Plots the results"""
    a, r = set_blank(isPlate, a, r)

    fig = plt.figure(figsize=(10, 10))

    # Ambas placas na mesma projeção
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    ax1.plot_surface(np.reshape(r[:len(r)//2, 0], (N, N)), np.reshape(r[:len(r)//2, 1], (N, N)), 
                     np.reshape(np.array(a[:len(a)//2]), (N, N)), cmap=cm.coolwarm, linewidth=0, antialiased=True)
    ax1.plot_surface(np.reshape(r[:len(r)//2, 0], (N, N)), np.reshape(r[:len(r)//2, 1], (N, N)), 
                     np.reshape(np.array(a[len(a)//2:len(a)]), (N, N)), cmap=cm.coolwarm, linewidth=0, antialiased=True)
    ax1.set_title('Projeção 3D Capacitor')

    # Projeção em duas dimensões os coeficientes
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(a)
    ax2.set_xlabel('N divisões')
    ax2.set_ylabel('Pulsos')

    # Projeção individual da placa superior
    ax3 = fig.add_subplot(2, 2, 3, projection='3d')
    ax3.plot_surface(np.reshape(r[:len(r)//2, 0], (N, N)), np.reshape(r[:len(r)//2, 1], (N, N)),
                     np.reshape(np.array(a[:len(a)//2]), (N, N)), cmap=cm.coolwarm, linewidth=0, antialiased=True)
    ax3.set_title('Projeção 3D (Parte superior)')

    # Projeção individual da placa inferior
    ax4 = fig.add_subplot(2, 2, 4, projection='3d')
    ax4.plot_surface(np.reshape(r[:len(r)//2, 0], (N, N)), np.reshape(r[:len(r)//2, 1], (N, N)), 
                     np.reshape(np.array(a[len(a)//2:len(a)]), (N, N)), cmap=cm.coolwarm, linewidth=0, antialiased=True)
    ax4.set_title('Projeção 3D (Parte inferior)')

    plt.show()
    

def paralell_plate_capacitor(L, d, V0, N):
    """Calculates the sup_dist of a paralell plate capacitor"""
    delta = L/N

    r, V, iP = set_initial_values(N, V0, delta)
    Z = impedance_matrix(N, r, delta, iP)
    A = get_pulses(Z, V)

    sup_dist = sum(A)
    print('Aproximação =', sup_dist)

    plot_results(N, r, A, V, iP)
    

def plot_cap_results(cap, acap):
    dif = np.subtract(cap, acap)
    
    plt.plot(cap, label='Valor obtido')
    plt.axhline(y=acap, xmin = 0.044, xmax = 0.956, color='r', linestyle='-', label='Valor analítico')
    plt.plot(dif, label='Diferença')
    plt.grid(True)
    plt.legend()

    plt.show()
    
    
def get_capacitances(L, d, V0):
    """Calculates the capacitances of the capacitor"""
    e0 = 8.85e-12
    anCap = e0*L**2/d
    
    C = []
    N = 5
    while N <= 45:
        delta = L/N
        r, V, iP = set_initial_values(N, V0, delta)
        Z = impedance_matrix(N, r, delta, iP)
        A = get_pulses(Z, V)
        C.append( sum(A[N**2: 2*N**2])*delta**2/V0 )
        N += 5
        
    plot_cap_results(C, anCap) 


if __name__ == "__main__":
    L = 0.1
    d = 0.002
    V0 = 10
    N = 30

    # paralell_plate_capacitor(L, d, V0, N)
    # for N in [8, 16, 32, 48, 64, 82]:
    #      paralell_plate_capacitor(L, d, V0, N)
    
    get_capacitances(L, d, V0)