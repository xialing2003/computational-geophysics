import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py

def plot(u, v1, v2, T, file_name, ite, step, dt):

    loc = 0.1 * np.arange(1000 + 1)
    plt.figure(figsize=(9,15))

    plt.subplot(1,3,1)
    for i in range(ite):
        plt.plot(loc, u[i*step, :]*2 + i*step*dt,'k')
    plt.title('displacement(second order)')
    plt.xlabel('Location(m)')
    plt.ylabel('Time(s)')

    plt.subplot(1,3,2)
    for i in range(ite):
        plt.plot(loc, v2[i*step, :]*2 + i*step*dt,'b')
    for i in range(ite):
        plt.plot(loc, v1[i*step, :]*2 + i*step*dt,'r')
    plt.title('velocity')
    plt.xlabel('Location(m)')

    plt.subplot(1,3,3)
    for i in range(ite):
        plt.plot(loc, T[i*step, :]*2 + i*step*dt,'k')
    plt.title('stress(first order)')
    plt.xlabel('Location(m)')

    plt.tight_layout()
    plt.savefig(file_name + '.png')

if __name__ == "__main__":

    with h5py.File("homo_data.h5", "r") as f:
        u_d = f["2order_dirich_u"][:]
        u_n = f["2order_neumann_u"][:]
        v2_d = f["2order_dirich_v"][:]
        v2_n = f["2order_neumann_v"][:]

        v1_d = f["1order_dirich_v"][:]
        v1_n = f["1order_neumann_v"][:]
        T_d = f["1order_dirich_T"][:]
        T_n = f["1order_neumann_T"][:]

    plot(u_d, v1_d, v2_d, T_d, 'homo_dirich', 20, 50, 0.1)
    plot(u_n, v1_n, v2_n, T_n, 'homo_neumann', 20, 50, 0.1)

    with h5py.File("hete_data.h5", "r") as f:
        u = f["2order_u"][:]
        v2 = f["2order_v"][:]
        v1 = f["1order_v"][:]
        T = f["1order_T"][:]
    
    print(u.shape)

    plot(u, v1, v2, T, 'hete', 20, 100, 0.05)