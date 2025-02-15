import numpy
import pandas as pd
import matplotlib.pyplot as plt
import h5py

def plot(u, v1, v2, T, file_name, ite):
    plt.figure(figsize=(15,15))
    plt.subplot(1,3,1)
    for i in range(ite):
        plt.plot(u[i*10, :] + i,'k')
    plt.title('displacement u(second order)')
    plt.subplot(1,3,2)
    for i in range(ite):
        plt.plot(v2[i*10, :] + i,'b')
    for i in range(ite):
        plt.plot(v1[i*10, :] + i,'r')
    plt.title('velocity v')
    plt.subplot(1,3,3)
    for i in range(ite):
        plt.plot(T[i*10, :] + i,'k')
    plt.title('stress T(first order)')

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

    # plot(u_d, v1_d, v2_d, T_d, 'homo_dirich')
    # plot(u_n, v1_n, v2_n, T_n, 'homo_neumann')

    with h5py.File("hete_data.h5", "r") as f:
        u = f["2order_u"][:]
        v2 = f["2order_v"][:]
        v1 = f["1order_v"][:]
        T = f["1order_T"][:]
    
    print(u.shape)

    plot(u, v1, v2, T, 'hete', 50)