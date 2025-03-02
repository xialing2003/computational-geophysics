import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py

def plot(a, i, n):

    plt.subplot(n,1,i)
    plt.plot(range(101),a,'k')
    plt.text(1,1,'time='+str(i*0.4))
    plt.tight_layout()
    


if __name__ == "__main__":

    with h5py.File("T_forward.h5", "r") as f:
        plt.figure(figsize=(9,15))
        for i in range(0,100,10):
            a = f["Time_"+str(i)+".00"][:]
            plot(a, int(i/10)+1, 20) 
        plt.xlabel('Location(m)')
        plt.ylabel('Temperature')
        
        plt.savefig('1.png')

    