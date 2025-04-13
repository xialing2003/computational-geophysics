import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

folder = 'results/dirichlet/'
timesteps = range(0,  360 + 1, 40)

time_interval = 0.25

plt.figure(figsize=(10, 15))
for i, itime in enumerate(timesteps):
    plt.subplot(10,1,i+1)
    filename = folder + f'snapshot{itime:05d}'  # e.g., snapshot00001
    data = np.loadtxt(filename)
    x = data[:, 0]
    s = data[:, 1]
    plt.plot(x, s, label=f't={itime * time_interval:.2f} s')
    plt.ylim(-1, 1)
    plt.xlim(0, 100)
    plt.legend(loc='lower left', fontsize='small', ncol=2)

plt.xlabel('X(m)')
# plt.title('Displacement Evolution over Time')
plt.tight_layout()
plt.savefig(folder + 'all_snapshots.png', dpi=300)
