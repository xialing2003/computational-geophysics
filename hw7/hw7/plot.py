import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

folder = 'results/case2/'

## model 1
timesteps = range(1, 180001 + 1, 20000)

## model 2
# timesteps = range(1, 36001 + 1, 4000)

time_interval = 10000/60/60/24/365
cmap = cm.get_cmap('rainbow_r', len(timesteps))

plt.figure(figsize=(10, 6))
for i, itime in enumerate(timesteps):
    filename = folder + f'snapshot{itime:06d}'  # e.g., snapshot00001
    try:
        data = np.loadtxt(filename)
        x = data[:, 0]
        T = data[:, 1]
        plt.plot(x, T, color=cmap(i), label = f't={itime * time_interval:.2f} kyr')
        plt.scatter(x, T, color=cmap(i), s=5)
    except Exception as e:
        print(f"Could not read file: {filename}, error: {e}")

itime = 480001
filename = folder + f'snapshot{itime:06d}'  # e.g., snapshot00001
data = np.loadtxt(filename)
x = data[:, 0]
T = data[:, 1]
plt.plot(x, T, color='violet', label = f't={itime * time_interval:.2f} kyr')
plt.scatter(x, T, color='violet', s=5)

## model 1
T = 10*x/3000
plt.plot(x, T, color='black', linestyle='dashed', label='analytical')

## model 2
# N = 1000
# x = np.linspace(0, 3000, N)
# T = np.zeros_like(x)
# mid = N // 2
# T[:mid+1] = (10/6) * x[:mid+1] / 1500
# T[mid+1:] = (10/6) + (10 - 10/6) * (x[mid+1:] - 1500) / 1500
# plt.plot(x, T, color='black', linestyle='dashed', label='analytical')

plt.xlabel('X(m)')
plt.ylabel('Temperature(K)')
plt.title('Temperature Evolution over Time')
plt.legend(loc='best', fontsize='small', ncol=2) 
plt.tight_layout()

plt.savefig(folder + 'all_snapshots.png', dpi=300)
