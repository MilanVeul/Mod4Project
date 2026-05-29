import numpy as np
import matplotlib.pyplot as plt
from data.meridonia_data import import_meridonia_data
from forward_euler import forward_euler_sir


def get_peak(I):
    return np.max(I)


# As approximated in assignment 3
dt = 0.1
df = import_meridonia_data()
N_init = 17_500_000
R_init = 0
I_init = 1
S_init = N_init - I_init

beta = 0.420
gamma = 0.265

step = 0.01
beta_range = np.arange(0.32, 0.53, step)
gamma_range = np.arange(0.165, 0.375, step)
B, G = np.meshgrid(beta_range, gamma_range)
Z = np.zeros(B.shape)

# Fill Z matrix
for i in range(len(gamma_range)):
    for j in range(len(beta_range)):
        _, I, _ = forward_euler_sir(S_init, I_init, R_init, B[i,j], G[i,j], dt, 365)
        Z[i,j] = get_peak(I)

# Plotting the Heatmap
# plt.figure(figsize=(10, 7))
# plt.pcolormesh(B, G, Z, shading='auto', cmap='magma')

# plt.colorbar(label='Peak Infections')
# plt.title("Peak Infections Sensitivity Analysis")
# plt.xlabel("Beta")
# plt.ylabel("Gamma")
# plt.savefig("diagrams/sensitivity.png")

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection="3d")

# Plot the 3D surface
surf = ax.plot_surface(B, G, Z, cmap="magma", edgecolor="none", antialiased=True)
# Add a color bar for reference
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label="Peak Infections")

# Labels and Titles
ax.set_title("Peak Infections Sensitivity Analysis (3D)")
ax.set_xlabel("Beta")
ax.set_ylabel("Gamma")
ax.set_zlabel("Peak Infections")

beta_idx = np.argmin(np.abs(beta_range - beta))
gamma_idx = np.argmin(np.abs(gamma_range - gamma))
peak = Z[gamma_idx, beta_idx]
print(peak)
scat = ax.scatter(beta, gamma, peak, color="red", s=100, marker="o", label="Estimation", alpha=1.0, zorder=999)
scat.set_depthshade(False)
ax.legend()



# Adjust the viewing angle if needed
ax.view_init(elev=15, azim=45)

plt.tight_layout()
plt.savefig("diagrams/sensitivity_3d.png")

