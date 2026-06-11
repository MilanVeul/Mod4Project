import matplotlib.pyplot as plt
import numpy as np

from data.covid_data import import_covid_data
from forward_euler import forward_euler_sir, forward_euler_seir


# variables in both sir and seir
sigma = 0.19
N_init = 60_000_000
R_init = 0
dt = 0.1

# step/range of beta and gamma
bt_step = 0.010
bt_rnge = 0.100
gm_step = 0.001
gm_rnge = 0.010


def init_ranges(beta, gamma):
    beta_range = np.arange(beta - bt_rnge, beta + bt_rnge + bt_step, bt_step)
    gamma_range = np.arange(gamma - gm_rnge, gamma + gm_rnge + gm_step, gm_step)
    B, G = np.meshgrid(beta_range, gamma_range)
    Z = np.zeros(B.shape)
    return B, G, Z, beta_range, gamma_range


def sir_sensitivity(duration):
    I_init = 0.0023 * N_init
    S_init = N_init - I_init

    beta = 0.1831
    gamma = 0.1000
    scale = 0.00202

    B, G, Z, beta_range, gamma_range = init_ranges(beta, gamma)

    # Fill Z matrix
    for i in range(len(gamma_range)):
        for j in range(len(beta_range)):
            _, I, _ = forward_euler_sir(S_init, I_init, R_init, B[i,j], G[i,j], dt, duration) * scale
            Z[i,j] = np.max(I)

    plot_3d(beta_range, gamma_range, Z, beta, gamma, "diagrams/sensitivity_sir.png")

def seir_sensitivity(duration):
    E_init = 0.0008 * N_init
    I_init = 0.0015 * N_init
    S_init = N_init - I_init - E_init

    beta = 0.1999
    gamma = 0.0370
    scale = 0.00333

    B, G, Z, beta_range, gamma_range = init_ranges(beta, gamma)

    # Fill Z matrix
    for i in range(len(gamma_range)):
        for j in range(len(beta_range)):
            _, _, I, _ = forward_euler_seir(S_init, E_init, I_init, R_init, B[i,j], sigma, G[i,j], dt, duration) * scale
            Z[i,j] = np.max(I)

    plot_3d(beta_range, gamma_range, Z, beta, gamma, "diagrams/sensitivity_seir.png")


def plot_3d(beta_range, gamma_range, Z, beta, gamma, path):
    B, G = np.meshgrid(beta_range, gamma_range)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(B, G, Z, cmap="magma", edgecolor="none", antialiased=True)
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label="Peak Infections")

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
    ax.view_init(elev=15, azim=75)

    plt.tight_layout()
    plt.savefig(path)
    print(f"Figure saved to {path}")

if __name__ == "__main__":
    df = import_covid_data()
    train_duration = int(.35 * len(df))
    sir_sensitivity(train_duration)
    seir_sensitivity(train_duration)
