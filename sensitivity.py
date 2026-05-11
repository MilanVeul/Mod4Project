import numpy as np
import matplotlib.pyplot as plt
from common import import_data,INFECTED_COL, TOTAL_RECOVERIES_COL, DATE_COL
import pandas as pd


def forward_euler(S_init, I_init, R_init, beta, gamma, dt, duration):
    """
    Solves the SIR model using Forward Euler method.
    """
    N = S_init + I_init + R_init
    
    num_steps = int(duration / dt)
    S = np.zeros(num_steps)
    I = np.zeros(num_steps)
    R = np.zeros(num_steps)
    
    # Set initial values
    S[0], I[0], R[0] = S_init, I_init, R_init
    
    for t in range(num_steps - 1):
        S_rate = -(beta * S[t] * I[t]) / N
        I_rate = (beta * S[t] * I[t]) / N - (gamma * I[t])
        R_rate = gamma * I[t]
        
        S[t+1] = S[t] + dt*S_rate
        I[t+1] = I[t] + dt*I_rate
        R[t+1] = R[t] + dt*R_rate
        
    return S, I, R


def get_peak(I):
    return np.max(I)


# As approximated in assignment 3
dt = 0.1
df = import_data()
N_init = 17_500_000
R_init = 0
I_init = 1
S_init = N_init - I_init

beta = 0.420
gamma = 0.265

step = 0.002
beta_range = np.arange(0.32, 0.53, step)
gamma_range = np.arange(0.165, 0.375, step)
B, G = np.meshgrid(beta_range, gamma_range)
Z = np.zeros(B.shape)

# Fill Z matrix
for i in range(len(gamma_range)):
    for j in range(len(beta_range)):
        _, I, _ = forward_euler(S_init, I_init, R_init, B[i,j], G[i,j], dt, 365)
        Z[i,j] = get_peak(I)

# Plotting the Heatmap
plt.figure(figsize=(10, 7))
plt.pcolormesh(B, G, Z, shading='auto', cmap='magma')

plt.colorbar(label='Peak Infections')
plt.title("Peak Infections Sensitivity Analysis")
plt.xlabel("Beta")
plt.ylabel("Gamma")
plt.savefig("diagrams/sensitivity.png")

