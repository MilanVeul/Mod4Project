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
dt = 0.01
df = import_data()
N_init = 17_500_000
R_init = 0
I_init = 1
S_init = N_init - I_init

beta = 0.420
gamma = 0.265

betas = np.arange(0.32, 0.53, 0.005)
beta_peaks = []
for beta_2 in betas:
    S, I, R = forward_euler(S_init, I_init, R_init, beta_2, gamma, dt, 365)
    beta_peaks.append(get_peak(I))

gammas = np.arange(0.165, 0.375, 0.005)
gamma_peaks = []
for gamma_2 in gammas:
    S, I, R = forward_euler(S_init, I_init, R_init, beta, gamma_2, dt, 365)
    gamma_peaks.append(get_peak(I))

print(beta_peaks[20] - beta_peaks[19])
print(gamma_peaks[1] - gamma_peaks[0])
    
# beta
plt.figure(figsize=(10, 6))
plt.plot(betas, beta_peaks, marker="o", linestyle="-")
plt.title("Value of Beta vs Peak Infections ")
plt.xlabel("Beta")
plt.ylabel("Peak")
plt.tight_layout()
plt.grid()
plt.savefig("diagrams/beta_sensitivity.png")


# gamma
plt.figure(figsize=(10, 6))
plt.plot(gammas, gamma_peaks, marker="o", linestyle="-")
plt.title("Value of Gamma vs Peak Infections ")
plt.xlabel("Gamma")
plt.ylabel("Peak")
plt.tight_layout()
plt.grid()
plt.savefig("diagrams/gamma_sensitivity.png")

