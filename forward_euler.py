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

# As approximated in assignment 3
beta = 0.4195460684174166
gamma = 0.2647898334098834
dt = 0.01

df = import_data()
N_init = 17_500_000
R_init = df[TOTAL_RECOVERIES_COL][0]
I_init = df[INFECTED_COL][0]
S_init = N_init - I_init

S, I, R = forward_euler(S_init, I_init, R_init, beta, gamma, dt, 20)
N = S + I + R

# Generate Date Axis
start_date = pd.to_datetime(df[DATE_COL].iloc[0])
# Map simulation steps to days relative to start_date
simulation_dates = start_date + pd.to_timedelta(np.arange(len(S)) * dt, unit='D')

# Print all
plt.figure(figsize=(10, 6))
plt.plot(simulation_dates, N, marker="", linestyle="-")
plt.title("Total Population Size of Simulation")
plt.xlabel("Time")
plt.ylabel("Individuals")
plt.ylim(N[0] * 0.99, N[0] * 1.01)
plt.xticks(rotation=45)
plt.tight_layout()
plt.grid()
plt.savefig("diagrams/population_size.png")


plt.figure(figsize=(10, 6))
plt.plot(simulation_dates, I, marker="", linestyle="-")
plt.plot(simulation_dates, R, marker="", linestyle="-")
plt.plot(df["Date"], df[INFECTED_COL], marker="o", linestyle="-")
plt.plot(df["Date"], df[TOTAL_RECOVERIES_COL], marker="o", linestyle="-")
plt.plot()
plt.legend(["Simulated Infected", "Simulated Recovered", "Actual Infected", "Actual Recovered"])
plt.title("Susceptable, Recovered and Infected Individuals over Time")
plt.xlabel("Time")
plt.ylabel("Individuals")
# plt.ylim(0,  * 1.1)  # Sets the bottom to 0 and top to 10% above N
plt.xticks(rotation=45)
plt.tight_layout()
plt.grid()
plt.savefig("diagrams/simulation.png")