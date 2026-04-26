import numpy as np
import matplotlib.pyplot as plt
from common import import_data, TOTAL_INFECTIONS_COL, TOTAL_RECOVERIES_COL, INFECTED_COL



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
beta = 0.8826642165086095
gamma = 0.2501883071115894

df = import_data()
N = max(df[TOTAL_INFECTIONS_COL])
R_init = df[TOTAL_RECOVERIES_COL][0]
I_init = df[TOTAL_INFECTIONS_COL][0]
S_init = N - I_init

print(R_init, I_init, S_init)

S, I, R = forward_euler(S_init, I_init, R_init, beta, gamma, 0.001, 20)
N = S + I + R

# Print all 
# TODO: Show correct time in x axis
plt.figure(figsize=(10, 6))
plt.plot(S, marker="", linestyle="-")
plt.plot(I, marker="", linestyle="-")
plt.plot(R, marker="", linestyle="-")
plt.legend(["Susceptable", "Infected", "Recovered"])
plt.title("Total Population Size of Simulation")
plt.xlabel("Time")
plt.ylabel("Individuals")
plt.ylim(0, N[0] * 1.1)  # Sets the bottom to 0 and top to 10% above N
plt.xticks(rotation=45)
plt.grid()
plt.savefig("diagrams/population_size.png")