from typing import Callable
import numpy as np
import matplotlib.pyplot as plt
from data.meridonia_data import import_meridonia_data, INFECTED_COL, TOTAL_RECOVERIES_COL, DATE_COL
import pandas as pd

def general_forward_euler(initial: np.ndarray, design: Callable[[np.ndarray], np.ndarray], dt: float, duration: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Solves a general compartmental model using the Forward Euler method.
    """
    num_steps = int(duration / dt)
    states = np.zeros((num_steps, len(initial)))
    rates = np.zeros((num_steps - 1, len(initial)))

    for j in range(len(initial)):
        states[0][j] = initial[j]
    
    for t in range(num_steps - 1):
        derivative = design(states[t])
        rates[t] = derivative
        states[t + 1] = states[t] + dt * derivative
        
    return states, rates

def forward_euler_sir(S_init: float, I_init: float, R_init: float, beta: float, gamma: float, dt: float, duration: float) -> np.ndarray:
    """
    Solves the SIR model using Forward Euler method.
    """

    N = S_init + I_init + R_init
    initial = np.array([S_init, I_init, R_init])
    
    # state = (S, I, R)^T
    def design(state: np.ndarray):
        S, I, R = state
        S_rate = -(beta * S * I) / N
        I_rate = (beta * S * I) / N - (gamma * I)
        R_rate = gamma * I
        return np.array([S_rate, I_rate, R_rate])
        
    return general_forward_euler(initial, design, dt, duration)[0].T

if __name__ == "__main__":
    # As approximated in assignment 3
    beta = 0.4195464445256309
    gamma = 0.26479002144073366
    dt = 0.1
    
    df = import_meridonia_data()
    # N_init = 17_500_000
    N_init = 1_000_000
    R_init = df[TOTAL_RECOVERIES_COL][0]
    I_init = df[INFECTED_COL][0]
    S_init = N_init - I_init
    
    S, I, R = forward_euler_sir(S_init, I_init, R_init, beta, gamma, dt, 20)
    N = S + I + R
    
    # Generate Date Axis
    start_date = pd.to_datetime(df[DATE_COL].iloc[0])
    # Map simulation steps to days relative to start_date
    simulation_dates = start_date + pd.to_timedelta(np.arange(len(S)) * dt, unit='D')
    N_trend = N - np.mean(N)
    
    # Print all
    # plt.figure(figsize=(10, 6))
    # plt.plot(simulation_dates, N_trend, marker="", linestyle="-")
    # plt.title("Population Change during Simulatioin")
    # plt.xlabel("Time")
    # plt.ylabel("Individuals")
    # plt.ylim(-0.2, 0.2)
    # plt.xticks(rotation=45)
    # plt.tight_layout()
    # plt.grid()
    # plt.savefig("diagrams/population_size.png")
    
    
    plt.figure(figsize=(10, 6))
    plt.plot(simulation_dates, I, color="C0", marker="", linestyle="-")
    plt.plot(simulation_dates, R, color="C1", marker="", linestyle="-")
    plt.plot(df["Date"], df[INFECTED_COL], color="C0", marker="o", linestyle="-")
    plt.plot(df["Date"], df[TOTAL_RECOVERIES_COL], color="C1", marker="o", linestyle="-")
    plt.plot()
    plt.legend(["Simulated Infected", "Simulated Recovered", "Actual Infected", "Actual Recovered"])
    # plt.legend(["Simulated Infected", "Simulated Recovered"])
    plt.title("Susceptable, Recovered and Infected Individuals over Time")
    plt.xlabel("Time")
    plt.ylabel("Individuals")
    # plt.ylim(0,  * 1.1)  # Sets the bottom to 0 and top to 10% above N
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.grid()
    plt.savefig("diagrams/simulation_vs_actual_v2.png")