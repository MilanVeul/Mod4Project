# Model: 
# dS/dt = -beta * S * I / N
# dE/dt = beta * S * I / N - sigma * E
# dI/dt = sigma * E - gamma * I
# dR/dt = gamma * I

import numpy as np
from forward_euler import general_forward_euler

initial_beta = 0.3
initial_gamma = 1 / 1.5 # per week
initial_sigma = 1 / (1/7) # per week

# returns (states, derivative) where states[t][i] for variable i at time t
# derivative[t][i] is the derivative for variable i at time t
def forward_euler_seir(
    S_init: float, E_init: float, I_init: float, R_init: float,
    beta: float, gamma: float, sigma: float,
    dt: float, duration: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Solves the SEIR model using Forward Euler method.
    """

    N = S_init + E_init + I_init + R_init
    initial = np.array([S_init, E_init, I_init, R_init])
    
    # state = (S, E, I, R)^T
    def design(state: np.ndarray):
        S, E, I, R = state
        S_rate = -(beta * S * I) / N
        E_rate = (beta * S * I) / N - (sigma * I)
        I_rate = sigma * E - gamma * I
        R_rate = gamma * I
        return np.array([S_rate, E_rate, I_rate, R_rate])
        
    return general_forward_euler(initial, design, dt, duration)

def improve(infected_actual: np.ndarray, S, E, I, R, beta: float, gamma: float, sigma: float):
    N = S + E + I + R
    dt = 0.1
    states, derivatives = forward_euler_seir(S, E, I, R, beta, gamma, sigma, dt, infected_actual.shape[0])
    
    S, E, I, R = states
    infected_actual = infected_actual - R
    S_rate, E_rate, I_rate, R_rate = derivatives

    beta = np.mean(np.reciprocal(S * infected_actual) * -S_rate * N)
    gamma = np.mean(R_rate * np.reciprocal(infected_actual))
    sigma = np.mean(-(S_rate + E_rate) * np.reciprocal(infected_actual))

    return beta, gamma, sigma