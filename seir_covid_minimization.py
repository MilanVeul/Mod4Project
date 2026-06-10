import numpy as np
from forward_euler import forward_euler_seir
from data.covid_data import import_covid_data, INFECTED_COL, TOTAL_RECOVERIES_COL, DATE_COL
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def fit_seir_minimization(df, train_fraction, dt):

    sigma = 1/5.1

    N = 60_000_000
    I_init = 0.0015 * N
    E_init = 0.0008 * N
    R_init = 0.0
    S_init = N - I_init - E_init - R_init

    initial_beta = 0.420
    initial_gamma = 0.08
    initial_scale = 0.001

    duration = len(df)
    train_duration = int(train_fraction * len(df))

    train_infections = df[INFECTED_COL][:train_duration]
    train_recoveries = df[TOTAL_RECOVERIES_COL][:train_duration]


    def loss_function(params: np.ndarray) -> float:
        beta, gamma, scale = params
        
        S, E, I, R = forward_euler_seir(S_init, E_init, I_init, R_init, beta, sigma, gamma, dt, train_duration) * scale
        rmse_infections = np.sqrt(np.mean((train_infections - I[::int(1 / dt)] - E[::int(1 / dt)]) ** 2))
        rmse_recoveries = np.sqrt(np.mean((train_recoveries - R[::int(1 / dt)]) ** 2))

        loss = rmse_infections + rmse_recoveries
        print(f"loss infections={rmse_infections:.4f}, loss recoveries={rmse_recoveries:.4f} beta={beta}, gamma={gamma}, scale={scale}")
        return loss

    # The optimization result represented as a OptimizeResult object.
    # Important attributes are:
    # x the solution array,
    # success a Boolean flag indicating if the optimizer exited successfully and message which describes the cause of the termination.
    # See OptimizeResult for a description of other attributes.
    result = minimize(
        loss_function,
        x0=[initial_beta, initial_gamma, initial_scale],
        method="Nelder-Mead",
        options={"maxiter": 400, "xatol": 1e-4, "fatol": 1e-5, "disp": True},
        bounds=(
            (0.1, 1.0),
            (0.0, 0.1),
            (0.0, 0.1),
        )
    )

    beta, gamma, scale = result.x

    S, E, I, R = forward_euler_seir(S_init, E_init, I_init, R_init, beta, sigma, gamma, dt, duration) * scale
    # A, I_for_scale = df[INFECTED_COL], I[::int(1 / dt)]
    # scale = np.dot(A, I_for_scale) / (np.linalg.norm(I_for_scale) ** 2)

    S, I, R, E = S[::int(1 / dt)], I[::int(1 / dt)], R[::int(1 / dt)], E[::int(1 / dt)]

    return S, E, I, R, beta, gamma, scale

