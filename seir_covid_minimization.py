import numpy as np
from forward_euler import forward_euler_seir
from data.covid_data import import_covid_data, INFECTED_COL, TOTAL_RECOVERIES_COL, DATE_COL
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def fit_seir_minimization(df):

    sigma = 1/5.1

    N = 60_000_000
    I_init = 0.0015 * N
    E_init = 0.0008 * N
    R_init = 0.0
    S_init = N - I_init - E_init - R_init

    initial_beta = 0.420
    initial_gamma = 0.08
    initial_scale = 0.001

    dt = 0.1
    duration = len(df)

    def loss_function(params: np.ndarray) -> float:
        beta, gamma, scale = params
        
        S, E, I, R = forward_euler_seir(S_init, E_init, I_init, R_init, beta, sigma, gamma, dt, duration) * scale
        rmse_infections = np.sqrt(np.mean((df[INFECTED_COL] - I[::int(1 / dt)]) ** 2))
        rmse_recoveries = np.sqrt(np.mean((df[TOTAL_RECOVERIES_COL] - R[::int(1 / dt)]) ** 2))

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

    S, E, I, R = forward_euler_seir(S_init, E_init, I_init, R_init, beta, sigma, gamma, dt, duration)
    # A, I_for_scale = df[INFECTED_COL], I[::int(1 / dt)]
    # scale = np.dot(A, I_for_scale) / (np.linalg.norm(I_for_scale) ** 2)

    I, R = I[::int(1 / dt)], R[::int(1 / dt)]

    return S, E, I, R, beta, gamma, scale


if __name__ == "__main__":
    df = import_covid_data()

    S, E, I, R, beta, gamma, sigma, scale = fit_seir_minimization(df)
    print(f"beta={beta:.4f}, gamma={gamma:.4f}, scale={scale:.5f}")
    print(f"Actual peak infections   : {np.max(df[INFECTED_COL])} on {df[DATE_COL][np.argmax(df[INFECTED_COL])]}")
    print(f"Simulated peak infections: {np.max(I)} on {df[DATE_COL][np.argmax(I)]}")

    plt.figure(figsize=(10, 6))
    plt.plot(df[DATE_COL], df[INFECTED_COL], color="C0",  linestyle="--", linewidth=2)
    plt.plot(df[DATE_COL], I * scale, color="C0", linestyle="-", linewidth=2)
    plt.plot(df[DATE_COL], df[TOTAL_RECOVERIES_COL], color="C1", linestyle="--", linewidth=2)
    plt.plot(df[DATE_COL], R* scale, color="C1", linestyle="-", linewidth=2)
    plt.legend(["Actual Infections", "Simulated Infections", "Actual Recoveries", "Simulated Recoveries"])
    plt.title("SEIR Fit to Covid Data (Italy)")
    plt.xlabel("Date")
    plt.ylabel("Current Infections")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.grid()
    plt.savefig("diagrams/seir_fit_covid.png")
    print("Saved diagrams/seir_fit_covid.png", flush=True)
