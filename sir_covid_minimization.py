import numpy as np
from forward_euler import forward_euler_sir
from data.covid_data import import_covid_data, INFECTED_COL, TOTAL_RECOVERIES_COL, DATE_COL
from scipy.optimize import minimize
from seir_covid_minimization import fit_seir_minimization
import matplotlib.pyplot as plt

def fit_sir_minimization(df):


    N = 60_000_000
    I_init = 0.0015 * N
    R_init = 0.0
    S_init = N - I_init - R_init

    initial_beta = 0.420
    initial_gamma = 0.08
    initial_scale = 0.001

    dt = 0.1
    duration = len(df)

    def loss_function(params: np.ndarray) -> float:
        beta, gamma, scale = params
        
        S, I, R = forward_euler_sir(S_init, I_init, R_init, beta, gamma, dt, duration) * scale
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

    S, I, R = forward_euler_sir(S_init, I_init, R_init, beta, gamma, dt, duration)
    # A, I_for_scale = df[INFECTED_COL], I[::int(1 / dt)]
    # scale = np.dot(A, I_for_scale) / (np.linalg.norm(I_for_scale) ** 2)

    I, R = I[::int(1 / dt)], R[::int(1 / dt)]

    return S, I, R, beta, gamma, scale

# Comparison
if __name__ == "__main__":
    df = import_covid_data()

    S_sir, I_sir, R_sir, beta_sir, gamma_sir, scale_sir = fit_sir_minimization(df)
    S_seir, E_seir, I_seir, R_seir, beta_seir, gamma_seir, scale_seir = fit_seir_minimization(df)
    print(f"  Actual peak infections   : {np.max(df[INFECTED_COL])} on {df[DATE_COL][np.argmax(df[INFECTED_COL])]}")
    print("SIR:")
    print(f"  beta={beta_sir:.4f}, gamma={gamma_sir:.4f}, scale={scale_sir:.5f}")
    print(f"  Simulated peak infections: {np.max(I_sir)} on {df[DATE_COL][np.argmax(I_sir)]}")
    print("\nSEIR")
    print(f" beta={beta_seir:.4f}, gamma={gamma_seir:.4f}, scale={scale_seir:.5f}")
    print(f"  Simulated peak infections: {np.max(I_seir)} on {df[DATE_COL][np.argmax(I_seir)]}")
    

    plt.figure(figsize=(10, 6))
    plt.plot(df[DATE_COL], df[INFECTED_COL], color="C0",  linestyle="--", linewidth=2)
    plt.plot(df[DATE_COL], I_seir * scale_seir, color="C1", linestyle="-", linewidth=2)
    plt.plot(df[DATE_COL], I_sir * scale_sir, color="C2", linestyle="-", linewidth=2)
    plt.legend(["Actual Infections", "Simulated Infections by SEIR", "Simulated Infections by SIR"])
    plt.title("Comparison of SEIR and SIR on Current Infections")
    plt.xlabel("Date")
    plt.ylabel("Current Infections")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.grid()
    plt.savefig("diagrams/seir_sir_comparison.png")
