from forward_euler import evaluate_seir_parameters
from data.meridonia_data import import_meridonia_data, TOTAL_RECOVERIES_COL, INFECTED_COL
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from diff import derivative

def evaluate_sir(
    df: pd.DataFrame, S: np.ndarray, I: np.ndarray, R: np.ndarray, N: float,
    S_init: float, I_init: float, R_init: float,
    dt: float, duration: float=200
):
    infection_rate = derivative(I)
    recovery_rate = derivative(R)
    susceptible_rate = derivative(S)

    plt.figure(figsize=(10, 6))
    plt.plot(df["Date"], infection_rate, marker="o", linestyle="-")
    plt.plot(df["Date"], recovery_rate, marker="o", linestyle="-")
    plt.title("Approximation of $\\frac{dI}{dt}$ and $\\frac{dR}{dt}$")
    plt.xlabel("Date")
    plt.ylabel("Rate of change")
    plt.xticks(rotation=45)
    plt.grid()
    plt.legend(["Infection rate", "Recovery rate"])
    plt.savefig("diagrams/approximation_derivative.png")

    A = np.vstack((
        np.column_stack((-S * I / N, np.zeros_like(I))),
        np.column_stack((S * I / N, -I)),
        np.column_stack((np.zeros_like(I), I))
    ))

    b = np.hstack((-infection_rate - recovery_rate, infection_rate, recovery_rate))
    res, residuals, _, _ = np.linalg.lstsq(A, b, rcond=None)
    beta, gamma = res

    print(f"Beta: {beta}")
    print(f"Gamma: {gamma}")
    print(f"Median approximation of beta: {np.median((infection_rate + recovery_rate) * N / S / I)}")
    print(f"Mean approximation of gamma: {np.mean(recovery_rate / I)}")
    print(f"R_0: {beta / gamma}")

    plt.figure(figsize=(10, 6))
    plt.plot(df["Date"], (A @ res - b)[:len(df)], marker="o", linestyle="-")
    # plt.plot(df["Date"], (A @ res - b)[len(df):], marker="o", linestyle="-")
    plt.title("Residuals of the least squares approximation")
    plt.xlabel("Date")
    plt.ylabel("Residual")
    plt.xticks(rotation=45)
    plt.grid()
    plt.legend(["Infection rate residual", "Recovery rate residual"])
    plt.savefig("diagrams/residuals.png")

    evaluate_seir_parameters(df, S_init, I_init, R_init, 0, beta, gamma, None, dt, duration)
    

if __name__ == "__main__":
    df = import_meridonia_data()

    N = 1_000_000
    I = df[INFECTED_COL].values
    R = df[TOTAL_RECOVERIES_COL].values
    S = N - I - R

    dt = 0.1
    I_init = df[INFECTED_COL].iloc[0]
    R_init = df[TOTAL_RECOVERIES_COL].iloc[0]
    S_init = N - I_init - R_init

    evaluate_sir(df, S, I, R, N, S_init, I_init, R_init, dt, 22)