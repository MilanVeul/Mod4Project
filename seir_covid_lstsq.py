from seir_covid_minimization import fit_seir_nelder_mead
from forward_euler import evaluate_seir_parameters
from data.covid_data import import_covid_data, TOTAL_RECOVERIES_COL, INFECTED_COL, DATE_COL, NEW_INFECTIONS_COL
from diff import derivative
from scipy.integrate import cumulative_trapezoid
import numpy as np
import matplotlib.pyplot as plt

df = import_covid_data(country="Italy")

sigma = 0.5
# p = 0.0001
# N = 1_000_000
N = 10 * np.max(df[NEW_INFECTIONS_COL])
I0 = df[INFECTED_COL].iloc[0]
R0 = df[TOTAL_RECOVERIES_COL].iloc[0]
E0 = I0 / 2
S0 = N - I0 - R0 - E0
# E0 = I0 = p * N
# S0 = N * (1 - 2 * p)


R = df[TOTAL_RECOVERIES_COL]
dRdt = derivative(df[TOTAL_RECOVERIES_COL])
I = df[INFECTED_COL]
dIdt = derivative(df[INFECTED_COL])

E = (dIdt + dRdt) / sigma
dEdt = derivative(E)
dSdt = -(dEdt + dIdt + dRdt)
S = S0 + cumulative_trapezoid(dSdt, dx=1, initial=0)

# A = np.vstack((
#     # np.column_stack((-S * I / N, np.zeros_like(I))),
#     np.column_stack((S * I / N, np.zeros_like(I))),
#     np.column_stack((np.zeros_like(I), -I)),
#     # np.column_stack((np.zeros_like(I), I))
# ))

# b = np.hstack((
#     # dSdt,
#     -dSdt,
#     -dRdt, # dIdt - (dIdt + dRdt),
#     # dRdt
# ))

# res, residuals, _, _ = np.linalg.lstsq(A, b, rcond=None)
# beta, gamma = res

A = np.vstack((
    np.column_stack((-S * I / N, np.zeros_like(I), np.zeros_like(I))),
    np.column_stack((S * I / N, np.zeros_like(I), -E)),
    np.column_stack((np.zeros_like(I), -I, E)),
    np.column_stack((np.zeros_like(I), I, np.zeros_like(I)))
))

b = np.hstack((
    dSdt,
    dEdt,
    dIdt, # dIdt - (dIdt + dRdt),
    dRdt
))

res, residuals, _, _ = np.linalg.lstsq(A, b, rcond=None)
beta, gamma, sigma = res

print(f"Beta: {beta}")
print(f"Gamma: {gamma}")
print(f"Sigma: {sigma}")
# print(f"Median approximation of beta: {np.median((infection_rate + recovery_rate) * N / S / I)}")
# print(f"Mean approximation of gamma: {np.mean(recovery_rate / I)}")
print(f"R_0: {beta / gamma}")
print(f"MSE: {np.mean((A @ res - b) ** 2)}")
print(f"RMSE: {np.sqrt(np.mean((A @ res - b) ** 2))}")

# plt.figure(figsize=(10, 6))
# plt.plot(df["Date"], (A @ res - b)[:len(df)], marker="o", linestyle="-")
# # plt.plot(df["Date"], (A @ res - b)[len(df):], marker="o", linestyle="-")
# plt.title("Residuals of the least squares approximation")
# plt.xlabel("Date")
# plt.ylabel("Residual")
# plt.xticks(rotation=45)
# plt.grid()
# plt.legend(["Infection rate residual", "Recovery rate residual"])
# plt.savefig("diagrams/residuals.png")

evaluate_seir_parameters(df, S0, I0, 0, E0, beta, gamma, sigma, 0.1)
fit_seir_nelder_mead(
    df=df,
    date_col=DATE_COL,
    observed_col=NEW_INFECTIONS_COL,
    sigma=sigma,
    initial_guess=(0.420, 0.265),
    # initial_guess=(beta, gamma),
    S_init=S0,
    E_init=E0,
    I_init=I0,
    R_init=R0,
    dt=0.1,
    time_scale_factor = 1,
    output_path="diagrams/seir_fit_covid.png",
    title="SEIR Fit to Covid Data",
    beta_bounds=(0.01, 1.5),
    gamma_bounds=(0.01, 0.5)
)