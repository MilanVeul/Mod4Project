from data.covid_data import import_covid_data, TOTAL_RECOVERIES_COL, INFECTED_COL
from diff import derivative
from scipy.integrate import cumulative_trapezoid
import numpy as np
import matplotlib.pyplot as plt

sigma = 0.2
p = 0.001
N = 60_000_000
E0 = I0 = p * N
S0 = N * (1 - 2 * p)

df = import_covid_data(country="Italy")

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