from common import import_data, TOTAL_INFECTIONS_COL, TOTAL_RECOVERIES_COL, INFECTED_COL
import numpy as np
import matplotlib.pyplot as plt

def forward_difference(data: np.ndarray, i: int):
    return data[i+1] - data[i]

def symmetric_difference(data: np.ndarray, i: int):
    return (data[i+1] - data[i-1]) / 2

def derivative(data: np.ndarray):
    derivative = np.zeros_like(data)
    for i in range(1, len(data) - 1):
        derivative[i] = symmetric_difference(data, i)

    derivative[len(data) - 1] = forward_difference(data, len(data) - 2)

    return derivative

df =import_data()
infection_rate = derivative(df[TOTAL_INFECTIONS_COL].values)
recovery_rate = derivative(df[TOTAL_RECOVERIES_COL].values)

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

N = max(df[TOTAL_INFECTIONS_COL])
I = df[INFECTED_COL].values
R = df[TOTAL_RECOVERIES_COL].values
S = N - I

# minimize |Ax - b|^2, where b = [infection_rate, recovery_rate]
# infection_rate = beta / N * S * I - gamma * I
# recovery_rate = gamma * I
# A = [[S * I / N, -I], [0, I]]
# x = [beta, gamma]

A = np.vstack(
    (np.column_stack((S * I / N, -I)),
    np.column_stack((np.zeros_like(I), I)))
)

b = np.hstack((infection_rate, recovery_rate))
res, residuals, _, _ = np.linalg.lstsq(A, b, rcond=None)
beta, gamma = res

print(f"Beta: {beta}")
print(f"Gamma: {gamma}")
print(f"Mean approximation of gamma: {np.mean(recovery_rate / I)}")
print(f"Median approximation of beta: {np.median((infection_rate + recovery_rate) * N / S / I)}")
print(f"R_0: {beta / gamma}")

plt.figure(figsize=(10, 6))
plt.plot(df["Date"], (A @ res - b)[:len(df)], marker="o", linestyle="-")
plt.plot(df["Date"], (A @ res - b)[len(df):], marker="o", linestyle="-")
plt.title("Residuals of the least squares approximation")
plt.xlabel("Date")
plt.ylabel("Residual")
plt.xticks(rotation=45)
plt.grid()
plt.legend(["Infection rate residual", "Recovery rate residual"])
plt.savefig("diagrams/residuals.png")