import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from forward_euler import general_forward_euler
from data.influenza_data import import_influenza_data, NEW_INFECTIONS_COL, DATE_COL


def forward_euler_seir(S_init: float, E_init: float, I_init: float, R_init: float,
                       beta: float, sigma: float, gamma: float,
                       dt: float, duration: float) -> np.ndarray:
    """
    Solves the SEIR model using Forward Euler method.
    """
    N = S_init + E_init + I_init + R_init
    initial = np.array([S_init, E_init, I_init, R_init])

    def design(state: np.ndarray):
        S, E, I, R = state
        S_rate = -(beta * S * I) / N
        E_rate = (beta * S * I) / N - sigma * E
        I_rate = sigma * E - gamma * I
        R_rate = gamma * I
        return np.array([S_rate, E_rate, I_rate, R_rate])

    return general_forward_euler(initial, design, dt, duration)[0].T


# Influenza incubation period ~ 2 days (according to the literature),
#  so sigma = 1/2
sigma = 1 / 2

# Initial guess for beta, gamma is taken from the SIR result from part 1
beta = 0.420
gamma = 0.265

# Now we are picking one full flu season in HHS Region 4 (Oct 2018 - May 2019)
print("Loading influenza data...", flush=True)
df = import_influenza_data(region=4, start_wk=40, start_yr=2018, end_wk=20, end_yr=2019)
df = df.reset_index(drop=True)

new_cases_raw = df[NEW_INFECTIONS_COL].values
# We are subtracting the background ILI level (week-0 baseline) so the season starts at approx 0
# which makes a SEIR outbreak start small instead of jumping in mid-season
baseline = new_cases_raw[0]
new_cases = new_cases_raw - baseline
# remove negatives due to weekly noise
new_cases = np.clip(new_cases, 0, None)  
num_weeks = len(new_cases)
print(f"Loaded {num_weeks} weeks of data. Baseline (week 0): {baseline:.0f}, max excess: {new_cases.max():.2f}", flush=True)

dt = 0.5
duration_days = num_weeks * 7 + 1 
steps_per_week = int(round(7 / dt))


# For a flu-like outbreak (R_0 ~ 1.5), peak weekly cases is roughly 0.15 * N,
# so N ~ 10 * peak_data gives a model curve of similar size.
peak_data = new_cases.max()
N = 10 * peak_data
print(f"Data peak (excess): {peak_data:.0f}, chosen N: {N:.0f}", flush=True)

# Initial conditions: epidemic just starting (small seed of infections)
I_init = 100.0
E_init = 100.0
R_init = 0.0
S_init = N - I_init - E_init - R_init


def predicted_weekly_cases(beta, gamma):
    S, E, I, R = forward_euler_seir(S_init, E_init, I_init, R_init,
                                     beta, sigma, gamma, dt, duration_days)
    # Everyone who has ever been infectious by time t is I(t) + R(t)
    # weekly new infections = jump in (I + R) over each week
    cumulative = I + R
    weekly_cumulative = cumulative[::steps_per_week][:num_weeks + 1]
    return np.diff(weekly_cumulative)


call_count = 0

def loss(params):
    global call_count
    call_count += 1
    beta, gamma = params
    # Biological bounds for beta and gamma:
    # beta in [0.1, 1.5]  (per-day infection rate)
    # gamma in [0.14, 0.5]  (infectious period 2 - 7 days, according to the literature)
    if not (0.1 <= beta <= 1.5 and 0.14 <= gamma <= 0.5):
        return 1e6
    predicted = predicted_weekly_cases(beta, gamma)
    L = np.mean(((new_cases - predicted) / peak_data) ** 2)
    print(f"  [call {call_count}] beta={beta:.4f}, gamma={gamma:.4f}, loss={L:.6f}", flush=True)
    return L

initial_pred = predicted_weekly_cases(beta, gamma)
print(f"Initial model peak: {initial_pred.max():.0f} (data peak: {peak_data:.0f})", flush=True)
print(f"Initial loss at beta={beta}, gamma={gamma}: {loss([beta, gamma]):.6f}", flush=True)
call_count = 0

result = minimize(
    loss,
    x0=[beta, gamma],
    method="Nelder-Mead",
    options={"maxiter": 400, "xatol": 1e-4, "fatol": 1e-5, "disp": True}
)

beta_fit, gamma_fit = result.x
print(f"Fitted beta:  {beta_fit}", flush=True)
print(f"Fitted gamma: {gamma_fit}", flush=True)
print(f"R_0:          {beta_fit / gamma_fit}", flush=True)

predicted = predicted_weekly_cases(beta_fit, gamma_fit)
predicted_with_baseline = predicted + baseline

plt.figure(figsize=(10, 6))
plt.plot(df[DATE_COL], new_cases_raw, color="C0", marker="o", linestyle="-")
plt.plot(df[DATE_COL], predicted_with_baseline, color="C1", marker="", linestyle="-")
plt.legend(["Actual New Infections (ILINet)", "SEIR Fit"])
plt.title("SEIR Fit to Influenza Data (2018-2019, Region 4)")
plt.xlabel("Date")
plt.ylabel("Weekly New Infections")
plt.xticks(rotation=45)
plt.tight_layout()
plt.grid()
plt.savefig("diagrams/seir_fit_influenza.png")
print("Saved diagrams/seir_fit_influenza.png", flush=True)
