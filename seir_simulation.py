from forward_euler import forward_euler_seir
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from data.influenza_data import import_influenza_data, NEW_INFECTIONS_COL, DATE_COL


def predicted_weekly_cases(beta, gamma, S_init, E_init, I_init, R_init, sigma, dt, duration_days, timescale, num_times):
    S, E, I, R = forward_euler_seir(S_init, E_init, I_init, R_init,
                                     beta, sigma, gamma, dt, duration_days)
    # Everyone who has ever been infectious by time t is I(t) + R(t)
    # weekly new infections = jump in (I + R) over each week
    cumulative = I + R
    weekly_cumulative = cumulative[::timescale][:num_times + 1]
    return np.diff(weekly_cumulative)

call_count = 0
def loss_function(params, S_init, E_init, I_init, R_init, sigma, dt, duration_days, 
                  timescale, num_times, new_cases, peak_data, beta_bounds, gamma_bounds):
    global call_count
    call_count += 1

    beta, gamma = params
    # Biological bounds for beta and gamma:
    # beta in [0.1, 1.5]  (per-day infection rate)
    # gamma in [0.14, 0.5]  (infectious period 2 - 7 days, according to the literature)
    if not (beta_bounds[0] <= beta <= beta_bounds[1] and gamma_bounds[0] <= gamma <= gamma_bounds[1]):
        return 1e6
    predicted = predicted_weekly_cases(beta, gamma, S_init, E_init, I_init, R_init, 
                                       sigma, dt, duration_days, timescale, num_times)
    L = np.mean(((new_cases - predicted) / peak_data) ** 2)
    print(f"  [call {call_count}] beta={beta:.4f}, gamma={gamma:.4f}, loss={L:.6f}", flush=True)
    return L


def fit_seir_nelder_mead(df, date_col, observed_col, sigma, initial_guess, S_init, E_init, I_init, R_init, dt, time_scale_factor: int = 7,
                      beta_bounds=(0.1, 1.5), gamma_bounds=(0.14, 0.5), output_path=None, title=None):
    df = df.reset_index(drop=True)
    new_cases_raw = df[observed_col].values
    baseline = new_cases_raw[0]
    new_cases = np.clip(new_cases_raw - baseline, 0, None)
    num_times = len(new_cases)
    
    dt_val = dt
    duration_days = num_times * time_scale_factor + 1
    timescale = int(round(time_scale_factor / dt_val))
    
    peak_data = new_cases.max()
    
    initial_pred = predicted_weekly_cases(initial_guess[0], initial_guess[1], S_init, E_init, I_init, R_init,
                                           sigma, dt_val, duration_days, timescale, num_times)
    print(f"Loaded {num_times} rows of data. Baseline: {baseline:.0f}, max excess: {new_cases.max():.2f}", flush=True)
    print(f"Initial model peak: {initial_pred.max():.0f} (data peak: {peak_data:.0f})", flush=True)
    print(f"Initial loss at beta={initial_guess[0]}, gamma={initial_guess[1]}: {loss_function(initial_guess, S_init, E_init, I_init, R_init, sigma, dt_val, duration_days, timescale, num_times, new_cases, peak_data, beta_bounds, gamma_bounds):.6f}", flush=True)
    
    result = minimize(
        lambda p: loss_function(p, S_init, E_init, I_init, R_init, sigma, dt_val, duration_days,
                               timescale, num_times, new_cases, peak_data, beta_bounds, gamma_bounds),
        x0=initial_guess,
        method="Nelder-Mead",
        options={"maxiter": 400, "xatol": 1e-4, "fatol": 1e-5, "disp": True}
    )
    
    beta_fit, gamma_fit = result.x
    print(f"Fitted beta:  {beta_fit}", flush=True)
    print(f"Fitted gamma: {gamma_fit}", flush=True)
    print(f"R_0:          {beta_fit / gamma_fit}", flush=True)
    
    predicted = predicted_weekly_cases(beta_fit, gamma_fit, S_init, E_init, I_init, R_init,
                                       sigma, dt_val, duration_days, timescale, num_times)
    predicted_with_baseline = predicted + baseline
    
    plt.figure(figsize=(10, 6))
    plt.plot(df[date_col], new_cases_raw, color="C0", marker="o", linestyle="-")
    plt.plot(df[date_col], predicted_with_baseline, color="C1", marker="", linestyle="-")
    plt.legend([f"Actual {observed_col}", "SEIR Fit"])
    plt.title(title or f"SEIR Fit to {observed_col}")
    plt.xlabel("Date")
    plt.ylabel(observed_col)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.grid()
    
    if output_path:
        plt.savefig(output_path)
        print(f"Saved {output_path}", flush=True)


if __name__ == "__main__":
    # Influenza incubation period ~ 2 days (according to the literature),
    #  so sigma = 1/2
    sigma = 1 / 2
    
    # Initial guess for beta, gamma is taken from the SIR result from part 1
    initial_guess = (0.420, 0.265)
    
    # Now we are picking one full flu season in HHS Region 4 (Oct 2018 - May 2019)
    print("Loading influenza data...", flush=True)
    df = import_influenza_data(region=4, start_wk=40, start_yr=2018, end_wk=20, end_yr=2019)
    
    # Compute N from peak data (for flu-like outbreak, peak ~ 0.15 * N, so N ~ 10 * peak)
    new_cases_raw = df[NEW_INFECTIONS_COL].values
    baseline = new_cases_raw[0]
    peak_data = np.clip(new_cases_raw - baseline, 0, None).max()
    N = 10 * peak_data
    print(f"Data peak (excess): {peak_data:.0f}, chosen N: {N:.0f}", flush=True)
    
    # Initial conditions: epidemic just starting (small seed of infections)
    I_init = 100.0
    E_init = 100.0
    R_init = 0.0
    S_init = N - I_init - E_init - R_init
    
    fit_seir_nelder_mead(
        df=df,
        date_col=DATE_COL,
        observed_col=NEW_INFECTIONS_COL,
        sigma=sigma,
        initial_guess=initial_guess,
        S_init=S_init,
        E_init=E_init,
        I_init=I_init,
        R_init=R_init,
        dt=0.5,
        output_path="diagrams/seir_fit_influenza.png",
        title="SEIR Fit to Influenza Data (2018-2019, Region 4)"
    )
