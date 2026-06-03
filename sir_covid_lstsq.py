from sir_meridonia_lstsq import evaluate_sir
from data.covid_data import import_covid_data, TOTAL_RECOVERIES_COL, INFECTED_COL, DATE_COL
from diff import derivative
import numpy as np
import matplotlib.pyplot as plt

sigma = 0.2
# N = 60_000_000
N = 1_000_000

df = import_covid_data(country="Italy", start_date="2020-02-25")
I0 = df[INFECTED_COL].iloc[0]
R0 = df[TOTAL_RECOVERIES_COL].iloc[0]
S0 = N - I0 - R0

R = df[TOTAL_RECOVERIES_COL]
I = df[INFECTED_COL]
S = N - I - R

evaluate_sir(df, S, I, R, N, S0, I0, R0, 0.1)