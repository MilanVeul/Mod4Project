from data.influenza_data import import_influenza_data, NEW_INFECTIONS_COL, TOTAL_INFECTIONS_COL, WEEK_COL, YEAR_COL, DATE_COL
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = import_influenza_data(region=1)

plt.figure(figsize=(10, 6))
plt.plot(df[DATE_COL], df[NEW_INFECTIONS_COL], marker="", linestyle="-")
plt.title("Daily new influenza")
plt.xlabel("Date")
plt.ylabel("Number of new infections")
# plt.xticks(rotation=45)
plt.grid()
plt.savefig("diagrams/influenza/new_infections.png")