from common import import_data, NEW_INFECTIONS_COL, NEW_RECOVERIES_COL, TOTAL_INFECTIONS_COL, TOTAL_RECOVERIES_COL, INFECTED_COL
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = import_data()

plt.figure(figsize=(10, 6))
plt.plot(df["Date"], df[NEW_INFECTIONS_COL], marker="o", linestyle="-")
plt.title("Daily new infections in Meridonia")
plt.xlabel("Date")
plt.ylabel("Number of new infections")
plt.xticks(rotation=45)
plt.grid()
plt.savefig("diagrams/new_infections.png")

plt.figure(figsize=(10, 6))
plt.plot(df["Date"], df[NEW_RECOVERIES_COL], marker="o", linestyle="-")
plt.title("Daily new recoveries in Meridonia")
plt.xlabel("Date")
plt.ylabel("Number of new recoveries")
plt.xticks(rotation=45)
plt.grid()
plt.savefig("diagrams/new_recoveries.png")

plt.figure(figsize=(10, 6))
plt.plot(df["Date"], df[NEW_INFECTIONS_COL], marker="o", linestyle="-")
plt.plot(df["Date"], df[NEW_RECOVERIES_COL], marker="o", linestyle="-")
plt.title("Daily new infections and recoveries in Meridonia")
plt.xlabel("Date")
plt.ylabel("Number of cases")
plt.xticks(rotation=45)
plt.grid()
plt.legend(["New infections", "New recoveries"])
plt.savefig("diagrams/new_infections_and_recoveries.png")

plt.figure(figsize=(10, 6))
plt.plot(df["Date"], df[TOTAL_INFECTIONS_COL], marker="o", linestyle="-")
plt.plot(df["Date"], df[TOTAL_RECOVERIES_COL], marker="o", linestyle="-")
plt.plot(df["Date"], df[INFECTED_COL], marker="o", linestyle="-")
plt.title("Total infections and recoveries in Meridonia")
plt.xlabel("Date")
plt.ylabel("Number of cases")
plt.xticks(rotation=45)
plt.grid()
plt.legend(["Total infections", "Total recoveries", "Infected"])
plt.savefig("diagrams/total_infections_and_recoveries.png")