import pandas as pd

DATE_COL = "Date"
NEW_INFECTIONS_COL = "New Infections"
NEW_RECOVERIES_COL = "New Recoveries"
TOTAL_INFECTIONS_COL = "Total Infections"
TOTAL_RECOVERIES_COL = "Total Recoveries"
INFECTED_COL = "Infected"

def import_data() -> pd.DataFrame:
    df = pd.read_csv("Meridonia_infections_data.csv")

    df[DATE_COL] = pd.to_datetime(df[DATE_COL])
    df[TOTAL_INFECTIONS_COL] = df[NEW_INFECTIONS_COL].cumsum()
    df[TOTAL_RECOVERIES_COL] = df[NEW_RECOVERIES_COL].cumsum()
    df[INFECTED_COL] = df[TOTAL_INFECTIONS_COL] - df[TOTAL_RECOVERIES_COL]

    return df
