import pandas as pd
from epiweeks import Week

DATE_COL = "DATE"
WEEK_COL = "WEEK"
YEAR_COL = "YEAR"
NEW_INFECTIONS_COL = "New Infections"
TOTAL_INFECTIONS_COL = "Total Infections"

def import_influenza_data(region) -> pd.DataFrame:
    cols = ["YEAR", "WEEK", "% WEIGHTED ILI", "ILITOTAL", "TOTAL PATIENTS", "REGION"]

    df = pd.read_csv("data/ILINet.csv", usecols=cols)
    df = df[df["REGION"] == f"Region {region}"]
    df = df.sort_values(by=["YEAR", "WEEK"]).reset_index(
        drop=True
    )
    # Date column
    # iso_week_string = region_df["YEAR"].astype(str) + "-W" + region_df["WEEK"].astype(str) + "-1"
    # df[DATE_COL] = pd.to_datetime(iso_week_string, format="%G-W%V-%u")

    df[DATE_COL] = df.apply(
        lambda row: pd.to_datetime(
            Week(int(row["YEAR"]), int(row["WEEK"])).startdate()
        ),
        axis=1,
    )

    df[NEW_INFECTIONS_COL] = df["% WEIGHTED ILI"] * df["TOTAL PATIENTS"]
    df[TOTAL_INFECTIONS_COL] = df[NEW_INFECTIONS_COL].cumsum()
    return df
    


