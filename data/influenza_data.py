import pandas as pd
from epiweeks import Week

DATE_COL = "DATE"
WEEK_COL = "WEEK"
YEAR_COL = "YEAR"
NEW_INFECTIONS_COL = "New Infections"
TOTAL_INFECTIONS_COL = "Total Infections"

def import_influenza_data(region, start_wk, start_yr, end_wk, end_yr) -> pd.DataFrame:
    cols = ["YEAR", "WEEK", "% WEIGHTED ILI", "ILITOTAL", "TOTAL PATIENTS", "REGION"]

    df = pd.read_csv("data/ILINet.csv", usecols=cols)
    df = df[df["REGION"] == f"Region {region}"]
    df = df.sort_values(by=["YEAR", "WEEK"]).reset_index(
        drop=True
    )

    df["EPI_SCORE"] = df["YEAR"] * 100 + df["WEEK"]
    start_score = start_yr * 100 + start_wk
    end_score = end_yr * 100 + end_wk
    df = df[
        (df["EPI_SCORE"] >= start_score)
        & (df["EPI_SCORE"] <= end_score)
    ]

    df[DATE_COL] = df.apply(
        lambda row: pd.to_datetime(
            Week(int(row["YEAR"]), int(row["WEEK"])).startdate()
        ),
        axis=1,
    )

    df[NEW_INFECTIONS_COL] = df["% WEIGHTED ILI"] * df["TOTAL PATIENTS"]
    df[TOTAL_INFECTIONS_COL] = df[NEW_INFECTIONS_COL].cumsum()
    return df
    


