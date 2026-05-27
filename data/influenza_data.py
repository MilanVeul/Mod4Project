import pandas as pd

DATE_COL = "DATE"
WEEK_COL = "WEEK"
YEAR_COL = "YEAR"
NEW_INFECTIONS_COL = "New Infections"
TOTAL_INFECTIONS_COL = "Total Infections"

def import_influenza_data(region) -> pd.DataFrame:
    cols = ["YEAR", "WEEK", "% WEIGHTED ILI", "ILITOTAL", "TOTAL PATIENTS", "REGION"]

    df = pd.read_csv("data/ILINet.csv", usecols=cols)
    region_df = df[df["REGION"] == f"Region {region}"]
    region_df = region_df.sort_values(by=["YEAR", "WEEK"]).reset_index(
        drop=True
    )
    # df[DATE_COL] = f"{df["YEAR"]-df["WEEK"]}"
    df[DATE_COL] = df["YEAR"] * 52 + df["WEEK"]
    df[NEW_INFECTIONS_COL] = df["% WEIGHTED ILI"] * df["TOTAL PATIENTS"]
    df[TOTAL_INFECTIONS_COL] = df[NEW_INFECTIONS_COL].cumsum()
    return df
    


