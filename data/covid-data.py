

# Date format: "DD-MM-YY"
def import_covid_data(country, start_date, end_date) -> pd.DataFrame:


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
