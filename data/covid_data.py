import os

import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent

DATE_COL = "Date"
NEW_INFECTIONS_COL = "New Infections"
NEW_RECOVERIES_COL = "New Recoveries"
TOTAL_INFECTIONS_COL = "Total Infections"
TOTAL_RECOVERIES_COL = "Total Recoveries"
INFECTED_COL = "Infected"
DATA_PATH = SCRIPT_DIR / "csse_covid_19_data/csse_covid_19_daily_reports/"

def import_covid_data(country="Italy", start_date="2020-01-22", end_date="2020-07-31") -> pd.DataFrame:
    dates = pd.date_range(start=start_date, end=end_date, freq="D")
    df_list = []
    combined_df = pd.DataFrame()

    for dt in dates:
        filename = f"{dt.strftime('%m-%d-%Y')}.csv"
        file_path = DATA_PATH / filename

        if not file_path.exists():
            print(f"Warning: {file_path} not found.")
            continue

        df = pd.read_csv(file_path)

        if "Country/Region" in df.columns:
            country_col = "Country/Region"
        elif "Country_Region" in df.columns:
            country_col = "Country_Region"
        else:
            print(f"Warning: no country column found in {filename}.")
            continue

        if "Province/State" in df.columns:
            state_col = "Province/State"
        elif "Province_State" in df.columns:
            state_col = "Province_State"
        else:
            print(f"Warning: no state column found in {filename}.")
            continue

        country_df = df[df[country_col] == country]
        if not country_df.empty:
            date_str = dt.strftime("%Y-%m-%d")
            country_df = country_df.copy()
            country_df[DATE_COL] = date_str

            # add all states if there are multiple rows for one country
            METRIC_COLS = ["Confirmed", "Deaths", "Recovered"]
            available_metrics = [col for col in METRIC_COLS if col in country_df.columns]
            country_df = country_df.groupby(
                [DATE_COL, country_col], as_index=False
            )[available_metrics].sum()

            # add column types
            country_df[TOTAL_INFECTIONS_COL] = country_df["Confirmed"].astype(int)
            country_df[TOTAL_RECOVERIES_COL] = country_df["Recovered"].astype(int) + country_df["Deaths"].astype(int)
            country_df[INFECTED_COL] = country_df[TOTAL_INFECTIONS_COL] - country_df[TOTAL_RECOVERIES_COL]

            # remove all columns except useful ones
            COLS_TO_KEEP = [DATE_COL, TOTAL_INFECTIONS_COL, TOTAL_RECOVERIES_COL, INFECTED_COL]
            country_df = country_df[
                [col for col in COLS_TO_KEEP if col in country_df.columns]
            ]

            df_list.append(country_df)

    if df_list:
        combined_df = pd.concat(df_list, ignore_index=True)
        combined_df[NEW_INFECTIONS_COL] = combined_df[TOTAL_INFECTIONS_COL].diff().fillna(0).astype(int)
        combined_df[NEW_RECOVERIES_COL] = combined_df[TOTAL_RECOVERIES_COL].diff().fillna(0).astype(int)

        print(f"Successfully combined {country} data! Total rows: {len(combined_df)}")
        combined_df.to_csv(str(SCRIPT_DIR / f"{country.lower()}_covid_data.csv"), index=False)
    else:
        print("No matching data found.")

    # Fix dates
    combined_df[DATE_COL] = pd.to_datetime(combined_df[DATE_COL])
    return combined_df


if __name__ == "__main__":
    df = import_covid_data()