import os

import pandas as pd

DATA_PATH = "csse_covid_19_data/csse_covid_19_daily_reports/"

def import_covid_data(country="Italy", start_date="2020-01-22", end_date="2020-07-31") -> pd.DataFrame:
    dates = pd.date_range(start=start_date, end=end_date, freq="D")
    df_list = []
    combined_df = pd.DataFrame()

    for dt in dates:
        filename = f"{dt.strftime('%m-%d-%Y')}.csv"
        file_path = os.path.join(DATA_PATH, filename)


        if not os.path.exists(file_path):
            print(f"Warning: {filename} not found.")
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
            country_df["Date"] = date_str

            COLS_TO_KEEP = ["Confirmed", "Deaths", "Recovered"]
            available_metrics = [
                col for col in COLS_TO_KEEP if col in country_df.columns
            ]

            # add all states if there are multiple rows for one country
            country_df = country_df.groupby(
                ["Date", country_col], as_index=False
            )[available_metrics].sum()

            # remove all columns except useful ones
            COLS_TO_KEEP = ["Date", "Confirmed", "Deaths", "Recovered"]
            country_df = country_df[
                [col for col in COLS_TO_KEEP if col in country_df.columns]
            ]

            # add removed col
            country_df["Removed"] = country_df["Recovered"] + country_df["Deaths"]

            df_list.append(country_df)

    if df_list:
        combined_df = pd.concat(df_list, ignore_index=True)
        print(f"Successfully combined {country} data! Total rows: {len(combined_df)}")
        print(combined_df)
        combined_df.to_csv(f"{country.lower()}_covid_data.csv", index=False)
    else:
        print("No matching data found.")
    return combined_df


if __name__ == "__main__":
    df = import_covid_data()