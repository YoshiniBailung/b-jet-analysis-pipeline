import pandas as pd

def load_parquet(path):
    print(f"Reading file {path}")
    return pd.read_parquet(path)

def cleanup(df):
    num_cols = df.select_dtypes(include="number").columns

    low = df[num_cols].quantile(0.001)
    high = df[num_cols].quantile(0.999)

    mask = ((df[num_cols] >= low) & (df[num_cols] <= high)).all(axis=1)

    df_clean = df[mask]

    outliers = df[~mask]

    # df_clean = df[~(df == -99999.).any(axis=1)]

    # print("Cleaned")

    # outliers = df[(df == -99999.).any(axis=1)] 
    # print(outliers.head())
    n_outliers = len(outliers)
    frac_outliers = n_outliers/len(df)
    print(f"Outliers: {n_outliers} ({frac_outliers:.2%})")

    return df_clean

    