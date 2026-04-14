import argparse

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", help="Output CSV file")
    args = parser.parse_args()

    rf_df = pd.read_csv("train_output/rf_predictions.csv")
    chemeleon_df = pd.read_csv("train_output/chemeleon_predictions.csv")
    chemprop_df = pd.read_csv("train_output/chemprop_predictions.csv")
    pred = (2*rf_df["rf_pred"] + chemeleon_df["pEC50"] + chemprop_df["pEC50"]) / 4  # weighted average
    test_df = pd.read_csv("test.csv")
    test_df["pEC50"] = pred
    test_df.to_csv(args.output, index=False)
