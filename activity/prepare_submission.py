import argparse

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error


def parity_plot(y_true, y_pred, outname):
    """Create a parity plot with hexbin-based density."""
    # Compute regression statistics
    r, _ = pearsonr(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = mean_absolute_error(y_true, y_pred)

    # Create plot
    plt.figure(figsize=(6, 6))

    hb = plt.hexbin(y_true, y_pred, gridsize=50, mincnt=1,)
    cb = plt.colorbar(
        hb,
        fraction=0.04,  # width relative to axes
        pad=0.02,  # gap between plot and colorbar
        shrink=0.85,  # height scaling
    )
    cb.set_label("# of Compounds", fontsize=11)

    # 1:1 line
    lims = [
        min(y_true.min(), y_pred.min()),
        max(y_true.max(), y_pred.max()),
    ]
    plt.plot(lims, lims, color="red", linestyle="--", linewidth=1.5)

    # Labels and title
    plt.xlabel("Measured log(solubility) [mol/L]", fontsize=12)
    plt.ylabel("Predicted log(solubility) [mol/L]", fontsize=12)
    plt.title("Parity Plot", fontsize=14)

    # Annotation with statistics
    stats_text = f"$r$ = {r:.3f}\n" f"RMSE = {rmse:.3f}\n" f"MAE = {mae:.3f}"
    plt.text(
        0.05,
        0.95,
        stats_text,
        transform=plt.gca().transAxes,
        fontsize=11,
        verticalalignment="top",
        horizontalalignment="left",
        bbox=dict(facecolor="white", edgecolor="gray", boxstyle="round,pad=0.3"),
    )

    # Aesthetics
    plt.gca().set_aspect("equal", adjustable="box")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(outname, dpi=300)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", help="Output CSV file")
    args = parser.parse_args()

    # assemble cv predictions and train a meta-model on the cv predictions, then use that meta-model to predict on the test set
    # taking care to also fold the exploded resonance structures back together by taking the median prediction for each molecule
    # note that chemprop-based models have already averaged across the 4 inner kmeans folds
    cv_df = pd.read_csv("train_augmented.csv")
    holdout_df = pd.read_csv("test_augmented.csv")
    for fold_number in range(5):
        subdir = f"fold_{fold_number}"
        rf_df = pd.read_csv(f"train_output/rf_{subdir}_cv_predictions.csv")
        chemeleon_df = pd.read_csv(f"train_output/chemeleon_{subdir}_cv_predictions.csv")
        chemprop_df = pd.read_csv(f"train_output/chemprop_{subdir}_cv_predictions.csv")
        # add predictions to cv_df
        cv_df.loc[cv_df["Molecule Name"].isin(rf_df["Molecule Name"]), "rf_pred"] = rf_df["rf_pred"].values
        cv_df.loc[cv_df["Molecule Name"].isin(chemeleon_df["Molecule Name"]), "chemeleon_pred"] = chemeleon_df["pEC50"].values
        cv_df.loc[cv_df["Molecule Name"].isin(chemprop_df["Molecule Name"]), "chemprop_pred"] = chemprop_df["pEC50"].values

       
        # repeat for holdout predictions
        rf_holdout_df = pd.read_csv(f"train_output/rf_{subdir}_holdout_predictions.csv")
        chemeleon_holdout_df = pd.read_csv(f"train_output/chemeleon_{subdir}_holdout_predictions.csv")
        chemprop_holdout_df = pd.read_csv(f"train_output/chemprop_{subdir}_holdout_predictions.csv")
        holdout_df.loc[holdout_df["Molecule Name"].isin(rf_holdout_df["Molecule Name"]), "rf_pred"] = rf_holdout_df["rf_pred"].values
        holdout_df.loc[holdout_df["Molecule Name"].isin(chemeleon_holdout_df["Molecule Name"]), "chemeleon_pred"] = chemeleon_holdout_df["pEC50"].values
        holdout_df.loc[holdout_df["Molecule Name"].isin(chemprop_holdout_df["Molecule Name"]), "chemprop_pred"] = chemprop_holdout_df["pEC50"].values
    
    # fold resonance structures back together by taking the median prediction for each molecule
    cv_df = cv_df.groupby("Molecule Name").median(numeric_only=True).reset_index()
    holdout_df = holdout_df.groupby("Molecule Name").median(numeric_only=True).reset_index()

    # random forest has one nan prediction - just fill with median of other predictions for that molecule
    for df in [cv_df, holdout_df]:
        nan_rows = df[df["rf_pred"].isna()]
        for i, row in nan_rows.iterrows():
            molecule_name = row["Molecule Name"]
            median_pred = df[df["Molecule Name"] == molecule_name][["chemeleon_pred", "chemprop_pred"]].median().median()
            df.at[i, "rf_pred"] = median_pred

    meta_model = RandomForestRegressor(n_estimators=100, random_state=42, max_depth=5)
    meta_model.fit(cv_df[["rf_pred", "chemeleon_pred", "chemprop_pred"]], cv_df["pEC50"])
    cv_df["meta_pred"] = meta_model.predict(cv_df[["rf_pred", "chemeleon_pred", "chemprop_pred"]])  # overfit, but for inspection purposes only
    holdout_pred = meta_model.predict(holdout_df[["rf_pred", "chemeleon_pred", "chemprop_pred"]])
    # massage to expected format: Molecule Name,SMILES,pEC50
    smiles = pd.read_csv("test.csv")["SMILES"]
    pd.DataFrame({"Molecule Name": holdout_df["Molecule Name"], "SMILES": smiles, "pEC50": holdout_pred}).to_csv(args.output, index=False)

    # make parity plot of cv predictions, including the meta model
    parity_plot(cv_df["pEC50"], cv_df["rf_pred"], "train_output/rf_cv_parity_plot.png")
    parity_plot(cv_df["pEC50"], cv_df["chemeleon_pred"], "train_output/chemeleon_cv_parity_plot.png")
    parity_plot(cv_df["pEC50"], cv_df["chemprop_pred"], "train_output/chemprop_cv_parity_plot.png")
    parity_plot(cv_df["pEC50"], cv_df["meta_pred"], "train_output/meta_model_cv_parity_plot.png")
