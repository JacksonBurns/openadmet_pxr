import pandas as pd
from molpipeline.predefined_pipelines import get_rf_regressor_baseline

if __name__ == "__main__":
    for fold_number in range(5):
        subdir = f"fold_{fold_number}"
        train_df = pd.read_csv(f"splits/{subdir}/train_val.csv")
        test_df = pd.read_csv(f"splits/{subdir}/test.csv")
        rf = get_rf_regressor_baseline(n_jobs=-1, random_state=42, error_handling=True)
        rf.fit(train_df["SMILES"], train_df["pEC50"])
        test_pred = rf.predict(test_df["SMILES"])
        test_df["rf_pred"] = test_pred
        test_df.to_csv(f"train_output/rf_{subdir}_cv_predictions.csv", index=False)
        holdout_df = pd.read_csv("test_augmented.csv")
        holdout_pred = rf.predict(holdout_df["SMILES"])
        holdout_df["rf_pred"] = holdout_pred
        holdout_df.to_csv(f"train_output/rf_{subdir}_holdout_predictions.csv", index=False)
