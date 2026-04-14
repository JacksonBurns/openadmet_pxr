import pandas as pd
from molpipeline.predefined_pipelines import get_rf_regressor_baseline

if __name__ == "__main__":
    train_df = pd.read_csv("train_augmented.csv")
    test_df = pd.read_csv("test.csv")
    rf = get_rf_regressor_baseline(n_jobs=-1, random_state=42, error_handling=True)
    rf.fit(train_df["SMILES"], train_df["pEC50"])
    test_pred = rf.predict(test_df["SMILES"])
    test_df["rf_pred"] = test_pred
    test_df.to_csv("train_output/rf_predictions.csv", index=False)
