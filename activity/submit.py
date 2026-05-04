import pandas as pd

if __name__ == "__main__":
    # Load the predictions
    df = pd.read_csv("train_output_denoised/predictions.csv")
    
    # undo augmentation
    pred = df.groupby("Molecule Name").median(numeric_only=True)["pEC50"]

    df_test = pd.read_csv("test.csv")
    df_test["pEC50"] = pred.reindex(df_test["Molecule Name"]).values
    df_test.to_csv("submission.csv", index=False)
