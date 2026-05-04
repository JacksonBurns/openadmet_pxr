import pandas as pd

if __name__ == "__main__":
    # load model predictions
    df_pred = pd.read_csv("train_output/training_predictions.csv")

    # load original training data
    df_train = pd.read_csv("train_augmented.csv")

    # undo augmentation
    df_train = df_train.groupby("Molecule Name").first()
    df_pred = df_pred.groupby("Molecule Name").median(numeric_only=True)

    # calculate error
    df_error = df_pred.sub(df_train["pEC50"], axis=0)

    # drop top 10% of largest errors
    df_error = df_error.drop(df_error.nlargest(int(len(df_error) * 0.1), columns="pEC50").index)

    # save only those molecules
    df_train = pd.read_csv("train.csv")
    df_denoised = df_train[df_train["Molecule Name"].isin(df_error.index)]
    df_denoised.to_csv("train_denoised.csv", index=False)
