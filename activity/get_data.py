if __name__ == "__main__":
    import pandas as pd

    train_df         = pd.read_csv("hf://datasets/openadmet/pxr-challenge-train-test/pxr-challenge_TRAIN.csv")
    test_df          = pd.read_csv("hf://datasets/openadmet/pxr-challenge-train-test/pxr-challenge_TEST_BLINDED.csv")
    train_counter_df = pd.read_csv("hf://datasets/openadmet/pxr-challenge-train-test/pxr-challenge_counter-assay_TRAIN.csv")
    train_single_df  = pd.read_csv("hf://datasets/openadmet/pxr-challenge-train-test/pxr-challenge_single_concentration_TRAIN.csv")
    train_df = pd.merge(train_df, train_counter_df, on="SMILES", how="outer", suffixes=("", "_counter"))
    train_df = train_df[['Molecule Name', 'SMILES', 'pEC50', 'pEC50_ci.lower (-log10(molarity))',  'pEC50_ci.upper (-log10(molarity))',  'Emax_estimate (log2FC vs. baseline)', 'Emax_ci.lower (log2FC vs. baseline)', 'Emax_ci.upper (log2FC vs. baseline)',        'pEC50_std.error (-log10(molarity))',   'Emax_std.error (log2FC vs. baseline)', 'pEC50_counter',       'pEC50_ci.lower (-log10(molarity))_counter',    'pEC50_ci.upper (-log10(molarity))_counter',       'Emax_estimate (log2FC vs. baseline)_counter',    'Emax_ci.lower (log2FC vs. baseline)_counter',       'Emax_ci.upper (log2FC vs. baseline)_counter',    'pEC50_std.error (-log10(molarity))_counter', 'Emax_std.error (log2FC vs. baseline)_counter']]
    train_df["pEC50_weight"] = 1 / train_df["pEC50_std.error (-log10(molarity))"]
    q_min, q_max = train_df["pEC50_weight"].quantile([0.05, 0.95])
    train_df["pEC50_weight"] = train_df["pEC50_weight"].clip(lower=q_min, upper=q_max)
    train_df.to_csv("train.csv", index=False)
    test_df.to_csv("test.csv", index=False)
