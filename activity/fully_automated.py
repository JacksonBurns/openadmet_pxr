# Required: pip install pandas autogluon 'setuptools==81.0.0' autogluon.tabular[tabpfn] tabdpt chemprop torch scikit-learn
import pandas as pd
from autogluon.tabular import TabularPredictor
import torch
from sklearn.decomposition import PCA
from chemeleon_fingerprint import CheMeleonFingerprint


if __name__ == "__main__":
    chemeleon_fingerprint = CheMeleonFingerprint(device="cuda" if torch.cuda.is_available() else "cpu")
    num_features = 500  # max for tabpfn
    df = pd.read_csv("train.csv")
    features = chemeleon_fingerprint(df["SMILES"].tolist())
    pca = PCA(n_components=num_features)
    reduced_features = pca.fit_transform(features)
    train_df = pd.DataFrame(data=reduced_features, columns=[f"chemeleon_{i}" for i in range(num_features)])
    train_df['pEC50'] = df['pEC50']
    train_df['pEC50_weight'] = df["pEC50_weight"]
    predictor = TabularPredictor(
        label='pEC50', 
        problem_type='regression',
        eval_metric='mean_absolute_error',
        sample_weight="pEC50_weight",
    ).fit(
        train_df,
        # extreme_quality -> GPU required, trains tabular foundation models, limited to 500 features
        # best_quality -> CPU only, trains classical models, can use any number of features
        presets='extreme_quality',
        time_limit=60 * 60,  # 1 hour time limit for training
    )

    test_df = pd.read_csv("test.csv")
    features = chemeleon_fingerprint(test_df["SMILES"].tolist())
    reduced_features = pca.transform(features)
    feature_df = pd.DataFrame(data=reduced_features, columns=[f"chemeleon_{i}" for i in range(num_features)])
    preds = predictor.predict(feature_df)
    test_df["pEC50"] = preds
    test_df[["Molecule Name", "SMILES", "pEC50"]].to_csv("submission.csv", index=False)
    print("Inference complete. Results saved to submission.csv")
