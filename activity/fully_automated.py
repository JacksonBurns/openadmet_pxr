# Required: pip install pandas autogluon 'setuptools==81.0.0' autogluon.tabular[tabpfn] tabdpt chemprop torch
import pandas as pd
from autogluon.tabular import TabularPredictor
import torch

from chemeleon_fingerprint import CheMeleonFingerprint


def extract_chemeleon_features(smiles_series):
    chemeleon_fingerprint = CheMeleonFingerprint(device="cuda" if torch.cuda.is_available() else "cpu")
    return pd.DataFrame(data=chemeleon_fingerprint(smiles_series), columns=[f"chemeleon_{i}" for i in range(chemeleon_fingerprint.model.message_passing.output_dim)])

def prepare_and_train(train_path):
    df = pd.read_csv(train_path)
    train_df = extract_chemeleon_features(df["SMILES"])
    train_df['pEC50'] = df['pEC50']
    train_df['pEC50_weight'] = df["pEC50_weight"]
    predictor = TabularPredictor(
        label='pEC50', 
        problem_type='regression',
        eval_metric='mean_absolute_error',
        sample_weight="pEC50_weight",
    ).fit(
        train_df,
        presets='best_quality',  # extreme_quality -> GPU required, trains tabular foundation models, limited to 500 features
        time_limit=60 * 60,  # 1 hour time limit for training
    )
    
    return predictor

def predict_pEC50(predictor, smiles_list):
    feature_df = extract_chemeleon_features(smiles_list)
    return predictor.predict(feature_df)

if __name__ == "__main__":
    predictor = prepare_and_train("train.csv")
    test_df = pd.read_csv("test.csv")
    preds = predict_pEC50(predictor, test_df["SMILES"].tolist())
    test_df["pEC50"] = preds
    test_df[["Molecule Name", "SMILES", "pEC50"]].to_csv("submission.csv", index=False)
    print("Inference complete. Results saved to submission.csv")
