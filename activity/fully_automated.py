# Required: pip install pandas rdkit autogluon 'setuptools==81.0.0' scikit-learn mordredcommunity autogluon.tabular[tabpfn] tabdpt
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
from sklearn.ensemble import RandomForestRegressor
from autogluon.tabular import TabularPredictor

# 1. Mordred Feature Extraction
def extract_mordred_features(smiles_series):
    print(f"Calculating Mordred descriptors for {len(smiles_series)} molecules...")
    # Initialize Mordred Calculator with all descriptors (2D and 3D)
    calc = Calculator(descriptors, ignore_3D=False)
    
    mols = []
    for s in smiles_series:
        mol = Chem.MolFromSmiles(s)
        mol = Chem.AddHs(mol)
        mols.append(mol)
    
    # Calculate and convert to DataFrame
    df_mordred = calc.pandas(mols, quiet=False)
    
    # Mordred returns custom Error objects for failed descriptors; convert to NaN
    df_mordred = df_mordred.apply(pd.to_numeric, errors='coerce')
    # Fill NaNs with 0 for the Random Forest selector
    df_mordred = df_mordred.fillna(0)
    
    return df_mordred

# 2. Feature Selection via Random Forest
def select_top_features(X, y, n_top=128):
    print(f"Selecting top {n_top} features using Random Forest...")
    # Using a fast RF to rank importance
    rf = RandomForestRegressor(n_estimators=100, n_jobs=-1, random_state=42)
    rf.fit(X, y)
    
    importances = pd.Series(rf.feature_importances_, index=X.columns)
    top_features = importances.sort_values(ascending=False).head(n_top).index.tolist()
    return top_features

# 3. Pipeline logic
def prepare_and_train(train_path):
    df = pd.read_csv(train_path)
    
    # Extract all ~1800 Mordred features
    mordred_df = extract_mordred_features(df['SMILES'])
    
    # Select 128 most important
    # Note: Using pEC50 for importance ranking
    target = df['pEC50']
    top_cols = select_top_features(mordred_df, target, n_top=128)
    
    # Build final training set: Top Descriptors + Target + Weights
    final_train = pd.concat([
        df[['pEC50', 'pEC50_weight']], 
        mordred_df[top_cols]
    ], axis=1)

    print("Starting AutoGluon Training...")
    predictor = TabularPredictor(
        label='pEC50', 
        problem_type='regression',
        eval_metric='mean_absolute_error',
        sample_weight="pEC50_weight",
    ).fit(
        final_train,
        presets='extreme_quality',  # GPU required, trains tabular foundation models
        time_limit=600 * 2,  # 20 minutes for training
    )
    
    return predictor, top_cols

def predict_pEC50(predictor, selected_features, smiles_list):
    """
    Inference needs to calculate the SAME Mordred features selected during training.
    """
    # Calculate Mordred for new SMILES
    mordred_df = extract_mordred_features(pd.Series(smiles_list))
    
    # Filter to only the 128 features the model expects
    inference_df = mordred_df[selected_features].copy()
    inference_df['SMILES'] = smiles_list
    
    return predictor.predict(inference_df)

if __name__ == "__main__":
    # 1. Train and identify key features
    predictor, top_features = prepare_and_train("train.csv")
    
    # 2. Run inference using the same feature list
    test_df = pd.read_csv("test.csv")
    preds = predict_pEC50(predictor, top_features, test_df["SMILES"].tolist())
    
    # 3. Save results
    test_df["pEC50"] = preds
    test_df[["Molecule Name", "SMILES", "pEC50"]].to_csv("submission.csv", index=False)
    print("Inference complete. Results saved to submission.csv")
