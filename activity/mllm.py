import re
import logging
from typing import Literal

import numpy as np
import pandas as pd
import hdbscan
from pysr import PySRRegressor
from mordred import Calculator, descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import MolFromSmiles, Descriptors
from sklearn.preprocessing import StandardScaler

# Configure logging
logging.basicConfig(
    filename='pysr_pipeline.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filemode='w' # Overwrites the log file each run; use 'a' to append
)
logger = logging.getLogger(__name__)

def _add_features(df: pd.DataFrame, smiles_col: str = "SMILES", feature_set: Literal["rdkit", "mordred"] = "rdkit", means: np.ndarray | None = None):
    if feature_set == "mordred":
        calc = Calculator(descriptors, ignore_3D=True)
        descs = calc.pandas(mols=[MolFromSmiles(s) for s in df[smiles_col]]).fill_missing()
        # retain only Only alphanumeric characters, numbers, and underscores in column names for compatibility with PySR
        descs.columns = [re.sub(r'[^\w]+', '_', col) for col in descs.columns]
        # suffix with _mordred to avoid conflicts with julia vars
        descs.columns = [col + "_mordred" for col in descs.columns]
    else:
        names = [x[0] for x in Descriptors._descList]
        calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)
        data = [calc.CalcDescriptors(MolFromSmiles(smiles)) for smiles in df[smiles_col]]
        descs = pd.DataFrame(columns=calc.GetDescriptorNames(), data=data)

    descs = descs.replace([np.inf, -np.inf], np.nan)
    if means is None:
        # replace columns which are entirely nan with 0
        descs[descs.columns[descs.isna().all()]] = 0
        means = descs.mean(axis=0, skipna=True)
    descs = descs.fillna(means)
    return descs, means


if __name__ == "__main__":
    logger.info("Loading training data...")
    df = pd.read_csv("train.csv")
    
    logger.info("Extracting features...")
    train_features, means = _add_features(df, smiles_col="SMILES")
    
    # Upcast to float32 for training stability
    train_features_f32 = train_features.astype(np.float32)
    
    logger.info("Scaling and clustering with HDBSCAN...")
    scaler = StandardScaler()
    train_features_scaled = scaler.fit_transform(train_features_f32)
    
    # HDBSCAN clustering: prediction_data=True is required to use approximate_predict later
    # min_cluster_size dictates the smallest allowable cluster size.
    clusterer = hdbscan.HDBSCAN(min_cluster_size=15, prediction_data=True)
    train_labels = clusterer.fit_predict(train_features_scaled)
    df['cluster'] = train_labels
    
    logger.info(f"Discovered {len(set(train_labels))} clusters (including noise label -1).")
    
    models = {}
    for cluster_id in df['cluster'].unique():
        cluster_size = sum(df['cluster'] == cluster_id)
        logger.info(f"Fitting PySR model for cluster {cluster_id} (n={cluster_size})...")
        
        cluster_df = df[df['cluster'] == cluster_id]
        cluster_features = train_features_f32.loc[cluster_df.index]
        
        model = PySRRegressor(
            niterations=50,  
            populations=30,  
            population_size=100,  
            binary_operators=["+", "-", "*", "/"],
            loss="L2DistLoss()",  
            parsimony=0.002,  
            maxsize=20,  
            maxdepth=10,  
            turbo=True,  
            bumper=True,  
            precision=64,  
            random_state=42,
            parallelism='serial',
            deterministic=True,
            progress=False,  # Turn off progress bar to keep console/logs clean
            verbosity=0,
            temp_equation_file=True,
        )
        # Upcast the target variables to float32 as well 
        model.fit(cluster_features, cluster_df['pEC50'].astype(np.float32))
        
        # Store the equation representation
        models[cluster_id] = str(model.sympy())
        logger.info(f"Cluster {cluster_id} equation: {models[cluster_id]}")
        
    logger.info("Processing test data...")
    test_df = pd.read_csv("test.csv")
    test_features, _ = _add_features(test_df, smiles_col="SMILES", means=means)
    
    # Upcast to float32 for forward pass predictions
    test_features_f32 = test_features.astype(np.float32)
    test_features_scaled = scaler.transform(test_features_f32)
    
    logger.info("Assigning clusters to test data...")
    # approximate_predict uses the density tree to assign clusters or -1 for noise
    test_labels, probabilities = hdbscan.approximate_predict(clusterer, test_features_scaled)
    test_df['cluster'] = test_labels
    
    logger.info("Evaluating local equations for test data...")
    preds = []
    # Evaluate equation per row depending on its cluster assignment
    for idx, row in test_df.iterrows():
        cluster_id = row['cluster']
        eqn = models.get(cluster_id)
        
        # Fallback if somehow a test point gets a cluster id not seen in training
        if eqn is None:
            logger.error(f"No equation found for cluster {cluster_id} at row {idx}. Exiting for debugging.")
            exit(1)
            
        row_features = test_features_f32.iloc[[idx]]
        try:
            pred = row_features.eval(eqn).values[0]
        except Exception as e:
            logger.error(f"Failed to evaluate equation for row {idx} in cluster {cluster_id}. Error: {e}. Exiting for debugging.")
            exit(1)
        preds.append(pred)
        
    test_df["pEC50"] = preds
    test_df[["Molecule Name", "SMILES", "pEC50"]].to_csv("submission.csv", index=False)
    logger.info("Inference complete. Results saved to submission.csv")
