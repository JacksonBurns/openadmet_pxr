# preprocess_smiles.py
#
# applies preprocessing to the challenge SMILES using RDKit to get the 'parent graph'
# and RMG to perform resonance structure augmentation
#
# the former is taken as a whole from:
# https://github.com/JacksonBurns/openadmet_expansionrx/blob/34ff4bb5230e0fbc42a5638d194142ab189964aa/get_data.py
# the latter is inspired by, but not directly adapted from:
# https://github.com/akshatzalte/rigr/blob/e315d8ab412473afcd3a07f24c8839e9b4b93a99/notebooks/resonance_generation_and_augmentation.ipynb

import pandas as pd
from astartes.molecules import train_test_split_molecules
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from sklearn.model_selection import KFold

def clean_smiles(
    smiles: str, remove_hs: bool = True, strip_stereochem: bool = False, strip_salts: bool = True
) -> str:
    """Applies preprocessing to SMILES strings, seeking the 'parent' SMILES

    Note that this is different from simply _neutralizing_ the input SMILES - we attempt to get the parent molecule, analogous to a molecular skeleton.
    This is adapted in part from https://rdkit.org/docs/Cookbook.html#neutralizing-molecules

    Args:
        smiles (str): input SMILES
        remove_hs (bool, optional): Removes hydrogens. Defaults to True.
        strip_stereochem (bool, optional): Remove R/S and cis/trans stereochemistry. Defaults to False.
        strip_salts (bool, optional): Remove salt ions. Defaults to True.

    Returns:
        str: cleaned SMILES
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Could not parse SMILES {smiles}"
        if remove_hs:
            mol = Chem.RemoveHs(mol)
        if strip_stereochem:
            Chem.RemoveStereochemistry(mol)
        if strip_salts:
            remover = SaltRemover()  # use default saltremover
            mol = remover.StripMol(mol)  # strip salts

        pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()
        out_smi = Chem.MolToSmiles(mol, kekuleSmiles=True)  # this also canonicalizes the input
        assert len(out_smi) > 0, f"Could not convert molecule to SMILES {smiles}"
        return out_smi
    except Exception as e:
        print(f"Failed to clean SMILES {smiles} due to {e}")
        return None

from rmgpy.molecule import Molecule

def resonate(smiles):
    mol = Molecule().from_smiles(smiles)
    structs = mol.generate_resonance_structures()
    return [s.to_smiles() for s in structs]

def explode_and_reweight(df, smiles_col="SMILES", weight_col="pEC50_weight"):
    """explodes on resonance structures and divides weights by the number of resonance structures for each molecule, so that the total weight for each molecule is unchanged by resonance augmentation"""
    exploded_df = df.explode(smiles_col).reset_index(drop=True)
    exploded_df[weight_col] = exploded_df[weight_col] / exploded_df.groupby("Molecule Name")[smiles_col].transform("count")
    return exploded_df

if __name__ == "__main__":
    from pathlib import Path

    from tqdm import tqdm

    test_df = pd.read_csv("test.csv")
    test_df["SMILES"] = test_df["SMILES"].astype(object)
    for i in tqdm(range(test_df.shape[0]), desc="Preprocessing SMILES"):
        og_smiles = test_df.iloc[i]['SMILES']
        try:
            clean_smi = clean_smiles(og_smiles)
        except Exception as e:
            print(f"Skipping {og_smiles}, failed initial clean")
            print(e)
            continue
        try:
            resonance_smiles = resonate(clean_smi)
            test_df.at[i, "SMILES"] = resonance_smiles
        except Exception as e:
            print(f"Skipping resonance generation for smiles {og_smiles}")
            print(e)
    test_df.explode("SMILES").to_csv("test_augmented.csv", index=False)

    train_df = pd.read_csv("train.csv")
    train_df["SMILES"] = train_df["SMILES"].astype(object)
    for i in tqdm(range(train_df.shape[0]), desc="Preprocessing SMILES"):
        og_smiles = train_df.iloc[i]['SMILES']
        try:
            clean_smi = clean_smiles(og_smiles)
        except Exception as e:
            print(f"Skipping {og_smiles}, failed initial clean")
            print(e)
            continue
        try:
            resonance_smiles = resonate(clean_smi)
            train_df.at[i, "SMILES"] = resonance_smiles
        except Exception as e:
            print(f"Skipping resonance generation for smiles {og_smiles}")
            print(e)
    train_df.explode("SMILES").to_csv("train_augmented.csv", index=False)

    outdir = Path("splits")
    outdir.mkdir(exist_ok=True)
    for fold_number, fold in enumerate(KFold(n_splits=5, shuffle=True, random_state=42).split(train_df)):
        subdir = outdir / f"fold_{fold_number}"
        subdir.mkdir(exist_ok=True)
        explode_and_reweight(train_df.iloc[fold[1]].reset_index(drop=True)).to_csv(subdir / "test.csv", index=False)
        subdf = train_df.iloc[fold[0]].reset_index(drop=True)
        explode_and_reweight(subdf).to_csv(subdir / "train_val.csv", index=False)  # for models that don't need a separate validation set, they can just use this whole fold for training
        for i in range(4):
            subsubdir = subdir / f"split_{i}"
            subsubdir.mkdir(exist_ok=True)
            *_, train_idxs, val_idxs = train_test_split_molecules([s if isinstance(s, str) else s[0] for s in subdf["SMILES"].to_list()], train_size=0.80, test_size=0.20, sampler="kmeans", random_state=i + 42, return_indices=True)
            i_train_df = subdf.iloc[train_idxs].reset_index(drop=True)
            i_val_df = subdf.iloc[val_idxs].reset_index(drop=True)
            i_val_df = explode_and_reweight(i_val_df)
            i_val_df.to_csv(subsubdir / "val.csv", index=False)
            i_train_df = explode_and_reweight(i_train_df)
            i_train_df.to_csv(subsubdir / "train.csv", index=False)
