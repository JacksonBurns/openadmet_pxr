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
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover

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

if __name__ == "__main__":
    from tqdm import tqdm


    train_df = pd.read_csv("train.csv")
    train_df["SMILES"] = train_df["SMILES"].astype(object)  # allows writing lists of SMILES here
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
    train_df = train_df.explode("SMILES")
    train_df.to_csv("train_augmented.csv", index=False)
