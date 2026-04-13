#!/usr/bin/env python3
import sys
import argparse
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from dimorphite_dl import protonate_smiles

def standardize_and_protonate(smiles: str, target_ph: float = 7.4) -> str:
    """Standardizes a SMILES string and applies pH-specific protonation."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("RDKit could not parse the SMILES string.")

    # 1. Keep the largest fragment (removes salts/solvents)
    mol = rdMolStandardize.FragmentParent(mol)

    # 2. Canonicalize the tautomer before protonation
    enumerator = rdMolStandardize.TautomerEnumerator()
    mol = enumerator.Canonicalize(mol)
    
    # Generate an intermediate clean SMILES for Dimorphite
    intermediate_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

    # 3. Apply physiological protonation using Dimorphite-DL
    # run_with_smiles returns a list of possible states in the pH range.
    # We set min and max to the same value to force a specific physiological window.
    try:
        protonated_smiles_list = protonate_smiles(
            intermediate_smiles, 
            ph_min=target_ph, 
            ph_max=target_ph,
        )
    except Exception as e:
        raise RuntimeError(f"Dimorphite-DL execution failed: {e}")

    if not protonated_smiles_list:
        # Fallback to the intermediate neutral/canonical state if the pKa model fails
        return intermediate_smiles
        
    # Dimorphite sorts the outputs; the first is the most probable state at this pH
    final_smiles = protonated_smiles_list[0]
    
    return final_smiles

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Standardize and protonate SMILES for Boltz.")
    parser.add_argument("--smiles", type=str, required=True, help="Input SMILES string")
    parser.add_argument("--ph", type=float, default=7.4, help="Target physiological pH")
    args = parser.parse_args()

    try:
        clean_smiles = standardize_and_protonate(args.smiles, args.ph)
        # Print only the resulting SMILES to stdout for bash capture
        print(clean_smiles)
    except Exception as e:
        # Print error to stderr and fail explicitly so the bash script catches it
        print(f"Error processing {args.smiles}: {e}", file=sys.stderr)
        sys.exit(1)
