import argparse

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input CSV file")
    parser.add_argument("--output", help="Output CSV file")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df[["Molecule Name", "SMILES", "pEC50"]].to_csv(args.output, index=False)
