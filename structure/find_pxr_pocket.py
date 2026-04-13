#
#
# requires simplefold and p2rank
#
import os
import subprocess
import pandas as pd
import glob

# ================= Configuration =================
# Use absolute paths to ensure the subprocess commands locate the files correctly
FASTA_PATH = os.path.abspath("PXR_protein_sequence.fasta")
OUTPUT_DIR = os.path.abspath("./pxr_pocket_prediction")

# SimpleFold configuration
# Options: simplefold_100M, simplefold_360M, simplefold_700M, simplefold_1.1B, simplefold_1.6B, simplefold_3B
SIMPLEFOLD_MODEL = "simplefold_1.1B" 

# UPDATE THIS to the path where you extracted P2Rank
P2RANK_DIR = "p2rank_2.5.1" 
P2RANK_BIN = os.path.join(P2RANK_DIR, "prank")
# =================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

def main():
    # 1. Generate 3D Structure using SimpleFold
    print(f"Folding PXR sequence with Apple's SimpleFold ({SIMPLEFOLD_MODEL})...")
    
    # SimpleFold CLI command. It utilizes flow-matching on a transformer backbone without MSAs.
    simplefold_cmd = (
        f"simplefold --simplefold_model {SIMPLEFOLD_MODEL} "
        f"--num_steps 500 --tau 0.01 --nsample_per_protein 1 "
        f"--fasta_path {FASTA_PATH} --output_dir {OUTPUT_DIR} "
        f"--backend torch"
    )
    
    # We run the command inside the OUTPUT_DIR so any default outputs land there
    subprocess.run(simplefold_cmd, shell=True, check=True, cwd=OUTPUT_DIR)

    # SimpleFold may name the file based on the FASTA header or append indices (e.g., _0.pdb)
    # We dynamically locate the generated structure file
    generated_files = glob.glob(os.path.join(OUTPUT_DIR, "**/*.pdb"), recursive=True) + \
                      glob.glob(os.path.join(OUTPUT_DIR, "**/*.cif"), recursive=True)
    
    if not generated_files:
        raise FileNotFoundError(f"Could not find a generated PDB/CIF from SimpleFold in {OUTPUT_DIR}")
        
    structure_file = generated_files[0]
    print(f"Structure successfully generated: {structure_file}")

    # 2. Predict Binding Pockets using P2Rank
    print("\nRunning P2Rank geometry analysis...")
    p2rank_out = os.path.join(OUTPUT_DIR, 'p2rank_results')
    p2rank_cmd = f"{P2RANK_BIN} predict -f {structure_file} -o {p2rank_out}"
    subprocess.run(p2rank_cmd, shell=True, check=True)

    # 3. Parse the output and format for Umol
    # P2Rank names the CSV based on the input structure file name
    base_name = os.path.basename(structure_file)
    csv_path = os.path.join(p2rank_out, f"{base_name}_predictions.csv")
    
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"P2Rank output not found at {csv_path}")

    df = pd.read_csv(csv_path)
    
    # P2Rank ranks pockets by probability. The top pocket is index 0.
    top_pocket = df.iloc[0]
    raw_residues = top_pocket[' residue_ids'] # Format is usually like "A_10 A_11 A_12"

    # Umol requires a comma-separated list of 0-indexed residues
    clean_res_ids = []
    for res in raw_residues.split():
        # Extract the integer, ignoring the chain prefix (e.g. "A_142" -> 142)
        res_num = int(res.split('_')[1])
        # Convert to 0-index for Umol
        clean_res_ids.append(str(res_num - 1))

    formatted_pocket = ",".join(clean_res_ids)

    print("\n=== POCKET PREDICTION COMPLETE ===")
    print(f"Top Pocket Probability Score: {top_pocket[' probability']:.3f}")
    print(f"Top Pocket Center (x,y,z):    {top_pocket['   center_x']}, {top_pocket['   center_y']}, {top_pocket['   center_z']}")
    print(f"Umol Target Positions (0-indexed):\n{formatted_pocket}")

if __name__ == "__main__":
    main()
