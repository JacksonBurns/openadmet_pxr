#!/bin/bash
#
# Runs Boltz with Precomputed MSA and Pocket Constraints
#
# Usage: ./boltz_inference.sh 

set -euo pipefail

# ================= Configuration =================
FASTA_PATH="PXR_protein_sequence.fasta"
SMILES_CSV="pxr-challenge_structure_TEST_BLINDED.csv"
OUTDIR_BASE="boltz_outputs"
MSA_PATH="$(realpath mmseqs2_pxr.a3m)"  # Use absolute path for the MSA

# Pocket residues from your previous analysis
# POCKET_RESIDUES="101,102,104,105,109,139,143,146,157,164,181,185,265,268,269,272,278,283,287,63,64,67,68,69,97,98"

mkdir -p "$OUTDIR_BASE"

# Convert comma-separated residues into YAML list format: [['A', 101], ['A', 102], ...]
# This assumes Chain A is the protein
# YAML_CONTACTS=$(echo "$POCKET_RESIDUES" | sed "s/,/], ['A', /g; s/^/[['A', /; s/$/]]/")

# Extract the protein sequence
PROTEIN_SEQ=$(grep -v "^>" "$FASTA_PATH" | tr -d '\n' | tr -d '\r')

tail -n +2 "$SMILES_CSV" | while IFS=',' read -r ligand_id raw_smiles; do
    start_time=$(date +%s)

    ligand_id=$(echo "$ligand_id" | tr -d '\r' | xargs)
    raw_smiles=$(echo "$raw_smiles" | tr -d '\r' | xargs)

    COMPLEX_OUTDIR="${OUTDIR_BASE}/${ligand_id}"
    mkdir -p "$COMPLEX_OUTDIR"

    echo -e "\n--- Processing Ligand: $ligand_id ---"

    # A. Ligand Preparation
    if ! CLEAN_SMILES=$(python3 prepare_ligand.py --smiles "$raw_smiles" --ph 7.4); then
        echo "WARNING: Ligand prep failed for $ligand_id. Skipping..."
        continue
    fi

    # B. Generate Boltz YAML input with Constraints and MSA
    INPUT_YAML="${COMPLEX_OUTDIR}/${ligand_id}_input.yaml"
    
    cat <<EOF > "$INPUT_YAML"
version: 1
sequences:
  - protein:
      id: A
      sequence: $PROTEIN_SEQ
      msa: "$MSA_PATH"
  - ligand:
      id: B
      smiles: '$CLEAN_SMILES'
EOF
# constraints:
#   - pocket:
#       binder: B
#       contacts: $YAML_CONTACTS
#       max_distance: 6.0
#       force: false  # Set to true if you want to enforce the constraint more strongly

    # C. Run Boltz Prediction
    # We remove --use_msa_server because the YAML now points to a local file.
    # We add --use_potentials to enable the physics-based refinement for the pocket constraint.
    boltz predict "$INPUT_YAML" \
        --out_dir "$COMPLEX_OUTDIR" \
        --output_format pdb \
        --method "x-ray diffraction"
        # --use_potentials \

    end_time=$(date +%s)
    duration=$((end_time - start_time))
    echo "Finished $ligand_id in $duration seconds."
done
