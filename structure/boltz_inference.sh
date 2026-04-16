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

mkdir -p "$OUTDIR_BASE"


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
      templates:
        - pdb: "$(realpath ./processed_templates/pdb9fzj_chainA.pdb)"
        - pdb: "$(realpath ./processed_templates/pdb8r00_chainA.pdb)"
        - pdb: "$(realpath ./processed_templates/pdb9fzj_chainA.pdb)"
        - pdb: "$(realpath ./processed_templates/pdb7axe_chainA.pdb)"
        - pdb: "$(realpath ./processed_templates/pdb4ny9_chainA.pdb)"
        - pdb: "$(realpath ./processed_templates/pdb4xhd_chainA.pdb)"
  - ligand:
      id: B
      smiles: '$CLEAN_SMILES'
# from prepare_templates_and_constraints.py
constraints:
  - pocket:
      binder: B
      contacts: [['A', 106], ['A', 144], ['A', 266]]
      max_distance: 6.0
      force: true
EOF

    # C. Run Boltz Prediction
    # We can add --use_potentials to enable the physics-based refinement for the pocket constraint.
    boltz predict "$INPUT_YAML" \
        --out_dir "$COMPLEX_OUTDIR" \
        --output_format pdb \
        --method "x-ray diffraction" \
        --diffusion_samples 10 \
        --max_parallel_samples 10 \
        --use_potentials

    end_time=$(date +%s)
    duration=$((end_time - start_time))
    echo "Finished $ligand_id in $duration seconds."
done
