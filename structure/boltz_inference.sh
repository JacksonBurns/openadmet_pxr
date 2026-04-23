#!/bin/bash
#
# Runs Boltz with Precomputed MSA, Dynamic Templates, and Constraints
#
# Usage: ./boltz_inference.sh 

set -euo pipefail

# ================= Configuration =================
FASTA_PATH="PXR_protein_sequence.fasta"
SMILES_CSV="pxr-challenge_structure_TEST_BLINDED.csv"
OUTDIR_BASE="boltz_outputs"
PROCESSED_TEMPLATES="processed_templates"
RAW_TEMPLATES="raw_templates"
MSA_PATH="$(realpath mmseqs2_pxr.a3m)"  

mkdir -p "$OUTDIR_BASE"
mkdir -p "$PROCESSED_TEMPLATES"

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

    # B. Generate the dynamic YAML blocks (Templates and Constraints)
    echo "Calculating optimal structural templates..."
    INPUT_YAML="${COMPLEX_OUTDIR}/${ligand_id}_input.yaml"
    
    # 1. Write the static headers
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

    # 2. Append the dynamic blocks directly to the YAML
    if ! python3 dynamic_template_selector.py \
        --query_smiles "$CLEAN_SMILES" \
        --fasta "$FASTA_PATH" \
        --raw_dir "$RAW_TEMPLATES" \
        --out_dir "$PROCESSED_TEMPLATES" \
        --chain "A" \
        --p2rank "p2rank_2.5.1/prank" >> "$INPUT_YAML"; then
        echo "WARNING: Dynamic template generation failed. Skipping..."
        continue
    fi

    # C. Run Boltz Prediction
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
