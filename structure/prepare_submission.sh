#!/usr/bin/env bash

set -euo pipefail

mkdir -p submission_tmp
rm -f submission_tmp/*.pdb

echo "Evaluating ensembles for highest ligand_iptm..."

for ligand_dir in boltz_outputs/x*; do
    [ -d "$ligand_dir" ] || continue
    
    ligand_id=$(basename "$ligand_dir")
    
    best_iptm="-1.0"
    best_pdb=""

    # Use find to locate all confidence JSONs anywhere inside this ligand's folder
    while IFS= read -r json_file; do
        
        # Extract the score: grep the line, split by colon, remove spaces and commas
        current_iptm=$(grep '"ligand_iptm"' "$json_file" | awk -F':' '{print $2}' | tr -d ' ,')
        
        # Skip if grep failed to find the line
        [ -z "$current_iptm" ] && continue
        
        # Compare floating point numbers using awk
        is_better=$(awk -v curr="$current_iptm" -v best="$best_iptm" 'BEGIN {if (curr > best) print 1; else print 0}')
        
        if [ "$is_better" -eq 1 ]; then
            best_iptm="$current_iptm"
            
            # The PDB is in the exact same directory as the JSON
            json_dir=$(dirname "$json_file")
            base_name=$(basename "$json_file")
            
            # Swap prefix and extension
            pdb_name=$(echo "$base_name" | sed 's/confidence_//; s/\.json$/.pdb/')
            best_pdb="$json_dir/$pdb_name"
        fi
        
    done < <(find "$ligand_dir" -type f -name "confidence_*.json")

    # Copy and rename if we found a winner
    if [ -n "$best_pdb" ] && [ -f "$best_pdb" ]; then
        cp "$best_pdb" "submission_tmp/${ligand_id}.pdb"
        echo "Selected top pose for $ligand_id -> $(basename "$best_pdb") (score: $best_iptm)"
    else
        # This will correctly catch cases like x00035-1 where the prediction failed/crashed
        echo "WARNING: Could not process $ligand_id (no output files found)"
    fi
done

echo -e "\nPackaging submission..."

cd submission_tmp
zip -q ../submission.zip *.pdb
cd ..
rm -rf submission_tmp

echo "Done! Check submission.zip"
