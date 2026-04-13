#!/usr/bin/env bash

set -euo pipefail

# 1. Create a temporary staging directory and clear it if it already exists
mkdir -p submission_tmp
rm -f submission_tmp/*.pdb

# 2. Find files spepdbically in boltz_outputs/x*
find boltz_outputs/ -type f -path "*/x*/*_model_0.pdb" | while read -r filepath; do
    # Extract directory name (e.g., x00011-1)
    dirname=$(echo "$filepath" | grep -oE 'x[0-9]+-[0-9]+' | head -n1)
    
    # Define new filename
    newname="${dirname}.pdb"
    
    # Copy and rename into the staging directory
    cp "$filepath" "submission_tmp/$newname"
done

# 3. Navigate INSIDE the directory to avoid zipping the folder itself
cd submission_tmp

# Create the archive in the parent directory containing ONLY the .pdb files
# Use this for a ZIP file:
zip ../submission.zip *.pdb

# 4. Clean up (Optional but recommended)
cd ..
rm -rf submission_tmp

echo "Done! Check submission.zip"
