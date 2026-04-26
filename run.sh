# run.sh
#
# Usage: bash -i run.sh
#
# Runs the activity training and inference, followed by zero-shot inference for the structure challenge.

set -euo pipefail

cd activity
conda activate autogluon
python fully_automated.py
cd ../structure
conda activate boltz
./boltz_inference.sh
./prepare_submission.sh
