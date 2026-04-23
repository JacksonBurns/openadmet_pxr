# run.sh
#
# Usage: bash -i run.sh
#
# Runs the activity training and inference, followed by zero-shot inference for the structure challenge.

set -euo pipefail

cd activity
conda activate chemprop_live
./chemprop.sh
./chemeleon.sh
python random_forest.py
python prepare_submission.py --output submission.csv
# cd ../structure
# conda activate boltz
# ./boltz_inference.sh
# ./prepare_submission.sh
