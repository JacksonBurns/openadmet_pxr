#!/usr/bin/env bash

set -euo pipefail

mkdir -p train_output
mkdir train_output/chemprop


for i in {0..3}; do
    chemprop train \
        --output-dir train_output/chemprop/split_${i} \
        --logfile train_output/chemprop/split_${i}/log.txt \
        --data-path \
            splits/split_${i}_train.csv \
            splits/split_${i}_val.csv \
            splits/split_${i}_val.csv \
        --pytorch-seed 42 \
        --smiles-columns SMILES \
        --target-columns pEC50 \
        --weight-column pEC50_weight \
        --task-type regression \
        --patience 5 \
        --init-lr 0.00001 \
        --warmup-epochs 5 \
        --loss mse \
        --metrics rmse r2 mse mae \
        --show-individual-scores \
        --ffn-num-layers 2 \
        --ffn-hidden-dim 512 \
        --message-hidden-dim 512 \
        --batch-size 32 \
        --epochs 50
done

chemprop predict \
    --model-paths train_output/chemprop \
    --test-path test.csv \
    --preds-path train_output/chemprop_predictions.csv \
    --smiles-columns SMILES
