#!/usr/bin/env bash

set -euo pipefail

mkdir -p train_output
mkdir train_output/chemprop_pretrain
mkdir train_output/chemprop

chemprop train \
    --output-dir train_output/chemprop_pretrain \
    --logfile train_output/chemprop_pretrain/log.txt \
    --data-path pretrain.csv \
    --split-sizes 0.8 0.1 0.1 \
    --descriptors-columns concentration_M \
    --pytorch-seed 42 \
    --smiles-columns SMILES \
    --target-columns log2_fc_estimate \
    --task-type regression \
    --patience 5 \
    --init-lr 0.0001 \
    --warmup-epochs 5 \
    --loss mse \
    --metrics rmse r2 mse mae \
    --show-individual-scores \
    --message-hidden-dim 512 \
    --ffn-num-layers 2 \
    --ffn-hidden-dim 512 \
    --batch-size 32 \
    --epochs 50

for i in {0..2}; do
    chemprop train \
        --output-dir train_output/chemprop/split_${i} \
        --logfile train_output/chemprop/split_${i}/log.txt \
        --data-path \
            splits/split_${i}_train.csv \
            splits/split_${i}_val.csv \
            splits/split_${i}_val.csv \
        --from-foundation train_output/chemprop_pretrain/model_0/best.pt \
        --pytorch-seed 42 \
        --smiles-columns SMILES \
        --target-columns pEC50 \
        --task-weights 10 1 1 1 \
        --task-type regression \
        --patience 5 \
        --init-lr 0.00001 \
        --warmup-epochs 5 \
        --loss mse \
        --metrics rmse r2 mse mae \
        --show-individual-scores \
        --ffn-num-layers 2 \
        --ffn-hidden-dim 512 \
        --batch-size 32 \
        --epochs 50
done

chemprop predict \
    --model-paths train_output/chemprop \
    --test-path test.csv \
    --preds-path train_output/chemprop_predictions.csv \
    --smiles-columns SMILES
