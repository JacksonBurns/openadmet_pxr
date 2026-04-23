#!/usr/bin/env bash

set -euo pipefail

mkdir -p train_output
mkdir train_output/chemeleon

for fold in {0..4}; do
    subdir="fold_${fold}"
    mkdir -p train_output/chemeleon/${subdir}
    for i in {0..3}; do
        chemprop train \
            --output-dir train_output/chemeleon/${subdir}/split_${i} \
            --logfile train_output/chemeleon/${subdir}/split_${i}/log.txt \
            --data-path \
                splits/${subdir}/split_${i}/train.csv \
                splits/${subdir}/split_${i}/val.csv \
                splits/${subdir}/split_${i}/val.csv \
            --from-foundation CheMeleon \
            --pytorch-seed 42 \
            --smiles-columns SMILES \
            --target-columns \
                pEC50 \
                "Emax_estimate (log2FC vs. baseline)" \
                pEC50_counter \
                "Emax_estimate (log2FC vs. baseline)_counter" \
            --weight-column pEC50_weight \
            --task-type regression \
            --patience 3 \
            --max-lr 0.0001 \
            --init-lr 0.000001 \
            --warmup-epochs 2 \
            --loss mse \
            --metrics rmse r2 mse mae \
            --show-individual-scores \
            --ffn-num-layers 1 \
            --ffn-hidden-dim 512 \
            --batch-size 32 \
            --epochs 20
    done
    chemprop predict \
        --model-paths train_output/chemeleon/${subdir} \
        --test-path splits/${subdir}/test.csv \
        --preds-path train_output/chemeleon_${subdir}_cv_predictions.csv \
        --smiles-columns SMILES
    
    chemprop predict \
        --model-paths train_output/chemeleon/${subdir} \
        --test-path test_augmented.csv \
        --preds-path train_output/chemeleon_${subdir}_holdout_predictions.csv \
        --smiles-columns SMILES
done
