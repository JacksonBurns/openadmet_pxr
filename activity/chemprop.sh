#!/usr/bin/env bash

set -euo pipefail

mkdir -p train_output
mkdir train_output/chemprop

for fold in {0..4}; do
    subdir="fold_${fold}"
    mkdir -p train_output/chemprop/${subdir}
    for i in {0..3}; do
        chemprop train \
            --output-dir train_output/chemprop/${subdir}/split_${i} \
            --logfile train_output/chemprop/${subdir}/split_${i}/log.txt \
            --data-path \
                splits/${subdir}/split_${i}/train.csv \
                splits/${subdir}/split_${i}/val.csv \
                splits/${subdir}/split_${i}/val.csv \
            --pytorch-seed 42 \
            --smiles-columns SMILES \
            --target-columns \
                pEC50 \
                "Emax_estimate (log2FC vs. baseline)" \
                pEC50_counter \
                "Emax_estimate (log2FC vs. baseline)_counter" \
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
        --model-paths train_output/chemprop/${subdir} \
        --test-path splits/${subdir}/test.csv \
        --preds-path train_output/chemprop_${subdir}_cv_predictions.csv \
        --smiles-columns SMILES
    
    chemprop predict \
        --model-paths train_output/chemprop/${subdir} \
        --test-path test_augmented.csv \
        --preds-path train_output/chemprop_${subdir}_holdout_predictions.csv \
        --smiles-columns SMILES
done
