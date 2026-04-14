for i in {0..2}; do
    chemprop train \
        --output-dir train_output/chemeleon/split_${i} \
        --logfile train_output/chemeleon/split_${i}/log.txt \
        --data-path \
            splits/split_${i}_train.csv \
            splits/split_${i}_val.csv \
            splits/split_${i}_val.csv \
        --from-foundation CheMeleon \
        --pytorch-seed 42 \
        --smiles-columns SMILES \
        --target-columns "pEC50" "Emax_estimate (log2FC vs. baseline)" "pEC50_counter" "Emax_estimate (log2FC vs. baseline)_counter" \
        --task-weights 10 1 1 1 \
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
    --model-paths train_output/chemeleon \
    --test-path test.csv \
    --preds-path train_output/chemeleon_predictions.csv \
    --smiles-columns SMILES
