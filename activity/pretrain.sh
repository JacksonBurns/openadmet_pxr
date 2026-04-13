chemprop train \
    --output-dir pretrain_output \
    --logfile pretrain_output/log.txt \
    --data-path pretrain.csv \
    --split-sizes 0.8 0.1 0.1 \
    --descriptors-columns concentration_M \
    --from-foundation CheMeleon \
    --pytorch-seed 42 \
    --smiles-columns SMILES \
    --target-columns log2_fc_estimate log2_fc_stderr \
    --task-type regression \
    --patience 5 \
    --init-lr 0.0001 \
    --warmup-epochs 5 \
    --loss mse \
    --metrics rmse r2 mse mae \
    --show-individual-scores \
    --ffn-num-layers 2 \
    --ffn-hidden-dim 512 \
    --batch-size 32 \
    --epochs 50

# if stepping away from CheMeleon, change the FNN settings...
