alias cmln_env='conda run --no-capture-output -n chemprop_live'
alias rmgpy_env='conda run --no-capture-output -n rmg_env'

cmln_env python get_data.py
rmgpy_env python preprocess_smiles.py --do-train --do-test

# fit on the initial training data
for fold in {0..4}; do
    subdir="fold_${fold}"
    cmln_env chemprop train \
        --from-foundation CheMeleon \
        --ffn-num-layers 2 \
        --ffn-hidden-dim 512 \
        --batch-size 64 \
        --epochs 50 \
        --patience 5 \
        --pytorch-seed 42 \
        --data-path \
            splits/${subdir}/train.csv \
            splits/${subdir}/val.csv \
            splits/${subdir}/val.csv \
        --smiles-columns SMILES \
        --target-columns pEC50 \
        --loss mse \
        --metrics rmse r2 mse mae \
        --task-type regression \
        --weight-column pEC50_weight \
        --show-individual-scores \
        --output-dir train_output/fold_${fold} \
        --init-lr 0.00001 \
        --max-lr 0.0001 \
        --final-lr 0.00001 \
        --warmup-epochs 2 \
        --num-workers 8 \
        --dropout 0.50 \
        --batch-norm \
        --molecule-featurizers morgan_count rdkit_2d
done

# run prediction on the training data
cmln_env chemprop predict \
    --model-paths train_output \
    --test-path train_augmented.csv \
    --preds-path train_output/training_predictions.csv \
    --molecule-featurizers morgan_count rdkit_2d \
    --smiles-columns SMILES

# generate denoised training data
cmln_env python denoise.py

# tell the user to run preprocess again, waiting for them to confirm it is done
rmgpy_env python preprocess_smiles.py --do-train-denoised

# fit on the denoised training data
for fold in {0..4}; do
    subdir="fold_${fold}"
    cmln_env chemprop train \
        --from-foundation CheMeleon \
        --ffn-num-layers 2 \
        --ffn-hidden-dim 512 \
        --batch-size 64 \
        --epochs 50 \
        --patience 5 \
        --pytorch-seed 42 \
        --data-path \
            splits_denoised/${subdir}/train.csv \
            splits_denoised/${subdir}/val.csv \
            splits_denoised/${subdir}/val.csv \
        --smiles-columns SMILES \
        --target-columns pEC50 \
        --loss mse \
        --metrics rmse r2 mse mae \
        --task-type regression \
        --weight-column pEC50_weight \
        --show-individual-scores \
        --output-dir train_output_denoised/fold_${fold} \
        --init-lr 0.00001 \
        --max-lr 0.0001 \
        --final-lr 0.00001 \
        --warmup-epochs 2 \
        --num-workers 8 \
        --dropout 0.50 \
        --batch-norm \
        --molecule-featurizers morgan_count rdkit_2d
done

# final inference
cmln_env chemprop predict \
    --model-paths train_output_denoised \
    --test-path test_augmented.csv \
    --preds-path train_output_denoised/predictions.csv \
    --molecule-featurizers morgan_count rdkit_2d \
    --smiles-columns SMILES

cmln_env python submit.py
