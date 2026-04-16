This is the model currently implemented in this repo:

Ensemble of the below models:

A. [Chemprop](https://doi.org/10.1021/acs.jcim.5c02332) (`chemprop.sh`)

Requires: `chemprop>=2.2.2`

1. pre-train chemprop on the primary screen with concentration as an input descriptor target is log2_fc_estimate
2. fine tune (2) on compounds put through primary assay

B. [CheMeleon](https://doi.org/10.48550/arXiv.2506.15792) (`chemeleon.sh`)

Requires: `chemprop>=2.2.2`

Direct training on the pec50 task

C. Physicochemical Random Forest (`random_forest.py`)

Requires: `molpipeline>=0.13.0`

From [`molpipeline`](https://doi.org/10.1021/acs.jcim.4c00863)'s `predefined_baselines`, on pEC50 only since multitask doesn't help in this architecture.


I use data augmentation inspired by the [RIGR paper](https://doi.org/10.1021/acs.jcim.5c00495) to teach the models to be resonance-invariant - see `preprocess_smiles.py`.
You will need a working installation of [RMG-Py](https://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/index.html) to run this script.


Later considerations:

 - can add uncertainty to the loss function with the known uncertainty values, possibly
