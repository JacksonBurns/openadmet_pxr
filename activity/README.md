This repository implements a [CheMeleon](https://doi.org/10.48550/arXiv.2506.15792) model with extra descriptors inspired by [`molpipeline`](https://doi.org/10.1021/acs.jcim.4c00863)'s `predefined_baselines`.

I use data augmentation inspired by the [RIGR paper](https://doi.org/10.1021/acs.jcim.5c00495) to teach the models to be resonance-invariant - see `preprocess_smiles.py`.
You will need a working installation of [RMG-Py](https://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/index.html) to run this script.

I also use denoising inspired by [this paper](doi.org/10.1021/acs.jcim.4c00639).

The main driver script is in `fit.sh` - you will need two `conda` environments, one with `rmg` and one with `Chemprop`.
Getting RMG installed is more involved, so follow the install tutorial linked above.

Once that's running, you just need to run `fit.sh`.
It takes care of training, denoising, retraining, and inference.
If you're curious about how the model is configured in greater detail than I mention above, I suggest poking around in this and the other scripts.
