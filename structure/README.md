Disclaimer - I am very new to co-folding and am basically just following tutorials and LLMing my way to the finish line.
There are entire sections of this which may make no sense, or worse yet, be actively the totally wrong thing to do.
Please tell me about them!

To run this code, install `boltz`:

```bash
pip install boltz[cuda]
pip uninstall nvidia-cuda-nvrtc
conda install -c nvidia cuda-nvrtc
```

Currently unused:

I initially wrote up some code to try and get the protein pocket as a constraint for Boltz, I might go back to this in the future:

`pip install git+https://github.com/apple/ml-simplefold.git`

`python find_pxr_pocket.py`
