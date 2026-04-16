Disclaimer - I am very new to co-folding and am basically just following tutorials and LLMing my way to the finish line.
There are entire sections of this which may make no sense, or worse yet, be actively the totally wrong thing to do.
Please tell me about them!

To run this code, install `boltz`:

```bash
pip install boltz[cuda]
pip uninstall nvidia-cuda-nvrtc
conda install -c nvidia cuda-nvrtc
```

I found this structure that is similar to the input FASTA:

https://www.uniprot.org/uniprotkb/O75469/entry

using this FASTA similarity search tool:

https://www.ebi.ac.uk/jdispatcher/sss

With `biopython` installed run `prepare_templates_and_constraints.py` (to reproduce from scratch, I may just include these in the repo to avoid this).
