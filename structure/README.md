Disclaimer - I am very new to co-folding and am basically just following tutorials and LLMing my way to the finish line.
There are entire sections of this which may make no sense, or worse yet, be actively the totally wrong thing to do.
Please tell me about them!

My general approach is to use Boltz for structure generation, using templates of experimentally measured PXR complexes and `p2rank` to suggest which residues should be used as binding constraints for the ligands.
The appropriate templates are identified based on their ligands' similarity to the ligands in the test data.

To run the structure generation pipeline (`boltz_inference.sh`), install `boltz` (I have an extra workaround included to get it working on my machine, may not be needed for you):

```bash
pip install boltz[cuda]
pip uninstall nvidia-cuda-nvrtc
conda install -c nvidia cuda-nvrtc
```

On top of that, you will need to install `dimorphite-dl` and BioPython, and download the `p2rank` binary into this directory.

All of the other required plaintext files are already in the repository.
To identify useful starting templates for Boltz, I found this structure that is similar to the input FASTA: https://www.uniprot.org/uniprotkb/O75469/entry using this FASTA similarity search tool: https://www.ebi.ac.uk/jdispatcher/sss

The MSA was generated using the MMseqs2 server.
