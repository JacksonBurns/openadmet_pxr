import os
import subprocess
from pathlib import Path
from typing import List, Tuple

from Bio import pairwise2
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1


# ================= CONFIG =================
TEMPLATE_DIR = "raw_templates"      # input ENT/PDB files
OUTDIR = "processed_templates"
CHAIN_ID = "A"
FASTA_PATH = "PXR_protein_sequence.fasta"
P2RANK_PATH = "p2rank_2.5.1/prank"
TOP_N_POCKETS = 1
DISTANCE_THRESHOLD = 6.0


# ================= UTIL =================
def read_fasta_sequence(path: str) -> str:
    with open(path) as f:
        return "".join([line.strip() for line in f if not line.startswith(">")])

class ProteinOnlySelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        # Keep only standard protein residues
        return residue.id[0] == " "

    def accept_atom(self, atom):
        # Optionally filter altlocs (keep A or blank)
        altloc = atom.get_altloc()
        return altloc in (" ", "A")


def clean_and_extract_chain(pdb_path: str, out_path: str, chain_id: str):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)

    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, ProteinOnlySelect(chain_id))


def extract_sequence_from_pdb(pdb_path: str) -> Tuple[str, List[int]]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)

    seq = []
    res_ids = []

    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                try:
                    aa = seq1(res.resname)
                    seq.append(aa)
                    res_ids.append(res.id[1])
                except Exception:
                    continue
        break

    return "".join(seq), res_ids


def align_sequences(template_seq: str, target_seq: str):
    alignments = pairwise2.align.globalxx(template_seq, target_seq)
    return alignments[0]


def map_residues(template_seq, template_ids, target_seq):
    aln = align_sequences(template_seq, target_seq)
    t_aln, s_aln, _, _, _ = aln

    mapping = {}
    t_idx = 0
    s_idx = 0

    for t_char, s_char in zip(t_aln, s_aln):
        if t_char != "-":
            t_res_id = template_ids[t_idx]
            t_idx += 1
        else:
            t_res_id = None

        if s_char != "-":
            s_res_id = s_idx + 1  # FASTA indexing (1-based)
            s_idx += 1
        else:
            s_res_id = None

        if t_res_id and s_res_id:
            mapping[t_res_id] = s_res_id

    return mapping


def run_p2rank(pdb_path: str, out_dir: str):
    cmd = [
        P2RANK_PATH,
        "predict",
        "-f",
        pdb_path,
        "-o", out_dir
    ]
    subprocess.run(cmd, check=True)


def parse_p2rank_residues(prediction_dir: str) -> List[int]:
    # Reads residues from P2Rank output CSV
    csv_files = list(Path(prediction_dir).glob("*_residues.csv"))
    if not csv_files:
        raise RuntimeError("No P2Rank residue output found")

    residues = []
    with open(csv_files[0]) as f:
        next(f)
        for line in f:
            parts = line.strip().split(",")
            if int(parts[6]) == 1:
                resnum = int(parts[1])
                residues.append(resnum)

    return residues


# ================= MAIN =================
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    target_seq = read_fasta_sequence(FASTA_PATH)

    all_constraints = set()

    for pdb_file in Path(TEMPLATE_DIR).glob("*"):
        name = pdb_file.stem
        print(f"Processing {name}")

        clean_pdb = os.path.join(OUTDIR, f"{name}_chain{CHAIN_ID}.pdb")

        # Step 1: clean + chain select
        clean_and_extract_chain(str(pdb_file), clean_pdb, CHAIN_ID)

        # Step 2: extract sequence
        template_seq, template_ids = extract_sequence_from_pdb(clean_pdb)

        # Step 3: align to FASTA
        mapping = map_residues(template_seq, template_ids, target_seq)

        # Step 4: run P2Rank
        p2rank_out = os.path.join(OUTDIR, f"{name}_p2rank")
        os.makedirs(p2rank_out, exist_ok=True)
        run_p2rank(clean_pdb, p2rank_out)

        # Step 5: get pocket residues
        pocket_res = parse_p2rank_residues(p2rank_out)

        # Step 6: map to FASTA indices
        mapped = [mapping[r] for r in pocket_res if r in mapping]

        print(f"{name}: {len(mapped)} mapped pocket residues")

        # all_constraints.update(mapped)  # <-- keep all
        all_constraints.intersection_update(mapped)  # <-- keep only overlapping

    # Step 7: build YAML constraint
    sorted_res = sorted(all_constraints)

    yaml_contacts = "[[" + "], [".join([f"'A', {r}" for r in sorted_res]) + "]]"

    print("\n=== FINAL CONSTRAINT BLOCK ===\n")
    print(f"""constraints:
  - pocket:
      binder: B
      contacts: {yaml_contacts}
      max_distance: {DISTANCE_THRESHOLD}
      force: false
""")


if __name__ == "__main__":
    main()
