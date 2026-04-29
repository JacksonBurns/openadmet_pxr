#!/usr/bin/env python3
import os
import sys
import json
import argparse
import subprocess
import warnings
from pathlib import Path
from typing import List, Tuple
from collections import Counter

from Bio import BiopythonWarning, BiopythonDeprecationWarning
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

from Bio import pairwise2
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

def log(msg: str):
    print(msg, file=sys.stderr)

class ProteinOnlySelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id
    def accept_residue(self, residue):
        return residue.id[0] == " "
    def accept_atom(self, atom):
        return atom.get_altloc() in (" ", "A")

class LigandOnlySelect(Select):
    def accept_residue(self, residue):
        return residue.id[0].startswith("H_") and residue.resname not in ["HOH", "WAT", "DOD", "NA", "CL", "MG", "ZN"]

def read_fasta_sequence(path: str) -> str:
    with open(path) as f:
        return "".join([line.strip() for line in f if not line.startswith(">")])

def clean_and_extract_chain(pdb_path: str, out_path: str, chain_id: str):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, ProteinOnlySelect(chain_id))

def extract_sequence_from_pdb(pdb_path: str) -> Tuple[str, List[int]]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)
    seq, res_ids = [], []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] == " ":
                    try:
                        seq.append(seq1(res.resname))
                        res_ids.append(res.id[1])
                    except Exception:
                        continue
        break
    return "".join(seq), res_ids

def map_residues(template_seq, template_ids, target_seq):
    aln = pairwise2.align.globalxx(template_seq, target_seq)[0]
    t_aln, s_aln = aln[0], aln[1]
    mapping, t_idx, s_idx = {}, 0, 0
    for t_char, s_char in zip(t_aln, s_aln):
        t_res_id = template_ids[t_idx] if t_char != "-" else None
        if t_char != "-": t_idx += 1
        s_res_id = s_idx + 1 if s_char != "-" else None
        if s_char != "-": s_idx += 1
        if t_res_id and s_res_id:
            mapping[t_res_id] = s_res_id
    return mapping

def run_p2rank(p2rank_path: str, pdb_path: str, out_dir: str):
    subprocess.run([p2rank_path, "predict", "-f", pdb_path, "-o", out_dir], check=True, stdout=subprocess.DEVNULL)

def parse_p2rank_residues(prediction_dir: str) -> List[int]:
    csv_files = list(Path(prediction_dir).glob("*_residues.csv"))
    if not csv_files: return []
    residues = []
    with open(csv_files[0]) as f:
        next(f)
        for line in f:
            parts = line.strip().split(",")
            if int(parts[6]) == 1:
                residues.append(int(parts[1]))
    return residues

def get_template_smiles(pdb_path: str, outdir: str) -> str:
    """Caches the template ligand as a SMILES string to avoid repeated PDB parsing."""
    cache_path = os.path.join(outdir, "template_smiles.json")
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            cache = json.load(f)
            
    name = Path(pdb_path).stem
    if name in cache:
        return cache[name]
        
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)
    io = PDBIO()
    io.set_structure(structure)
    
    tmp_ligand = os.path.join(outdir, f"{name}_tmp_ligand.pdb")
    io.save(tmp_ligand, LigandOnlySelect())
    
    mol = Chem.MolFromPDBFile(tmp_ligand, sanitize=False)
    if os.path.exists(tmp_ligand):
        os.remove(tmp_ligand)
        
    if mol:
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.GetSSSR(mol) # Ensures non-symmetrized ring perception is locked
            smi = Chem.MolToSmiles(mol)
            cache[name] = smi
            with open(cache_path, "w") as f:
                json.dump(cache, f)
            return smi
        except Exception:
            pass
            
    cache[name] = None
    with open(cache_path, "w") as f:
        json.dump(cache, f)
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query_smiles", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--raw_dir", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--p2rank", required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    
    query_mol = Chem.MolFromSmiles(args.query_smiles)
    if not query_mol:
        log("ERROR: Could not parse query SMILES.")
        sys.exit(1)
    
    Chem.GetSSSR(query_mol)
    query_heavy_atoms = query_mol.GetNumHeavyAtoms()
    
    scores = []
    for pdb_file in Path(args.raw_dir).glob("*"):
        if not pdb_file.is_file() or pdb_file.name.startswith("."): 
            continue
        
        tmpl_smi = get_template_smiles(str(pdb_file), args.out_dir)
        if tmpl_smi:
            tmpl_mol = Chem.MolFromSmiles(tmpl_smi)
            if tmpl_mol:
                Chem.GetSSSR(tmpl_mol)
                
                # Compute Maximum Common Substructure (MCS)
                # timeout prevents hangs on exceptionally complex graphs
                mcs_res = rdFMCS.FindMCS([query_mol, tmpl_mol], timeout=2)
                
                # Score is the fraction of the query's heavy atoms present in the template ligand
                sim = mcs_res.numAtoms / query_heavy_atoms
                scores.append((sim, pdb_file))
                continue
                
        scores.append((0.0, pdb_file))
            
    if not scores:
        log(f"ERROR: No structural files found in {args.raw_dir}")
        sys.exit(1)

    scores.sort(key=lambda x: x[0], reverse=True)
    # select up to 5 templates with highest MCS similarity, but only if they have at least 20% of the query's heavy atoms in common
    top_templates = [p for sim, p in scores if sim >= 0.2][:(5 if len(scores) >= 5 else len(scores))]
    log(f"Selected Top {len(top_templates)} Templates by MCS Overlap:")
    for sim, p in scores[:len(top_templates)]: log(f" - {p.stem} (Fraction of Ligand Matched: {sim:.3f})")

    target_seq = read_fasta_sequence(args.fasta)
    pocket_sets = []
    
    for pdb_file in top_templates:
        name = pdb_file.stem
        clean_pdb = os.path.join(args.out_dir, f"{name}_chain{args.chain}.pdb")
        p2rank_out = os.path.join(args.out_dir, f"{name}_p2rank")
        
        if not os.path.exists(clean_pdb):
            clean_and_extract_chain(str(pdb_file), clean_pdb, args.chain)
            
        template_seq, template_ids = extract_sequence_from_pdb(clean_pdb)
        mapping = map_residues(template_seq, template_ids, target_seq)
        
        csv_files = list(Path(p2rank_out).glob("*_residues.csv"))
        if not csv_files:
            os.makedirs(p2rank_out, exist_ok=True)
            run_p2rank(args.p2rank, clean_pdb, p2rank_out)
            csv_files = list(Path(p2rank_out).glob("*_residues.csv"))
            
        if csv_files:
            pocket_res = parse_p2rank_residues(p2rank_out)
            mapped = set([mapping[r] for r in pocket_res if r in mapping])
            if mapped:
                pocket_sets.append(mapped)

    if pocket_sets:
        all_pocket_res = [res for subset in pocket_sets for res in subset]
        res_counts = Counter(all_pocket_res)
        sorted_by_freq = sorted(res_counts.keys(), key=lambda r: (-res_counts[r], r))
        
        log("\nPocket Residue Frequencies:")
        for r in sorted_by_freq:
            log(f" - Residue {r} ({target_seq[r-1]}): {res_counts[r]} out of {len(pocket_sets)} templates")
    else:
        log("ERROR: No pockets found in selected templates.")
        sys.exit(1)

    polar_amino_acids = {'S', 'Q', 'H', 'R', 'T', 'Y', 'N', 'D', 'E', 'K'}
    anchor_residues = [r for r in sorted_by_freq if target_seq[r - 1] in polar_amino_acids]

    if not anchor_residues:
        log("WARNING: No polar residues found. Falling back to most frequent overall residues.")
        anchor_residues = sorted_by_freq
        
    final_targets = sorted(anchor_residues[:3])
    yaml_contacts = "[[" + "], [".join([f"'{args.chain}', {r}" for r in final_targets]) + "]]"

    print("templates:")
    for raw_path in top_templates:
        ext = raw_path.suffix.lower()
        file_type = "cif" if ext == ".cif" else "pdb"
        print(f"  - {file_type}: \"{os.path.abspath(raw_path)}\"")
        print(f"    chain_id: {args.chain}")
        
    print("constraints:")
    print("  - pocket:")
    print("      binder: B")
    print(f"      contacts: {yaml_contacts}")
    print("      max_distance: 6.0")
    print("      force: false")

if __name__ == "__main__":
    main()
