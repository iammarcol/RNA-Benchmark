#!/usr/bin/env python3
import os
import subprocess
import csv
import re
import argparse
from collections import defaultdict, Counter

# -------------------------------
# Chain type / SEQRES utilities
# -------------------------------

PROTEIN_RESIDUES = {
    # standard 20
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

NUCLEIC_RESIDUES = {
    # standard
    "A", "C", "G", "U", "I", "T",
    "DA", "DC", "DG", "DT", "DI",
    # modified
    "1MA", "1MG", "2MG", "M2G", "7MG", "MIA",
    "OMG", "OMC", "H2U", "5MC", "PSU", "5MU",
    "H2M", "M5C", "Y", "G7M", "M2A", "M2I",
}


def classify_chain_type_from_resnames(resnames):
    """
    Classify chain as protein / rna / other based on residue names.
    DNA is grouped into 'rna' (nucleic acid) here.
    """
    counts = Counter()
    for rn in resnames:
        rn = rn.upper()
        if rn in PROTEIN_RESIDUES:
            counts["protein"] += 1
        if rn in NUCLEIC_RESIDUES:
            counts["rna"] += 1
    if not counts:
        return "other"
    return counts.most_common(1)[0][0]


def parse_seqres(pdb_path):
    """
    Parse SEQRES records.
    Returns:
        seqres_resnames: dict[chain_id] -> list of residue names (3-letter) in sequence order
        seqres_type: dict[chain_id] -> 'protein' / 'rna' / 'other'
        seqres_length: dict[chain_id] -> integer length from SEQRES
    """
    seqres_resnames = defaultdict(list)

    with open(pdb_path, "r") as fh:
        for line in fh:
            if not line.startswith("SEQRES"):
                continue
            chain_id = line[11].strip() or " "
            parts = line[19:].split()
            seqres_resnames[chain_id].extend(r.upper() for r in parts)

    seqres_type = {}
    seqres_length = {}
    for chain_id, resnames in seqres_resnames.items():
        seqres_type[chain_id] = classify_chain_type_from_resnames(resnames)
        seqres_length[chain_id] = len(resnames)

    return seqres_resnames, seqres_type, seqres_length


def parse_polymer_residues_from_atoms(pdb_path):
    """
    Parse polymer residues from ATOM/HETATM to estimate chain types/lengths
    when SEQRES is missing.
    Returns:
        chain_type: dict[chain_id] -> 'protein'/'rna'/'other'
        chain_length: dict[chain_id] -> number of distinct residues in coordinates
    """
    chain_resnames = defaultdict(list)
    chain_positions = defaultdict(set)

    with open(pdb_path, "r") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            resname = line[17:20].strip().upper()
            chain_id = line[21].strip() or " "
            resseq_str = line[22:26]
            icode = line[26].strip() or " "

            try:
                resseq = int(resseq_str)
            except ValueError:
                continue

            # Only consider polymer residues
            if (resname not in PROTEIN_RESIDUES) and (resname not in NUCLEIC_RESIDUES):
                continue

            chain_resnames[chain_id].append(resname)
            chain_positions[chain_id].add((resseq, icode))

    chain_type = {}
    chain_length = {}
    for chain_id, resnames in chain_resnames.items():
        chain_type[chain_id] = classify_chain_type_from_resnames(resnames)
        chain_length[chain_id] = len(chain_positions[chain_id])

    return chain_type, chain_length


def match_chains(ref_types, ref_lengths, model_types, model_lengths):
    """
    Match reference chains to model chains based on type and length.
    Returns dict: ref_chain -> model_chain

    Strategy:
      - only match chains with same type (protein vs rna)
      - pick model chain with minimum |len_model - len_ref|
      - each model chain is used at most once
    """
    matches = {}
    used_model_chains = set()

    for ref_chain, ref_type in ref_types.items():
        if ref_type == "other":
            continue
        ref_len = ref_lengths.get(ref_chain, 0)
        candidates = []
        for mc, mt in model_types.items():
            if mc in used_model_chains:
                continue
            if mt != ref_type:
                continue
            m_len = model_lengths.get(mc, 0)
            diff = abs(m_len - ref_len)
            candidates.append((diff, mc))

        if not candidates:
            print(f"[WARN] No model chain candidate for ref chain {ref_chain} "
                  f"(type {ref_type}, length {ref_len}).")
            continue

        candidates.sort(key=lambda x: x[0])
        best_diff, best_mc = candidates[0]
        matches[ref_chain] = best_mc
        used_model_chains.add(best_mc)
        print(f"[INFO] Matched ref chain {ref_chain} (len {ref_len}, {ref_type}) "
              f"-> model chain {best_mc} (len {model_lengths.get(best_mc, 0)}, diff {best_diff}).")

    return matches


def build_mapping_string(ref_pdb, model_pdb):
    """
    Build DockQ --mapping string MODELCHAINS:NATIVECHAINS using SEQRES-based
    chain matching (with ATOM fallback).

    Returns mapping string like "CDBA:ABCD" or None if mapping fails.
    """
    # Reference
    ref_seqres_resnames, ref_seqres_type, ref_seqres_len = parse_seqres(ref_pdb)
    ref_atom_type, ref_atom_len = parse_polymer_residues_from_atoms(ref_pdb)

    if ref_seqres_type:
        ref_types = ref_seqres_type
        ref_lengths = ref_seqres_len
        print("[INFO] Using SEQRES for reference chain types/lengths.")
    else:
        ref_types = ref_atom_type
        ref_lengths = ref_atom_len
        print("[INFO] No SEQRES in reference; using ATOM-based types/lengths.")

    # Model
    model_seqres_resnames, model_seqres_type, model_seqres_len = parse_seqres(model_pdb)
    model_atom_type, model_atom_len = parse_polymer_residues_from_atoms(model_pdb)

    if model_seqres_type:
        model_types = model_seqres_type
        model_lengths = model_seqres_len
        print("[INFO] Using SEQRES for model chain types/lengths.")
    else:
        model_types = model_atom_type
        model_lengths = model_atom_len
        print("[INFO] No SEQRES in model; using ATOM-based types/lengths.")

    chain_matches = match_chains(ref_types, ref_lengths, model_types, model_lengths)
    if not chain_matches:
        print("[WARN] No chain matches; not using --mapping.")
        return None

    # preserve reference chain order as encountered
    native_chains = list(chain_matches.keys())
    model_chains = [chain_matches[c] for c in native_chains]

    native_str = "".join(native_chains)
    model_str = "".join(model_chains)

    if not native_str or not model_str or len(native_str) != len(model_str):
        print("[WARN] Invalid chain mapping; not using --mapping.")
        return None

    mapping = f"{model_str}:{native_str}"  # DockQ expects MODELCHAINS:NATIVECHAINS
    print(f"[INFO] DockQ mapping: --mapping {mapping}")
    return mapping


# -------------------------------
# DockQ parsing / running
# -------------------------------

def save_dockq_output(output, output_path):
    """Save the full DockQ output to a file."""
    with open(output_path, "w") as f:
        f.write(output)


def parse_dockq_output(file_path, pdb_id):
    """Parse the saved DockQ output file to extract DockQ scores, F1 scores, and chain pairs."""
    results = []
    with open(file_path, "r") as f:
        lines = f.readlines()
    
    # Extract total DockQ score
    total_dockq_match = re.search(r"Total DockQ over \d+ native interfaces: ([\d.]+)", ''.join(lines))
    total_dockq_score = float(total_dockq_match.group(1)) if total_dockq_match else None
    
    for i, line in enumerate(lines):
        if line.startswith("Native chains:"):
            native_chains = line.split(":", 1)[1].strip() if ":" in line else ""
            model_chains_line = lines[i + 1] if i + 1 < len(lines) else ""
            model_chains = model_chains_line.split(":", 1)[1].strip() if ":" in model_chains_line else ""
            
            dockq_match = re.search(r"DockQ:\s([\d.]+)", lines[i + 2] if i + 2 < len(lines) else "")
            dockq_score = float(dockq_match.group(1)) if dockq_match else None
            
            f1_score = None
            for j in range(i, min(i + 10, len(lines))):
                f1_match = re.search(r"F1:\s([\d.]+)", lines[j])
                if f1_match:
                    f1_score = float(f1_match.group(1))
                    break
            
            results.append({
                "PDB_ID": pdb_id,
                "Total_DockQ_Score": total_dockq_score,
                "DockQ_Score": dockq_score,
                "F1_Score": f1_score,
                "Native_chains": native_chains,
                "Model_chains": model_chains
            })
    return results


def run_dockq_once(model_path, reference_path, output_dir, pdb_id, mapping=None, suffix="nomap"):
    """
    Run DockQ a single time (optionally with --mapping).
    Returns (success_bool, output_path).
    """
    try:
        command = ["DockQ", model_path, reference_path]
        if mapping:
            command.extend(["--mapping", mapping])

        print(f"[INFO] Running DockQ ({suffix}) for {pdb_id}: {' '.join(command)}")

        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        
        os.makedirs(output_dir, exist_ok=True)
        dockq_output_path = os.path.join(output_dir, f"{pdb_id}_dockq_{suffix}.txt")
        save_dockq_output(result.stdout + "\n\nSTDERR:\n" + result.stderr, dockq_output_path)

        if result.returncode != 0:
            print(f"[WARN] DockQ ({suffix}) returned non-zero exit code for {pdb_id}. "
                  f"Return code: {result.returncode}")
            return False, dockq_output_path

        return True, dockq_output_path

    except FileNotFoundError:
        print("DockQ is not installed or not in the PATH.")
        return False, None


def run_dockq_with_fallback(model_path, reference_path, output_dir, pdb_id):
    """
    Try DockQ without mapping first.
    If that fails, build mapping and retry with --mapping.
    Returns path to the successful DockQ output file, or None if all attempts fail.
    """
    # 1) Try without mapping
    success, output_path = run_dockq_once(
        model_path, reference_path, output_dir, pdb_id,
        mapping=None, suffix="nomap"
    )
    if success:
        print(f"[INFO] DockQ succeeded without mapping for {pdb_id}.")
        return output_path

    print(f"[INFO] DockQ without mapping failed for {pdb_id}, trying with --mapping...")

    # 2) Build mapping and retry
    mapping = build_mapping_string(reference_path, model_path)
    if not mapping:
        print(f"[ERROR] Could not build a valid mapping for {pdb_id}. Skipping.")
        return None

    success, output_path = run_dockq_once(
        model_path, reference_path, output_dir, pdb_id,
        mapping=mapping, suffix="mapping"
    )
    if success:
        print(f"[INFO] DockQ succeeded with mapping for {pdb_id}.")
        return output_path

    print(f"[ERROR] DockQ failed even with mapping for {pdb_id}.")
    return None


# -------------------------------
# Main driver
# -------------------------------

def main(model_dir, reference_dir, output_csv, output_dir):
    """Main function to run DockQ on matched model-reference pairs."""
    model_files = {
        os.path.splitext(f)[0]: os.path.join(model_dir, f)
        for f in os.listdir(model_dir)
        if f.endswith(".pdb")
    }
    reference_files = {
        os.path.splitext(f)[0]: os.path.join(reference_dir, f)
        for f in os.listdir(reference_dir)
        if f.endswith(".pdb")
    }
    
    matched_ids = set(model_files.keys()).intersection(reference_files.keys())
    
    if not matched_ids:
        print("No matched PDB IDs found between the two directories.")
        return

    print(f"Found {len(matched_ids)} matched PDB files. Running DockQ...")

    csv_dir = os.path.dirname(output_csv)
    if csv_dir:
        os.makedirs(csv_dir, exist_ok=True)

    with open(output_csv, mode='w', newline='') as csvfile:
        fieldnames = ["PDB_ID", "Total_DockQ_Score", "DockQ_Score",
                      "F1_Score", "Native_chains", "Model_chains"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for pdb_id in sorted(matched_ids):
            model_path = model_files[pdb_id]
            reference_path = reference_files[pdb_id]
            
            print(f"\n[INFO] Processing PDB_ID: {pdb_id}")
            dockq_output_path = run_dockq_with_fallback(
                model_path, reference_path, output_dir, pdb_id
            )
            
            if dockq_output_path:
                parsed_results = parse_dockq_output(dockq_output_path, pdb_id)
                if parsed_results:
                    writer.writerows(parsed_results)
                else:
                    print(f"[WARN] No DockQ results parsed for {pdb_id}.")
            else:
                print(f"[ERROR] DockQ completely failed for {pdb_id}; skipping.")

    print(f"\nDockQ results saved to {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run DockQ on matched PDB files and extract scores. "
                    "Tries without --mapping first, then with custom --mapping if needed."
    )
    parser.add_argument("--model_dir", required=True, help="Directory containing model PDB files.")
    parser.add_argument("--reference_dir", required=True, help="Directory containing reference PDB files.")
    parser.add_argument("--output_csv", required=True, help="Path to output CSV file.")
    parser.add_argument("--output_dir", required=True, help="Directory to store full DockQ outputs.")
    
    args = parser.parse_args()
    main(args.model_dir, args.reference_dir, args.output_csv, args.output_dir)
