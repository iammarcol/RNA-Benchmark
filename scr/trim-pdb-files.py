#!/usr/bin/env python3
import os
import argparse
from collections import defaultdict, Counter

######### res description #########

PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

NUCLEIC_RESIDUES = {
    "A", "C", "G", "U", "I", "T",
    "DA", "DC", "DG", "DT", "DI",
    "1MA", "1MG", "2MG", "M2G", "7MG", "MIA",
    "OMG", "OMC", "H2U", "5MC", "PSU", "5MU",
    "H2M", "M5C", "Y", "G7M", "M2A", "M2I",
}

######### MAPPING #########
CANONICAL_MAP = {
    # prot
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",

    # rna
    "A": "A", "C": "C", "G": "G", "U": "U", "I": "I",

    # dna
    "DA": "A",
    "DC": "C",
    "DG": "G",
    "DT": "U",
    "DI": "I",

    # modified nts
    "1MA": "A",
    "1MG": "G",
    "2MG": "G",
    "M2G": "G",
    "7MG": "G",
    "MIA": "A",
    "OMG": "G",
    "OMC": "C",
    "H2U": "U",
    "5MC": "C",
    "PSU": "U",
    "5MU": "U",
    "Y":   "U",
    "G7M": "G",
    "M2A": "A",
    "M2I": "I",
}

def canonical_residue(resname: str) -> str:
    resname = resname.upper()
    return CANONICAL_MAP.get(resname, "X")


def classify_chain_type_from_resnames(resnames):
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



######### SEQRES parsing #########

def parse_seqres(pdb_path):
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


# --- ATOM/HETATM parsing ---

def parse_polymer_residues_from_atoms(pdb_path):
    chain_residues = defaultdict(list)
    chain_resnames = defaultdict(list)

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

            canon = canonical_residue(resname)
            if canon == "X":
                continue  # skip unknown residues

            chain_residues[chain_id].append({
                "resseq": resseq,
                "icode": icode,
                "resname": resname,
                "canonical": canon,
            })
            chain_resnames[chain_id].append(resname)

    chain_type = {}
    chain_length = {}
    for chain_id, residues in chain_residues.items():
        chain_type[chain_id] = classify_chain_type_from_resnames(chain_resnames[chain_id])
        unique_positions = {(r["resseq"], r["icode"]) for r in residues}
        chain_length[chain_id] = len(unique_positions)

    return chain_residues, chain_type, chain_length


def unique_residue_list(residue_list):
    seen = set()
    unique = []
    for r in residue_list:
        key = (r["resseq"], r["icode"])
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique


def get_chain_sequence(chain_residues, chain_id):
    """Return canonical one-letter sequence for a chain from residue dict."""
    res_list = unique_residue_list(chain_residues.get(chain_id, []))
    return [r["canonical"] for r in res_list]


def alignment_identity(seq1, seq2):
    """Compute fraction of identical matches in global alignment."""
    if not seq1 or not seq2:
        return 0.0

    aln = glocal_affine_align_indices(seq1, seq2)
    matches = 0
    aligned = 0
    for i_idx, j_idx in aln:
        if i_idx is None or j_idx is None:
            continue
        aligned += 1
        if seq1[i_idx] == seq2[j_idx]:
            matches += 1
    return matches / aligned if aligned > 0 else 0.0




def glocal_affine_align_indices(seq1, seq2,
                                match_score=2,
                                mismatch_score=-3,
                                gap_open=-5,
                                gap_extend=-1):
    """
    Global in seq1 (ref), local in seq2 (model) with affine gaps.

    Returns:
        alignment: list of (i_idx, j_idx) pairs as before
                   (i_idx or j_idx can be None for gaps).
    """
    n, m = len(seq1), len(seq2)
    NEG_INF = -10**9

    # M: match/mismatch; Ix: gap in seq1; Iy: gap in seq2
    M  = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Iy = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    # pointers: (state, prev_i, prev_j), state: 0=M, 1=Ix, 2=Iy
    ptrM  = [[None] * (m + 1) for _ in range(n + 1)]
    ptrIx = [[None] * (m + 1) for _ in range(n + 1)]
    ptrIy = [[None] * (m + 1) for _ in range(n + 1)]

    # initialization
    M[0][0] = 0
    Ix[0][0] = Iy[0][0] = NEG_INF

    # i>0, j=0: align all of seq1 => penalized gaps in seq2
    for i in range(1, n + 1):
        Ix[i][0] = gap_open + (i - 1) * gap_extend
        ptrIx[i][0] = (1, i - 1, 0)  # extend gap in seq2
        # M[i][0], Iy[i][0] stay -inf

    # i=0, j>0: local in seq2 => free leading gaps in seq2
    for j in range(1, m + 1):
        M[0][j] = 0
        ptrM[0][j] = (0, 0, j - 1)

    # fill
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # M: match/mismatch
            s = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            candM = [
                (M[i - 1][j - 1] + s,  (0, i - 1, j - 1)),
                (Ix[i - 1][j - 1] + s, (1, i - 1, j - 1)),
                (Iy[i - 1][j - 1] + s, (2, i - 1, j - 1)),
            ]
            best_val, best_ptr = max(candM, key=lambda x: x[0])
            M[i][j] = best_val
            ptrM[i][j] = best_ptr

            # Ix: gap in seq1 (insert in seq2)
            candIx = [
                (M[i - 1][j] + gap_open + gap_extend, (0, i - 1, j)),
                (Ix[i - 1][j] + gap_extend,           (1, i - 1, j)),
            ]
            best_val, best_ptr = max(candIx, key=lambda x: x[0])
            Ix[i][j] = best_val
            ptrIx[i][j] = best_ptr

            # Iy: gap in seq2 (insert in seq1)
            candIy = [
                (M[i][j - 1] + gap_open + gap_extend, (0, i, j - 1)),
                (Iy[i][j - 1] + gap_extend,           (2, i, j - 1)),
            ]
            best_val, best_ptr = max(candIy, key=lambda x: x[0])
            Iy[i][j] = best_val
            ptrIy[i][j] = best_ptr

    # finish: i = n, best over all j and states (local in seq2)
    best_score = NEG_INF
    best_state = 0
    best_j = 0
    for j in range(0, m + 1):
        for state, mat in enumerate((M, Ix, Iy)):
            if mat[n][j] > best_score or (
                mat[n][j] == best_score and j > best_j
            ):
                best_score = mat[n][j]
                best_state = state
                best_j = j

    # traceback
    alignment = []
    state = best_state
    i, j = n, best_j

    while i > 0 or j > 0:
        if i == 0 and j > 0:
            # leading/trailing gaps in seq2 are free
            alignment.append((None, j - 1))
            j -= 1
            continue
        if j == 0 and i > 0:
            alignment.append((i - 1, None))
            i -= 1
            continue

        if state == 0:
            prev = ptrM[i][j]
        elif state == 1:
            prev = ptrIx[i][j]
        else:
            prev = ptrIy[i][j]

        if prev is None:
            break

        prev_state, pi, pj = prev

        if pi == i - 1 and pj == j - 1:
            alignment.append((i - 1, j - 1))
        elif pi == i - 1 and pj == j:
            alignment.append((i - 1, None))    # gap in seq2
        elif pi == i and pj == j - 1:
            alignment.append((None, j - 1))    # gap in seq1

        state = prev_state
        i, j = pi, pj

    alignment.reverse()
    return alignment



######### chain matching #########

def match_chains(ref_types, ref_lengths, model_types, model_lengths,
                 ref_atoms, model_atoms):
    """
    match reference chains to model chains using both type and sequence similarity.
    returns dict: ref_chain -> model_chain
    """
    matches = {}
    used = set()

    # Optional: process longer ref chains first
    ref_chains_sorted = sorted(
        ref_types.keys(),
        key=lambda c: ref_lengths.get(c, 0),
        reverse=True
    )

    for ref_chain in ref_chains_sorted:
        ref_type = ref_types[ref_chain]
        if ref_type == "other":
            continue

        ref_seq = get_chain_sequence(ref_atoms, ref_chain)
        if not ref_seq:
            # fall back to length-only behavior if no sequence
            print(f"no residue sequence for ref chain {ref_chain}, skipping.")
            continue

        ref_len = len(ref_seq)

        best_candidate = None  # (identity, -len_diff, model_chain)
        for mc, mt in model_types.items():
            if mc in used:
                continue
            if mt != ref_type:
                continue

            model_seq = get_chain_sequence(model_atoms, mc)
            if not model_seq:
                continue

            model_len = len(model_seq)
            len_diff = abs(model_len - ref_len)

            ident = alignment_identity(ref_seq, model_seq)
            if ident <= 0.0:
                continue

            # higher identity is better; for tie, smaller length difference better
            candidate = (ident, -len_diff, mc)
            if best_candidate is None or candidate > best_candidate:
                best_candidate = candidate

        if best_candidate is None:
            print(f"no chain candidate (with sequence similarity) for {ref_chain}")
            continue

        ident, neg_len_diff, best_mc = best_candidate
        matches[ref_chain] = best_mc
        used.add(best_mc)
        print(f"matched {ref_chain} → {best_mc} (identity={ident:.3f})")

    return matches



######### Build allowed residues and optionally write alignments #########

def build_allowed_residues(ref_chain_residues,
                           model_chain_residues,
                           chain_matches,
                           aln_dir=None,
                           pid=None):
    """
    Build set of allowed residues in model based on alignment.
    If aln_dir is not None, write a FASTA alignment file per chain pair.
    """
    allowed = set()

    for ref_chain, model_chain in chain_matches.items():
        ref_list = unique_residue_list(ref_chain_residues.get(ref_chain, []))
        model_list = unique_residue_list(model_chain_residues.get(model_chain, []))

        ref_seq = [r["canonical"] for r in ref_list]
        model_seq = [r["canonical"] for r in model_list]

        if not ref_seq or not model_seq:
            continue

        alignment = glocal_affine_align_indices(ref_seq, model_seq)

        # reconstruct aligned strings + track which model indices to keep
        ref_aln = []
        model_aln = []
        keep_j = set()

        for i_idx, j_idx in alignment:
            # build aligned ref string
            if i_idx is None:
                ref_aln.append("-")
            else:
                ref_aln.append(ref_seq[i_idx])

            # build aligned model string
            if j_idx is None:
                model_aln.append("-")
            else:
                model_aln.append(model_seq[j_idx])

            # mark kept residues: only when both present and equal
            if (
                i_idx is not None
                and j_idx is not None
                and ref_seq[i_idx] == model_seq[j_idx]
            ):
                keep_j.add(j_idx)

        if not keep_j:
            print(f"[WARN] No matching residues in {ref_chain} → {model_chain}")
        else:
            for j in sorted(keep_j):
                r = model_list[j]
                allowed.add((model_chain, r["resseq"], r["icode"], r["canonical"]))

            print(f"[INFO] Keeping {len(keep_j)}/{len(model_list)} residues")

        # write alignment file if requested
        if aln_dir is not None and pid is not None:
            fname = f"{pid}_ref{ref_chain}_model{model_chain}.aln.fasta"
            aln_path = os.path.join(aln_dir, fname)

            ref_header = f">ref {pid}_{ref_chain}"
            model_header = f">model {pid}_{model_chain}"

            with open(aln_path, "w") as f:
                f.write(ref_header + "\n")
                f.write("".join(ref_aln) + "\n")
                f.write(model_header + "\n")
                f.write("".join(model_aln) + "\n")

            print(f"[INFO] Wrote alignment: {aln_path}")

    return allowed



######### Trim model PDB #########

def trim_model_pdb(model_path, out_path, allowed_set):
    with open(model_path, "r") as inp, open(out_path, "w") as out:
        atom_serial = 1
        for line in inp:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                out.write(line)
                continue

            resname = line[17:20].strip().upper()
            chain_id = line[21].strip() or " "
            resseq_str = line[22:26]
            icode = line[26].strip() or " "

            try:
                resseq = int(resseq_str)
            except ValueError:
                out.write(f"{rec:<6}{atom_serial:>5}" + line[11:])
                atom_serial += 1
                continue

            canon = canonical_residue(resname)
            key = (chain_id, resseq, icode, canon)

            if canon != "X" and key not in allowed_set:
                continue  # trim residue

            out.write(f"{rec:<6}{atom_serial:>5}" + line[11:])
            atom_serial += 1




###########################

def main():
    parser = argparse.ArgumentParser(
        description="Trim model PDB files based on reference PDB SEQRES + coordinate alignment."
    )
    parser.add_argument("--ref_dir", required=True)
    parser.add_argument("--model_dir", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--suffix", default="")
    parser.add_argument(
        "--aln_dir",
        help="Optional directory to write FASTA alignment files per chain pair.",
        default=None,
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    if args.aln_dir is not None:
        os.makedirs(args.aln_dir, exist_ok=True)

    ref_ids = {
        os.path.splitext(f)[0]
        for f in os.listdir(args.ref_dir)
        if f.lower().endswith(".pdb")
    }
    model_ids = {
        os.path.splitext(f)[0]
        for f in os.listdir(args.model_dir)
        if f.lower().endswith(".pdb")
    }
    common = sorted(ref_ids & model_ids)

    print(f"[INFO] {len(common)} common IDs.")

    for pid in common:
        print(f"\n[INFO] Processing {pid}")

        ref_path = os.path.join(args.ref_dir, pid + ".pdb")
        model_path = os.path.join(args.model_dir, pid + ".pdb")

        try:
            # parse reference
            ref_seqres, ref_seqtype, ref_seqlen = parse_seqres(ref_path)
            ref_atoms, ref_atomtype, ref_atomlen = parse_polymer_residues_from_atoms(ref_path)

            # parse model
            model_seqres, model_seqtype, model_seqlen = parse_seqres(model_path)
            model_atoms, model_atomtype, model_atomlen = parse_polymer_residues_from_atoms(model_path)

            ref_types = ref_seqtype if ref_seqtype else ref_atomtype
            ref_lengths = ref_seqlen if ref_seqtype else ref_atomlen
            model_types = model_seqtype if model_seqtype else model_atomtype
            model_lengths = model_seqlen if model_seqtype else model_atomlen

            chain_matches = match_chains(
                ref_types, ref_lengths, model_types, model_lengths,
                ref_atoms, model_atoms
            )
            if not chain_matches:
                continue

            allowed = build_allowed_residues(
                ref_atoms,
                model_atoms,
                chain_matches,
                aln_dir=args.aln_dir,
                pid=pid,
            )
            if not allowed:
                print(f"[WARN] No allowed residues for {pid}, skipping.")
                continue

            out_path = os.path.join(args.out_dir, pid + args.suffix + ".pdb")
            trim_model_pdb(model_path, out_path, allowed)
            print(f"[INFO] Wrote: {out_path}")

        except Exception as e:
            print(f"[ERROR] Failed {pid}: {e}")


if __name__ == "__main__":
    main()

# run as: 
# python trim-pdb-files.py --ref_dir --model_dir --out_dir
