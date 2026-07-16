#!/usr/bin/env python3

import argparse
import os
import warnings

import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
from Bio import pairwise2


RNA_CODES = {
    "A": "A", "U": "U", "G": "G", "C": "C",
    "RA": "A", "RU": "U", "RG": "G", "RC": "C",
    "ADE": "A", "URA": "U", "GUA": "G", "CYT": "C",
    "DA": "A", "DU": "U", "DG": "G", "DC": "C","DT": "T",

}

DNA_CODES = {
    "DA": "A", "DT": "T", "DG": "G", "DC": "C",
    "ADE": "A", "THY": "T", "GUA": "G", "CYT": "C",
}

PROTEIN_CODES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

VALID_TYPES = {"protein", "rna", "dna"}


def load_structure(path):
    if path.endswith(".cif") or path.endswith(".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return parser.get_structure(os.path.basename(path), path)


def parse_list_field(x):
    """
    Parses ex. "A, C" into ["A", "C"]. Bc this uses an external CSV where chain types are annotated already!
    """
    if pd.isna(x):
        return []

    return [v.strip().strip('"').strip("'") for v in str(x).split(",")]


def normalize_type(x):
    x = str(x).strip().lower()

    if x in {"protein", "prot", "peptide"}:
        return "protein"

    if x in {"rna", "ribonucleic_acid", "ribonucleic acid"}:
        return "rna"

    if x in {"dna", "deoxyribonucleic_acid", "deoxyribonucleic acid"}:
        return "dna"

    return x


def residue_key(residue):
    chain_id = residue.get_parent().id
    hetflag, resseq, icode = residue.id
    return (chain_id, int(resseq), icode.strip())


def residue_label(key):
    chain_id, resseq, icode = key
    return f"{chain_id}:{resseq}{icode}"


def residue_name(residue):
    return residue.get_resname().strip()


def is_rna_residue(residue):
    return residue_name(residue) in RNA_CODES


def is_dna_residue(residue):
    return residue_name(residue) in DNA_CODES


def is_protein_residue(residue):
    return residue_name(residue) in PROTEIN_CODES


def residue_to_letter(residue):
    name = residue_name(residue)

    if name in RNA_CODES:
        return RNA_CODES[name]

    if name in DNA_CODES:
        return DNA_CODES[name]

    if name in PROTEIN_CODES:
        return PROTEIN_CODES[name]

    return "X"


def get_chain_residues(structure, chain_id, molecule_type):
    try:
        chain = structure[0][chain_id]
    except KeyError:
        raise KeyError(f"Chain {chain_id} not found in {structure.id}")

    residues = []

    for res in chain:
        if molecule_type == "rna" and is_rna_residue(res):
            residues.append(res)
        elif molecule_type == "dna" and is_dna_residue(res):
            residues.append(res)
        elif molecule_type == "protein" and is_protein_residue(res):
            residues.append(res)

    return residues


def get_sequence_and_residues(structure, chain_id, molecule_type):
    residues = get_chain_residues(structure, chain_id, molecule_type)
    seq = "".join(residue_to_letter(r) for r in residues)
    return seq, residues


def map_residues_by_sequence(ref_seq, ref_residues, model_seq, model_residues):
    """
    Maps model residue keys onto reference residue keys using global sequence alignment.

    Returns a model_to_ref dict
    key in model structure -> corresponding key in reference structure
    """
    if len(ref_seq) == 0 or len(model_seq) == 0:
        return {}

    aln = pairwise2.align.globalms(
        ref_seq,
        model_seq,
        2,      # match
        -1,     # mismatch
        -10,    # gap open
        -1,     # gap extend
        one_alignment_only=True,
    )[0]

    ref_aln = aln.seqA
    model_aln = aln.seqB

    ref_i = 0
    model_i = 0
    model_to_ref = {}

    for a, b in zip(ref_aln, model_aln):
        ref_res = None
        model_res = None

        if a != "-":
            ref_res = ref_residues[ref_i]
            ref_i += 1

        if b != "-":
            model_res = model_residues[model_i]
            model_i += 1

        if ref_res is not None and model_res is not None:
            model_to_ref[residue_key(model_res)] = residue_key(ref_res)

    return model_to_ref


def collect_residues_by_chain_type(structure, chain_type_pairs):
    """
    chain_type_pairs is a list of (chain_id, molecule_type)
    Returns all residues from those chains matching their declared type.
    """
    residues = []

    for chain_id, molecule_type in chain_type_pairs:
        residues.extend(get_chain_residues(structure, chain_id, molecule_type))

    return residues


def find_interface_contacts_between_groups(
    structure,
    group1_chain_type_pairs,
    group2_chain_type_pairs,
    cutoff,
):
    """
    Finds residue-level contacts between two chain groups using any-atom distance.
    Atom distance is defined (flag --cutoff)

    Returns:
        group1_interface: set of residue keys from group 1
        group2_interface: set of residue keys from group 2
        contacts: set of (group1_residue_key, group2_residue_key)
    """
    group1_residues = collect_residues_by_chain_type(
        structure,
        group1_chain_type_pairs,
    )

    group2_residues = collect_residues_by_chain_type(
        structure,
        group2_chain_type_pairs,
    )

    group2_atoms = []
    atom_to_group2_residue = {}

    for res2 in group2_residues:
        for atom in res2.get_atoms():
            group2_atoms.append(atom)
            atom_to_group2_residue[atom] = res2

    if len(group1_residues) == 0 or len(group2_atoms) == 0:
        return set(), set(), set()

    ns = NeighborSearch(group2_atoms)

    group1_interface = set()
    group2_interface = set()
    contacts = set()

    for res1 in group1_residues:
        key1 = residue_key(res1)

        for atom1 in res1.get_atoms():
            nearby_atoms = ns.search(atom1.coord, cutoff, level="A")

            for atom2 in nearby_atoms:
                res2 = atom_to_group2_residue[atom2]
                key2 = residue_key(res2)

                group1_interface.add(key1)
                group2_interface.add(key2)
                contacts.add((key1, key2))

    return group1_interface, group2_interface, contacts


def remap_residue_set(model_set, model_to_ref):
    return {model_to_ref[x] for x in model_set if x in model_to_ref}


def remap_contacts(model_contacts, group1_model_to_ref, group2_model_to_ref):
    mapped_contacts = set()

    for key1, key2 in model_contacts:
        if key1 in group1_model_to_ref and key2 in group2_model_to_ref:
            mapped_contacts.add(
                (
                    group1_model_to_ref[key1],
                    group2_model_to_ref[key2],
                )
            )

    return mapped_contacts


def precision_recall_f1(pred, ref):
    pred = set(pred)
    ref = set(ref)

    tp = len(pred & ref)
    fp = len(pred - ref)
    fn = len(ref - pred)

    precision = tp / (tp + fp) if tp + fp > 0 else 0.0
    recall = tp / (tp + fn) if tp + fn > 0 else 0.0

    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall > 0
        else 0.0
    )

    return tp, fp, fn, precision, recall, f1


def interface_sequence_from_ref(ref_structure, ref_interface_keys, ref_chain_type_pairs):
    """
    Returns interface sequence ordered according to the reference chain residue order.
    """
    ordered_letters = []
    ordered_labels = []

    for chain_id, molecule_type in ref_chain_type_pairs:
        residues = get_chain_residues(ref_structure, chain_id, molecule_type)

        for res in residues:
            key = residue_key(res)

            if key in ref_interface_keys:
                ordered_letters.append(residue_to_letter(res))
                ordered_labels.append(residue_label(key))

    return "".join(ordered_letters), ";".join(ordered_labels)


def build_chain_mapping_from_row(row):
    native_chains = parse_list_field(row["Native_chains"])
    model_chains = parse_list_field(row["Model_chains"])
    chain_types = [normalize_type(x) for x in parse_list_field(row["Type"])]

    if not (len(native_chains) == len(model_chains) == len(chain_types)):
        raise ValueError(
            f"Length mismatch: Native_chains={native_chains}, "
            f"Model_chains={model_chains}, Type={chain_types}"
        )

    mappings = []

    for native_chain, model_chain, molecule_type in zip(
        native_chains,
        model_chains,
        chain_types,
    ):
        if molecule_type not in VALID_TYPES:
            raise ValueError(f"Unsupported molecule type: {molecule_type}")

        mappings.append(
            {
                "native_chain": native_chain,
                "model_chain": model_chain,
                "type": molecule_type,
            }
        )

    return mappings


def build_rna_centered_groups(mappings):
    """
    Builds RNA-centered interface groups.

    Cases:
        RNA-protein:
            group1 = all RNA chains
            group2 = all protein chains

        RNA-DNA:
            group1 = all RNA chains
            group2 = all DNA chains

        RNA-protein-DNA:
            group1 = all RNA chains
            group2 = all non-RNA chains

        RNA-RNA:
            group1 = first RNA chain
            group2 = remaining RNA chains

    Rows without RNA are skipped before this function.
    """
    rna_mappings = [m for m in mappings if m["type"] == "rna"]
    partner_mappings = [m for m in mappings if m["type"] != "rna"]

    if len(rna_mappings) == 0:
        raise ValueError("This row does not contain RNA.")

    if len(partner_mappings) == 0:
        # RNA-RNA case.
        if len(rna_mappings) < 2:
            raise ValueError("RNA-only row must contain at least two RNA chains.")

        group1 = rna_mappings[:1]
        group2 = rna_mappings[1:]
        interface_type = "rna-rna"
    else:
        # RNA-protein, RNA-DNA, or RNA with multiple partner types.
        group1 = rna_mappings
        group2 = partner_mappings

        partner_types = sorted({m["type"] for m in partner_mappings})
        interface_type = "rna-" + "+".join(partner_types)

    return group1, group2, interface_type


def mapping_to_chain_type_pairs(mappings, chain_key):
    """
    chain_key should be either:
        "native_chain"
        "model_chain"
    """
    return [(m[chain_key], m["type"]) for m in mappings]


def build_model_to_ref_maps(group_mappings, ref_structure, model_structure):
    """
    Builds residue mapping from model residue keys to reference residue keys
    for all chains in a group.
    """
    model_to_ref_all = {}

    for m in group_mappings:
        native_chain = m["native_chain"]
        model_chain = m["model_chain"]
        molecule_type = m["type"]

        ref_seq, ref_residues = get_sequence_and_residues(
            ref_structure,
            native_chain,
            molecule_type,
        )

        model_seq, model_residues = get_sequence_and_residues(
            model_structure,
            model_chain,
            molecule_type,
        )

        model_to_ref = map_residues_by_sequence(
            ref_seq,
            ref_residues,
            model_seq,
            model_residues,
        )

        model_to_ref_all.update(model_to_ref)

    return model_to_ref_all


def compare_one_interface(row, ref_structure, model_structure, cutoff):
    mappings = build_chain_mapping_from_row(row)

    group1_mappings, group2_mappings, interface_type = build_rna_centered_groups(mappings)

    ref_group1_pairs = mapping_to_chain_type_pairs(group1_mappings, "native_chain")
    ref_group2_pairs = mapping_to_chain_type_pairs(group2_mappings, "native_chain")

    model_group1_pairs = mapping_to_chain_type_pairs(group1_mappings, "model_chain")
    model_group2_pairs = mapping_to_chain_type_pairs(group2_mappings, "model_chain")

    ref_rna_interface, ref_partner_interface, ref_contacts = find_interface_contacts_between_groups(
        ref_structure,
        ref_group1_pairs,
        ref_group2_pairs,
        cutoff,
    )

    model_rna_interface, model_partner_interface, model_contacts = find_interface_contacts_between_groups(
        model_structure,
        model_group1_pairs,
        model_group2_pairs,
        cutoff,
    )

    rna_model_to_ref = build_model_to_ref_maps(
        group1_mappings,
        ref_structure,
        model_structure,
    )

    partner_model_to_ref = build_model_to_ref_maps(
        group2_mappings,
        ref_structure,
        model_structure,
    )

    mapped_model_rna_interface = remap_residue_set(
        model_rna_interface,
        rna_model_to_ref,
    )

    mapped_model_partner_interface = remap_residue_set(
        model_partner_interface,
        partner_model_to_ref,
    )

    mapped_model_contacts = remap_contacts(
        model_contacts,
        rna_model_to_ref,
        partner_model_to_ref,
    )

    r_tp, r_fp, r_fn, r_prec, r_rec, r_f1 = precision_recall_f1(
        mapped_model_rna_interface,
        ref_rna_interface,
    )

    p_tp, p_fp, p_fn, p_prec, p_rec, p_f1 = precision_recall_f1(
        mapped_model_partner_interface,
        ref_partner_interface,
    )

    c_tp, c_fp, c_fn, c_prec, c_rec, c_f1 = precision_recall_f1(
        mapped_model_contacts,
        ref_contacts,
    )

    ref_rna_seq, ref_rna_labels = interface_sequence_from_ref(
        ref_structure,
        ref_rna_interface,
        ref_group1_pairs,
    )

    model_rna_seq_mapped, model_rna_labels_mapped = interface_sequence_from_ref(
        ref_structure,
        mapped_model_rna_interface,
        ref_group1_pairs,
    )

    ref_partner_seq, ref_partner_labels = interface_sequence_from_ref(
        ref_structure,
        ref_partner_interface,
        ref_group2_pairs,
    )

    model_partner_seq_mapped, model_partner_labels_mapped = interface_sequence_from_ref(
        ref_structure,
        mapped_model_partner_interface,
        ref_group2_pairs,
    )

    return {
        "interface_type": interface_type,

        "ref_rna_interface_n": len(ref_rna_interface),
        "model_rna_interface_n_mapped": len(mapped_model_rna_interface),
        "rna_interface_tp": r_tp,
        "rna_interface_fp": r_fp,
        "rna_interface_fn": r_fn,
        "rna_interface_precision": r_prec,
        "rna_interface_recall": r_rec,
        "rna_interface_f1": r_f1,

        "ref_partner_interface_n": len(ref_partner_interface),
        "model_partner_interface_n_mapped": len(mapped_model_partner_interface),
        "partner_interface_tp": p_tp,
        "partner_interface_fp": p_fp,
        "partner_interface_fn": p_fn,
        "partner_interface_precision": p_prec,
        "partner_interface_recall": p_rec,
        "partner_interface_f1": p_f1,

        "ref_contact_n": len(ref_contacts),
        "model_contact_n_mapped": len(mapped_model_contacts),
        "contact_tp": c_tp,
        "contact_fp": c_fp,
        "contact_fn": c_fn,
        "contact_precision": c_prec,
        "contact_recall": c_rec,
        "contact_f1": c_f1,

        "ref_rna_interface_sequence": ref_rna_seq,
        "model_rna_interface_sequence_mapped_to_ref": model_rna_seq_mapped,
        "ref_rna_interface_residues": ref_rna_labels,
        "model_rna_interface_residues_mapped_to_ref": model_rna_labels_mapped,

        "ref_partner_interface_sequence": ref_partner_seq,
        "model_partner_interface_sequence_mapped_to_ref": model_partner_seq_mapped,
        "ref_partner_interface_residues": ref_partner_labels,
        "model_partner_interface_residues_mapped_to_ref": model_partner_labels_mapped,
    }


def find_structure_file(directory, pdb_id):
    candidates = [
        os.path.join(directory, f"{pdb_id}.pdb"),
        os.path.join(directory, f"{pdb_id.lower()}.pdb"),
        os.path.join(directory, f"{pdb_id.upper()}.pdb"),

        os.path.join(directory, f"{pdb_id}.cif"),
        os.path.join(directory, f"{pdb_id.lower()}.cif"),
        os.path.join(directory, f"{pdb_id.upper()}.cif"),

        os.path.join(directory, f"{pdb_id}.mmcif"),
        os.path.join(directory, f"{pdb_id.lower()}.mmcif"),
        os.path.join(directory, f"{pdb_id.upper()}.mmcif"),
    ]

    for path in candidates:
        if os.path.exists(path):
            return path

    return None


def row_has_rna(chain_types):
    return "rna" in chain_types


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Batch compare RNA-centered interface residue and contact agreement. "
            "Processes only rows where Type contains RNA. Supports RNA-protein, "
            "RNA-DNA, and RNA-RNA interfaces."
        )
    )

    parser.add_argument(
        "--model-dir",
        required=True,
        help="Directory containing predicted model files named {PDB_ID}.pdb/.cif/.mmcif",
    )

    parser.add_argument(
        "--ref-dir",
        required=True,
        help="Directory containing reference/native files named {PDB_ID}.pdb/.cif/.mmcif",
    )

    parser.add_argument(
        "--dockq-csv",
        required=True,
        help="DockQ CSV with PDB_ID, Native_chains, Model_chains, and Type columns",
    )

    parser.add_argument(
        "--out-csv",
        required=True,
        help="Output CSV path",
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        default=5.0,
        help="Any-atom interface contact cutoff in Angstrom. Default: 5.0",
    )

    parser.add_argument(
        "--only-rna",
        action="store_true",
        help=(
            "Compatibility flag. The script already only processes rows whose "
            "Type contains RNA."
        ),
    )

    args = parser.parse_args()

    df = pd.read_csv(args.dockq_csv)

    required_cols = {"PDB_ID", "Native_chains", "Model_chains", "Type"}
    missing = required_cols - set(df.columns)

    if missing:
        raise ValueError(f"Missing required columns in DockQ CSV: {missing}")

    structure_cache = {}
    results = []

    for idx, row in df.iterrows():
        pdb_id = str(row["PDB_ID"])

        chain_types = [normalize_type(x) for x in parse_list_field(row["Type"])]

        # Critical rule:
        # skip every row unless Type contains RNA.
        if not row_has_rna(chain_types):
            continue

        ref_path = find_structure_file(args.ref_dir, pdb_id)
        model_path = find_structure_file(args.model_dir, pdb_id)

        base_result = row.to_dict()
        base_result["dockq_csv_row_index"] = idx
        base_result["ref_file"] = ref_path
        base_result["model_file"] = model_path
        base_result["interface_error"] = ""

        if ref_path is None:
            base_result["interface_error"] = f"Reference file not found for {pdb_id}"
            results.append(base_result)
            continue

        if model_path is None:
            base_result["interface_error"] = f"Model file not found for {pdb_id}"
            results.append(base_result)
            continue

        try:
            if ref_path not in structure_cache:
                structure_cache[ref_path] = load_structure(ref_path)

            if model_path not in structure_cache:
                structure_cache[model_path] = load_structure(model_path)

            ref_structure = structure_cache[ref_path]
            model_structure = structure_cache[model_path]

            metrics = compare_one_interface(
                row,
                ref_structure,
                model_structure,
                args.cutoff,
            )

            base_result.update(metrics)

        except Exception as e:
            base_result["interface_error"] = repr(e)

        results.append(base_result)

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.out_csv, index=False)

    print(f"Wrote: {args.out_csv}")
    print(f"Rows written: {len(out_df)}")

    if "interface_error" in out_df.columns:
        n_errors = (out_df["interface_error"].fillna("") != "").sum()
        print(f"Rows with errors: {n_errors}")

        if n_errors > 0:
            print("\nMost common errors:")
            print(
                out_df.loc[
                    out_df["interface_error"].fillna("") != "",
                    "interface_error",
                ]
                .value_counts()
                .head(10)
                .to_string()
            )


if __name__ == "__main__":
    main()