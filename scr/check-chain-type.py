# execute 2nd

import os
import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser, is_aa, MMCIFParser

def determine_chain_type(chain):
    # 3-letter protein residue names
    protein_residues = {
        'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
        'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL',
        'TRP', 'TYR'
    }

    # RNA-like residues
    rna_u_like = {'U', 'PSU'} 

    # DNA residues
    dna_deoxy = {
        'DA', 'DC', 'DG', 'DT', 'DI',
        # methylated DNA variants
        '5MC', '5MDC', '5MDU'
    }

    # single-letter nucleotides (both RNA/DNA)
    single_nuc_ambiguous = {'A', 'C', 'G', 'I'}
    single_t = {'T'}  # single-letter T: DNA indicator if no U

    has_protein = False
    has_u_like = False
    has_deoxy = False
    has_t = False
    has_ambig_nuc = False

    for residue in chain:
        resname = residue.get_resname().strip()

        if resname in protein_residues:
            has_protein = True
            continue

        if resname in rna_u_like:
            has_u_like = True
            continue

        if resname in dna_deoxy:
            has_deoxy = True
            continue

        if resname in single_t:
            has_t = True
            continue

        if resname in single_nuc_ambiguous:
            has_ambig_nuc = True
            continue

    # classification logic 
    # if any protein + any nucleic → mixed
    if has_protein and (has_u_like or has_deoxy or has_t or has_ambig_nuc):
        return "protein_nucleic"
    if has_protein:
        return "protein"

    # no protein from here on → nucleic only
    # any U-like base → RNA
    if has_u_like:
        return "rna"
    # DNA markers (T or deoxy) and no U → DNA
    if has_deoxy or has_t:
        return "dna"
    # only A/C/G/I and no T/U/deoxy → default to RNA
    if has_ambig_nuc and not (has_deoxy or has_t or has_u_like):
        return "rna"
    # otherwise
    return "unknown"


######### PROCESS PDB FILES #########
def process_pdb_files(input_dir, output_csv):
    data = []
    parser = PDBParser(QUIET=True)

    for pdb_file in os.listdir(input_dir):
        if pdb_file.endswith(".pdb"):
            pdb_id = os.path.splitext(pdb_file)[0]
            file_path = os.path.join(input_dir, pdb_file)
            try:
                structure = parser.get_structure(pdb_id, file_path)
                for model in structure:
                    for chain in model:
                        chain_type = determine_chain_type(chain)
                        data.append({"PDB_ID": pdb_id, "Chain": chain.id, "Type": chain_type})
            except Exception as e:
                print(f"Error processing {pdb_file}: {e}")

    # save to CSV
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f"Output saved to {output_csv}")

# dirs
input_directory = "../complexes/method"  # input dir with model files for a specific method
output_csv_file = "../csv/chain-types-method.csv"  # output CSV file for chain types, for a specific method

process_pdb_files(input_directory, output_csv_file)
