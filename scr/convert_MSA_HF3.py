import os
from Bio import AlignIO

def convert_afa_to_sto(input_afa, output_sto):
    """
    Convert a .afa (FASTA alignment) file to .sto (Stockholm) format, replacing T with U.

    Args:
        input_afa (str): Path to input .afa file.
        output_sto (str): Path to output .sto file.
    """
    try:
        alignment = AlignIO.read(input_afa, "fasta")
        os.makedirs(os.path.dirname(output_sto), exist_ok=True)
        with open(output_sto, "w") as f:
            AlignIO.write(alignment, f, "stockholm")
        with open(output_sto, "r") as f:
            updated = f.read().replace("T", "U")
        with open(output_sto, "w") as f:
            f.write(updated)
        print(f"[✓] Converted: {input_afa} → {output_sto}")
    except Exception as e:
        print(f"[!] Error with {input_afa}: {e}")


# This creates an output directory specifically for running HelixFold3 and puts newly generated .sto files inside

def batch_convert_afa_to_sto(input_dir, output_dir):
    """
    Convert all .afa files in a directory to .sto format with correct nested structure.

    Args:
        input_dir (str): Directory containing .afa files.
        output_dir (str): Base directory for output .sto files.
    """
    for fname in os.listdir(input_dir):
        if fname.endswith(".afa"):
            id_ = os.path.splitext(fname)[0]
            input_path = os.path.join(input_dir, fname)
            output_path = os.path.join(output_dir, id_, "msas", "rna_A", "A", "rfam_hits.sto")
            convert_afa_to_sto(input_path, output_path)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert .afa alignment files to .sto format.")
    parser.add_argument("input_dir", help="Path to directory containing .afa files")
    parser.add_argument("output_dir", help="Path to directory for output .sto files")

    args = parser.parse_args()
    batch_convert_afa_to_sto(args.input_dir, args.output_dir)
