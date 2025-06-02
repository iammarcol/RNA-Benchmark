import os
import subprocess
import pandas as pd
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Compute TM-scores between experimental and model PDB files using USalign with complex comparison (-mm 1 -d 5 -ter 0 for complexes).')
parser.add_argument('--afm_dir', required=True, help='Directory containing predicted PDB files (AFM/AF3)')
parser.add_argument('--exp_dir', required=True, help='Directory containing experimental PDB files')
parser.add_argument('--output_csv', required=True, help='Path to output CSV file')
args = parser.parse_args()

afms_dir = args.afm_dir
exps_dir = args.exp_dir
output_csv = args.output_csv

results = []

def run_command(command):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result

def extract_usalign_tm_score(output):
    last_tm_score = 'N/A'
    for line in output.splitlines():
        if line.startswith('TM-score='):
            last_tm_score = line.split('=')[1].strip().split()[0]
    return last_tm_score

# iterate, compute tm-score, extract the score
for exp_pdb_file in os.listdir(exps_dir):
    if exp_pdb_file.endswith('.pdb'):
        exp_pdb_path = os.path.join(exps_dir, exp_pdb_file)
        model_pdb_path = os.path.join(afms_dir, exp_pdb_file)
        exp_pdb_id = exp_pdb_file.split('.')[0]

        if os.path.exists(model_pdb_path):
            command = f"USalign {model_pdb_path} {exp_pdb_path} -mm 1 -d 5 -ter 0"
            result = run_command(command)
            tm_score = extract_usalign_tm_score(result.stdout)
            results.append([exp_pdb_id, f'm_{exp_pdb_id}', tm_score])

# save
pd.DataFrame(results, columns=['PDB_ID', 'Model_ID', 'TM_Score']).to_csv(output_csv, index=False)
print(f'TM scores have been saved to {output_csv}')
