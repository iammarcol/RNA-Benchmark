import os
import subprocess
import pandas as pd

# Define directories
afm_pdbs_dir = '/proj/berzelius-2021-29/users/x_malud/tmscore/compl_final/rf2na/trimmed'
exp_pdbs_dir = '/proj/berzelius-2021-29/users/x_malud/rna_benchmark/dataset/after/pdb/complex'
output_csv = '/proj/berzelius-2021-29/users/x_malud/tmscore/compl_final/out-new/usalign/tmscore-rf2na.csv'

# Initialize a list to store results
results = []

# Function to run a shell command
def run_command(command):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result

def extract_usalign_tm_score(output):
    last_tm_score = 'N/A'
    for line in output.splitlines(): 
        if line.startswith('TM-score='):
            last_tm_score = line.split('=')[1].strip().split()[0]  # Extract the TM-score value
    return last_tm_score


# Iterate over all experimental PDB files
for exp_pdb_file in os.listdir(exp_pdbs_dir):
    if exp_pdb_file.endswith('.pdb'):
        exp_pdb_path = os.path.join(exp_pdbs_dir, exp_pdb_file)
        exp_pdb_id = exp_pdb_file.split('.')[0]  # Extract PDB ID (e.g., '7dm0')

        # Find all corresponding model PDB files in afm_pdbs_dir
        matching_models = [f for f in os.listdir(afm_pdbs_dir) if f.startswith(exp_pdb_id) and f.endswith('.pdb')]

        for model_pdb_file in matching_models:
            model_pdb_path = os.path.join(afm_pdbs_dir, model_pdb_file)

            # Construct the USalign command
            command = f"USalign {model_pdb_path} {exp_pdb_path} -mm 1 -ter 1"

            # Run the command
            result = run_command(command)

            # Extract TM-score from the output
            tm_score = extract_usalign_tm_score(result.stdout)

            # Append result to the list
            results.append([exp_pdb_id, model_pdb_file, tm_score])

# Create a DataFrame and save to CSV
df = pd.DataFrame(results, columns=['PDB_ID', 'Model_ID', 'TM_Score'])
df.to_csv(output_csv, index=False)

print(f'TM scores have been saved to {output_csv}')
