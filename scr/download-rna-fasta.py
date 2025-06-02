import os
import requests
import pandas as pd
from time import sleep

output_dir = "../FASTA/raw-fasta" # change it in case you want it saved to another directory
os.makedirs(output_dir, exist_ok=True)

url = "https://search.rcsb.org/rcsbsearch/v2/query"

# rna-only
query = {
  "query": {
    "type": "terminal",
    "service": "text",
    "parameters": {
      "attribute": "entity_poly.rcsb_entity_polymer_type",
      "operator": "exact_match",
      "value": "RNA"
    }
  },
  "return_type": "entry",
  "request_options": {
    "paginate": {
      "start": 0,
      "rows": 10000 
    },
    "results_content_type": ["experimental"]
  }
}

response = requests.post(url, json=query)
if response.status_code != 200:
    raise RuntimeError(f"Failed API call: {response.status_code} {response.text}")

pdb_ids = [entry["identifier"] for entry in response.json().get("result_set", [])]
print(f"Found {len(pdb_ids)} RNA PDB entries.")

def download_fasta(pdb_id, save_dir):
    fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    response = requests.get(fasta_url)
    if response.status_code == 200 and ">" in response.text:
        with open(os.path.join(save_dir, f"{pdb_id}.fasta"), "w") as f:
            f.write(response.text)
        print(f"Downloaded {pdb_id}.fasta")
        return True
    else:
        print(f"Failed to download {pdb_id}.fasta")
        return False

successful, failed = [], []
for i, pdb_id in enumerate(pdb_ids):
    if download_fasta(pdb_id, output_dir):
        successful.append(pdb_id)
    else:
        failed.append(pdb_id)
    if i % 100 == 0:
        sleep(0.5)

# save csv files
pd.DataFrame({"PDB_ID": successful}).to_csv(os.path.join(output_dir, "rna-fastas.csv"), index=False)
pd.DataFrame({"Failed_ID": failed}).to_csv(os.path.join(output_dir, "rna-fastas-failed.csv"), index=False)

print(f"\nDownloaded {len(successful)} FASTA files.")
print(f"{len(failed)} failed downloads logged.")
