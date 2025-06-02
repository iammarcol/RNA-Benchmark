# EXAMPLE file of downloading PDB files from the database
###########################################################

import os
import requests
import pandas as pd

# where to save the files
output_dir = "./rna_training_pdb"
os.makedirs(output_dir, exist_ok=True)

# cutoff date
#cutoff_date = "2021-09-30"

# define URL
url = "https://search.rcsb.org/rcsbsearch/v2/query"

# query - look only for RNA
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
            "rows": 10000  # increased from 25 to 10000 to match original behavior
        },
        "results_content_type": [
            "experimental"
        ],
        "sort": [
            {
                "sort_by": "score",
                "direction": "desc"
            }
        ],
        "scoring_strategy": "combined"
    }
}

# search query
response = requests.post(url, json=query)
if response.status_code != 200:
    print(f"Failed to retrieve data from RCSB API. Status code: {response.status_code}")
    print(f"Response: {response.text}")
    exit()

# extract list of pdbs
pdb_ids = [entry["identifier"] for entry in response.json().get("result_set", [])]

# print the number of results
print(f"Number of RNA structures: {len(pdb_ids)}")

# f-ion to download pdb files
def download_pdb(pdb_id, save_dir):
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb1"
    response = requests.get(pdb_url)
    if response.status_code == 200:
        with open(os.path.join(save_dir, f"{pdb_id}.pdb1"), "w") as f:
            f.write(response.text)
        print(f"Downloaded {pdb_id}.pdb")
    else:
        print(f"Failed to download {pdb_id}")

# downloading pdb files
for pdb_id in pdb_ids:
    download_pdb(pdb_id, output_dir)

# saving
df = pd.DataFrame({"PDB_ID": pdb_ids})
df.to_csv(os.path.join(output_dir, "rna_structures.csv"), index=False)
print(f"Saved PDB IDs to rna_structures.csv")

print(f"Downloaded {len(pdb_ids)} RNA structures to {output_dir}.")
