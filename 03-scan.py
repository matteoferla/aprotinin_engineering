"""
This was originally a Jupyter notebook,
but I copy-pasted here for clarity.
"""
from typing import List

import pyrosetta
import pyrosetta_help as ph
from types import ModuleType
from IPython.display import display, HTML
prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility
prcc: ModuleType = pyrosetta.rosetta.core.conformation
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_options: ModuleType = pyrosetta.rosetta.basic.options
pr_res: ModuleType = pyrosetta.rosetta.core.select.residue_selector

pyrosetta.init("-mute all")

# =============================================================================
# ## Define scan
# No point in scanning everything, only interface... While checking in PyMOL, I got these:

# print [atom.resi for atom in cmd.get_model('name CA and chain E and (byres chain B around 5)').atom]
# ['11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '32', '34', '36', '37', '38', '38', '45', '46']
# I could repeat that in PyRosetta, but why?

res_pidxs = ['11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '32', '34', '36', '37', '38', '38', '45', '46']

wt = pyrosetta.pose_from_file('3U1J_relaxed.pdb')
pi = wt.pdb_info()
res_idx1s =[pi.pdb2pose(res=int(pidx), chain='E') for pidx in res_pidxs]
aas = 'ADEFGHIKLMNPQRSTVWY' # no cys
mut_prefices: List[str] = [f'{pose.residue(pi.pdb2pose(res=int(pidx), chain="E")).name1()}{pidx}' for pidx in res_pidxs]
mutationgroups: List[List[str]] = [[f'{prefix}{aa}' for aa in aas] for prefix in mut_prefices]

# =============================================================================
# scan multithreaded

from pebble import ProcessPool
from concurrent.futures import as_completed
from typing import Dict
import os
import pandas as pd

def initialize_and_score_mutations(mutation_set: Dict) -> Dict:
    """
    Worker function to initialize PyRosetta, load the pose, and score mutations.
    mutation_set: Dictionary containing parameters for scoring mutations.
    """
    pyrosetta.init("-mute all")  # Initialize PyRosetta within the worker

    # Load the pose
    pose = pyrosetta.pose_from_file(mutation_set['pose_file'])

    # Initialize MutantScorer
    mscorer = ph.MutantScorer(pose, modelname=mutation_set['modelname'])

    # Score mutations
    scores = mscorer.score_mutations(
        mutations=mutation_set['mutations'],
        chains=mutation_set['chains'],
        interfaces=mutation_set['interfaces'],
        preminimize=mutation_set['preminimize'],
        distance=mutation_set['distance'],
        cycles=mutation_set['cycles']
    )
    return {"mutation_set": mutation_set, "scores": scores}

# Define the mutation sets for parallel processing
mutation_sets = [
    {
        "pose_file": "3U1J_relaxed.pdb",
        "modelname": "scan",
        "mutations": mutationgroup,
        "chains": "E",
        "interfaces": [("interface", "E_AB")],
        "preminimize": False,
        "distance": 8,
        "cycles": 1
    } for mutationgroup in mutationgroups
]

# Use Pebble's ProcessPool for parallel processing
results = []
with ProcessPool(max_workers=int(os.getenv("SLURM_CPUS_PER_TASK", 4))) as pool:  # Adjust number of workers as needed
    future_to_mutation = {
        pool.schedule(initialize_and_score_mutations, args=(mutation_set,)): mutation_set
        for mutation_set in mutation_sets
    }

    for future in as_completed(future_to_mutation):
        try:
            result = future.result()  # Get the result of the worker function
            results.append(result)
        except Exception as exc:
            print(f"Mutation set {future_to_mutation[future]} generated an exception: {exc}")

# =============================================================================
# ## Analyse

import pandas as pd

df = pd.DataFrame([r for result in results for r in result['scores']])
df.to_csv('scan_scores.csv')
