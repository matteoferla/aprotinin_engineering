"""
This was originally a Jupyter notebook,
but I copy-pasted here for clarity.
"""

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
# ## Load

wt = pyrosetta.pose_from_file('3U1J_relaxed.pdb')
pose = wt.clone()

# =============================================================================
# ## Check chain 3 is the aprotinin

idx1s_chain_begins = [wt.chain_begin(c) for c in range(1, wt.num_chains()+1)]
pi = wt.pdb_info()
chain2idx1 = {pi.chain( idx ): idx for idx in idx1s_chain_begins}
print(chain2idx1)
#  {'A': 1, 'B': 21, 'E': 187}

# =============================================================================
# ## Get idx1 of the cysteines

cyss = []
for idx1 in range(pose.chain_begin(3), pose.total_residue() + 1):
    if pose.residue(idx1).name3() == 'CYS':
        cyss.append(idx1)
print(cyss)
# [191, 200, 216, 224, 237, 241]


# =============================================================================
# ## Mutate

for idx1 in cyss:
    prp.simple_moves.MutateResidue(idx1, 'ALA').apply(pose)

print([f'{o}{i+1}{m}' for i, (o, m) in enumerate( zip(wt.chain_sequence(3), pose.chain_sequence(3)) ) if o!=m])
# ['C5A', 'C14A', 'C30A', 'C38A', 'C51A', 'C55A']

# =============================================================================
# ## FastDesign

def create_design_tf(pose:pyrosetta.Pose, design_sele: pr_res.ResidueSelector, distance:int) -> prc.pack.task.TaskFactory:
    """
    Create design task factory for relax.
    Designs the ``design_sele`` and repacks around ``distance`` of it.

    Remember to do

    ... code-block:: python

        relax.set_enable_design(True)
        relax.set_task_factory(task_factory)
    """
    #residues_to_design = design_sele.apply(pose)
    # this is default:
    # design_ops = prc.pack.task.operation.OperateOnResidueSubset(????, residues_to_design)
    no_cys = pru.vector1_std_string(1)
    no_cys[1] = 'CYS'
    no_cys_ops =  prc.pack.task.operation.ProhibitSpecifiedBaseResidueTypes(no_cys)
    # No design, but repack
    repack_sele = pr_res.NeighborhoodResidueSelector(design_sele, distance, False)
    residues_to_repack = repack_sele.apply(pose)
    repack_rtl = prc.pack.task.operation.RestrictToRepackingRLT()
    repack_ops = prc.pack.task.operation.OperateOnResidueSubset(repack_rtl, residues_to_repack)
    # No repack, no design
    frozen_sele = pr_res.NotResidueSelector(pr_res.OrResidueSelector(design_sele, repack_sele))
    residues_to_freeze = frozen_sele.apply(pose)
    prevent_rtl = prc.pack.task.operation.PreventRepackingRLT()
    frozen_ops = prc.pack.task.operation.OperateOnResidueSubset(prevent_rtl, residues_to_freeze)
    # pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT
    # pyrosetta.rosetta.core.pack.task.operation.PreserveCBetaRLT
    task_factory = prc.pack.task.TaskFactory()
    task_factory.push_back(no_cys_ops)
    task_factory.push_back(repack_ops)
    task_factory.push_back(frozen_ops)
    return task_factory

# -----------------------------------------------------------------------------
# ### Only cysteines first

idx_sele = pr_res.ResidueIndexSelector()
for idx1 in cyss:
    idx_sele.append_index(idx1)

task_factory: prc.pack.task.TaskFactory = create_design_tf(pose, design_sele=idx_sele, distance=6)
scorefxn = pyrosetta.get_fa_scorefxn()
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
relax.set_enable_design(True)
relax.set_task_factory(task_factory)
relax.apply(pose)
pose.dump_pdb('decysteinified.pdb')
print(scorefxn(wt), scorefxn(pose) )
# -745.3 -764.0
print([f'{o}{i+1}{m}' for i, (o, m) in enumerate( zip(wt.chain_sequence(3), pose.chain_sequence(3)) ) if o!=m])
# ['C5V', 'C14A', 'C30A', 'C38V', 'C51T', 'C55A']

# -----------------------------------------------------------------------------
# ### Neighbours

