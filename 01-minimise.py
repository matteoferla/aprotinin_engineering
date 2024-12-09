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

# ## Minimise

def relax(original: pyrosetta.Pose,
          constraint_weight: float=5,
          cycles: int=15,
          relax_to_start_coords:bool=True) -> pyrosetta.Pose:
    """
    Copypasted from ``https://github.com/matteoferla/Fragment-hit-follow-up-chemistry/blob/main/fragment_elaboration_scripts/pyrosetta_min.py#L24`
    """
    pose: pyrosetta.Pose = original.clone()
    # configure constraints
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, constraint_weight)
    scorefxn.set_weight(pr_scoring.ScoreType.angle_constraint, constraint_weight)
    scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, constraint_weight)
    pyrosetta.rosetta.basic.options.set_boolean_option('relax:constrain_relax_to_start_coords', relax_to_start_coords)
    pyrosetta.rosetta.basic.options.set_boolean_option('relax:coord_constrain_sidechains', relax_to_start_coords)
    # set up the relax sampler
    pyrosetta.rosetta.protocols.relax.FastRelax.register_options()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    relax.constrain_relax_to_start_coords(relax_to_start_coords)
    relax.apply(pose)
    return pose

pose = pyrosetta.toolbox.pose_from_rcsb('3U1J')
relaxed = relax(pose,cycles=5)
relaxed.dump_pdb('3U1J_relaxed.pdb')