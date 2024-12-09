"""
This was originally a Jupyter notebook,
but I copy-pasted here for clarity.
"""
from typing import Union, Sequence, Optional, List
import pandas as pd
from pathlib import Path

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

def score_interface(complex: Union[pyrosetta.Pose, Sequence[pyrosetta.Pose]], interface: str):
    if isinstance(complex, Sequence):
        _complex = complex[0].clone()
        for c in complex[1:]:
            add_chain(_complex, c)
        complex = _complex
    ia = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(interface)
    ia.apply(complex)
    return {'complex_energy': ia.get_complex_energy(),
            'separated_interface_energy': ia.get_separated_interface_energy(),
            'complexed_sasa': ia.get_complexed_sasa(),
            'crossterm_interface_energy': ia.get_crossterm_interface_energy(),
            'interface_dG': ia.get_interface_dG(),
            'interface_delta_sasa': ia.get_interface_delta_sasa()}

def assess_variant(variant, wt, interface: str='A_B', chain_idx=1,
                   scorefxn:Optional[pr_scoring.ScoreFunction]=None):
    if scorefxn is None:
        scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    info = {'dG_variant': scorefxn(variant), 'dG_wt': scorefxn(wt)}
    info['ddG'] = info['dG_variant'] - info['dG_wt']
    info['mutations'] =  ' '.join([f'{o}{i+1}{m}' for i, (o, m) in enumerate( zip(wt.chain_sequence(chain_idx), variant.chain_sequence(chain_idx)) ) if o!=m])
    info.update(score_interface(variant, interface))
    return info

# =============================================================================
# ## Load

wt = pyrosetta.pose_from_file('3U1J_relaxed.pdb')

# load prior variants
scores = []
for variant_path in [Path('3U1J_relaxed.pdb')] + list(Path('variants').glob('*.pdb')):
    variant = pyrosetta.pose_from_file(variant_path.as_posix())
    info = assess_variant(variant, wt, 'E_AB', 3)
    info['Name'] = variant_path.stem
    scores.append(info)

scores[0]['Name'] = 'wt'

for i in range(100):
    variant = wt.clone()
    #variants.append(variant)
    chainB_sele = pr_res.ChainSelector('B')
    chainE_sele = pr_res.ChainSelector('E')
    neigh_sele = pr_res.NeighborhoodResidueSelector(chainB_sele, 8, False)
    sele = pr_res.AndResidueSelector(neigh_sele, chainE_sele)
    task_factory: prc.pack.task.TaskFactory = create_design_tf(variant, design_sele=sele, distance=6)
    scorefxn = pyrosetta.get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
    relax.set_enable_design(True)
    relax.set_task_factory(task_factory)
    relax.apply(variant)
    info = assess_variant(variant, wt, 'E_AB', 3)
    info['Name'] = f'design_{i}'
    scores.append(info)
    variant.dump_pdb(f'variants/design_{i}.pdb')
    pd.DataFrame(scores).to_csv('all_scores.csv')

for i in range(500):
    variant = wt.clone()
    #variants.append(variant)
    chainB_sele = pr_res.ChainSelector('B')
    chainE_sele = pr_res.ChainSelector('E')
    neigh_sele = pr_res.NeighborhoodResidueSelector(chainB_sele, 8, False)
    sele = pr_res.AndResidueSelector(neigh_sele, chainE_sele)
    task_factory: prc.pack.task.TaskFactory = create_design_tf(variant, design_sele=sele, distance=6)
    scorefxn = pyrosetta.get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
    relax.set_enable_design(True)
    relax.set_task_factory(task_factory)
    relax.apply(variant)
    info = assess_variant(variant, wt, 'E_AB', 3)
    info['Name'] = f'minidesign_{i}'
    scores.append(info)
    variant.dump_pdb(f'variants/minidesign_{i}.pdb')
    pd.DataFrame(scores).to_csv('all_scores.csv')