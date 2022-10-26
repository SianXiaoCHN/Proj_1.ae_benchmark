# Sian Xiao, SMU
# 2022 Mar 4

"""align ubq simulation trajs to the first frame
Warning: must be used in env MDAnalysis, which has numpy=1.21 but main env has numpy=1.19
Usage: python 0.align_traj_ubq.py
Output: ../data/trajs_ref/ubq_refed.dcd
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align


# load trajectories
system = mda.Universe(
    "../data/trajs_ori/ubq.psf",
    "../data/trajs_ori/ubq.dcd")


# align to the first frame and save
# Atom selection language: https://userguide.mdanalysis.org/stable/selections.html
align.AlignTraj(
    system, system, select='protein and not name H*',
    filename="../data/trajs_ref/ubq_refed.dcd").run()
