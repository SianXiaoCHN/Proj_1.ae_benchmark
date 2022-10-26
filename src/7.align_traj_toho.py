# Sian Xiao, SMU
# 2022 Mar 7

"""align toho simulation trajs to the first frame and average the coordinates
Warning: must be used in env MDAnalysis, which has numpy=1.21 but main env has numpy=1.19
Usage: python 7.align_traj_toho.py
Output: ../data/trajs_ref/toho_refed.dcd
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align


# load trajectories
system = mda.Universe(
    "../data/trajs_ori/toho.psf",
    "../data/trajs_ori/toho.dcd")

average = align.AverageStructure(
    system, system, select='protein and not name H*',
    ref_frame=0).run()

ref = average.universe

align.AlignTraj(
    system, ref,
    select='protein and not name H*',
    filename="../data/trajs_ref/toho_refed.dcd").run()
