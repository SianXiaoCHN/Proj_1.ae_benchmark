# Sian Xiao, SMU
# 2022 Apr 22

"""align toho simulation trajs to the first frame and average the coordinates
Warning:    must be used in env MDAnalysis, which has numpy=1.21 but main env has numpy=1.19
Usage:      python 1.align_traj_toho_new.py
Output:     ../data/trajs_ref/toho_in_npt_refed.dcd
            ../data/trajs_ref/toho_in_nvt_refed.dcd
            ../data/trajs_ref/toho_out_npt_refed.dcd
            ../data/trajs_ref/toho_out_nvt_refed.dcd
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


for point in ["in", "out"]:
    for state in ["npt", "nvt"]:
        psf = f"../data/trajs_ori/{point}.psf"
        dcd = f"../data/trajs_ori/{point}_{state}.dcd"
        traj = mda.Universe(psf, dcd)
        align.AlignTraj(
            traj, ref, select='protein and not name H*',
            filename=f"../data/trajs_ref/{point}_{state}_refed.dcd").run()
