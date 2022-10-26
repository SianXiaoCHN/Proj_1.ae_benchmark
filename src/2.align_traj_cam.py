# Sian Xiao, SMU
# 2022 Mar 4

"""align cam simulation trajs
Warning: must be used in env MDAnalysis, which has numpy=1.21 but main env has numpy=1.19
Usage: python 2.align_traj_cam.py
Output: ../data/trajs_ref/cam{0-4}_refed.dcd
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align


# dcd list
psf_list = [f"../data/trajs_ori/cam{i}.psf" for i in range(5)]
dcd_list = [f"../data/trajs_ori/cam{i}.dcd" for i in range(5)]


traj_list = [mda.Universe(psf_list[i], dcd_list[i]) for i in range(5)]

for i in range(5):
    align.AlignTraj(
        traj_list[i], traj_list[0], select='protein and not name H*',
        filename=f"../data/trajs_ref/cam{i}_refed.dcd").run()
