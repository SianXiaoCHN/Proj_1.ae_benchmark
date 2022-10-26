# Sian Xiao, SMU
# 2022 Mar 4

"""extract xyz coordinates from toho traj
Usage: python 8.featurize_toho.py
Output: ../data/features/toho.npz containing coor
"""

import mdtraj as md
import numpy as np


# load trajectories, extract heavy atom indices, and save flattened xyz data
system = md.load(
    f"../data/trajs_ref/toho_refed.dcd",
    top=f"../data/trajs_ori/toho.psf")
system_selected_atoms = system.topology.select_atom_indices("heavy")


# coor shape (no. of frames , no. of protein heavy atoms, 3)
coor_system = system.atom_slice(system_selected_atoms).xyz
# coor shape (no. of frames , 3 * no. of protein heavy atoms)
coor_flatten_system = coor_system.reshape(coor_system.shape[0], -1)
print(coor_flatten_system.shape)


# save the original coordinate
np.savez(
    f"../data/features/toho.npz",
    coor=coor_flatten_system)
