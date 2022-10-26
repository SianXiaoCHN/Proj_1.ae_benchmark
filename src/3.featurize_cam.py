# Sian Xiao, SMU
# 2022 May 19

"""extract xyz coordinates from cam traj
Usage: python 3.featurize_cam.py
Output: ../data/features/cam.npz containing coor
"""

import mdtraj as md
import numpy as np

# load trajectories, extract heavy atom indices, and save flattened xyz data

# first system cam0
system = md.load(
    f"../data/trajs_ref/cam0_refed.dcd",
    top=f"../data/trajs_ori/cam0.psf")
system_selected_atoms = system.topology.select_atom_indices("heavy")

# coor shape (no. of frames , no. of protein heavy atoms, 3)
coor_system = system.atom_slice(system_selected_atoms).xyz
# coor shape (no. of frames , 3 * no. of protein heavy atoms)
coor_flatten_system = coor_system.reshape(coor_system.shape[0], -1)
print(coor_flatten_system.shape)

coor = coor_flatten_system


for i in range(1, 5):
    system = md.load(
        f"../data/trajs_ref/cam{i}_refed.dcd",
        top=f"../data/trajs_ori/cam{i}.psf")
    system_selected_atoms = system.topology.select_atom_indices("heavy")
    coor_system = system.atom_slice(system_selected_atoms).xyz
    coor_flatten_system = coor_system.reshape(coor_system.shape[0], -1)
    print(coor_flatten_system.shape)
    coor = np.concatenate(
        (coor, coor_flatten_system[:, :coor.shape[1]]), axis=0)

print(coor.shape)

# save the original coordinate
np.savez(
    f"../data/features/cam.npz",
    coor=coor)
