# Sian Xiao, SMU
# 2022 May 19

"""load aligned cam simulation trajs and calculate rmsd
Warning: must be used in env MDAnalysis, which has numpy=1.21 but main env has numpy=1.19
Usage: python 4.cam_rmsd.py
Output: 
"""

import MDAnalysis
import MDAnalysis.analysis.rms
import numpy as np


# psf, dcd, trajectories
psf_list = [f"../data/trajs_ori/cam{i}.psf" for i in range(5)]
dcd_list = [f"../data/trajs_ref/cam{i}_refed.dcd" for i in range(5)]
traj_list = [MDAnalysis.Universe(psf_list[i], dcd_list[i]) for i in range(5)]


# backbone RMSD
for i in range(5):
    u = traj_list[i]
    R = MDAnalysis.analysis.rms.RMSD(u, u,
           select="backbone")  # superimpose on whole backbone of the whole protein
    R.run()
    rmsd = R.rmsd.T  # transpose makes it easier for plotting
    time = rmsd[1] # frame, time, rmsd
    np.savez(f'cam{i}.npz', rmsd_=rmsd[2], time=time)


# load and save together
time = np.load('cam0.npz')['time']
cam0 = np.load('cam0.npz')['rmsd_']
cam1 = np.load('cam1.npz')['rmsd_']
cam2 = np.load('cam2.npz')['rmsd_']
cam3 = np.load('cam3.npz')['rmsd_']
cam4 = np.load('cam4.npz')['rmsd_']

# save rmsd
np.savez('cam_rmsd.npz', cam0=cam0, cam1=cam1, cam2=cam2, cam3=cam3, cam4=cam4, time=time)
