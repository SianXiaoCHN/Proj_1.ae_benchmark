# Sian Xiao, SMU
# 2021 Dec 20
# todo not finished

"""
model performance evaluation class
"""

import os
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from scipy.stats import spearmanr, pearsonr
from utils.PDE_process import Protein
from utils.DOPE_score import DOPE
import time

import torch as th
import mp_nerf

seq = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
atoms = {'R':11,
       'H':10,
       'K':9,
       'D':8,
       'E':9,
       'S':6,
       'T':7,
       'N':8,
       'Q':9,
       'C':6,
       'U':6,
       'G':4,
       'P':7,
       'A':5,
       'I':8,
       'L':8,
       'M':8,
       'F':11,
       'W':14,
       'Y':12,
       'V':7}

class Evaluation():
    def __init__(self, encoder, decoder, data_scale):
        '''
        data: testing data
        '''

        self.data_scale = data_scale
        self.encoder = encoder
        self.decoder = decoder
        self.spearman = None
        self.pearson = None
        self.dope = None
        # temporary file
        self.temp_file = f"./temp{time.time()}.pdb"
        self._process()

    def _scale_reverse(self, d):
        d[:,:,:3] *= 2 * np.pi
        d[:,:,:3] -= np.pi
        d[:,:,3:6] *= np.pi
        d[:,:,6:] *= 2 * np.pi
        d[:,:,6:] -= np.pi
        return d

    def _angToCoor(self, d):
        '''
        convert one frame to coor
        th.tensor(), shape (Frame, L, 12)
        '''
        coor = []
        for i in range(len(d)):
            ang_tensor = th.tensor(d[i])
            scaffolds_ori = mp_nerf.build_scaffolds_from_scn_angles(
                seq=seq,
                angles=ang_tensor,
                coords=None,
                device='auto'
            )
            coords, cloud_mask = mp_nerf.proteins.protein_fold(**scaffolds_ori)
            curCoor = []
            for i in range(len(coords)):
                atom = atoms[seq[i]]
                for a in range(atom):
                    curCoor.extend(coords[i][a].tolist())
            coor.append(curCoor)
        return coor


    def _process(self):
        '''
        data processing
        '''

        # inverse transform the scaled data back to original angles
        self.data = self._scale_reverse(self.data_scale.reshape(self.data_scale.shape[0], -1, 12))

        # encode structure
        self.data_latent = self.encoder.predict(self.data_scale)

        # VAE outputs
        if type(self.data_latent) == list:
            self.data_latent = self.data_latent[0]

        # decode structure
        decoded_structure = self.decoder(self.data_latent).numpy()
        decoded_structure = decoded_structure.reshape(
            self.data_latent.shape[0], -1)

        # scale back to angles
        self.decoded_structure = decoded_structure.reshape(decoded_structure.shape[0], -1, 12)
        self.decoded_structure = self._scale_reverse(self.decoded_structure)

        # convert to coords
        self.data = self._angToCoor(self.data)
        self.decoded_structure = self._angToCoor(self.decoded_structure)

        # calculate distances
        self.dist_ori = np.square(
            euclidean_distances(self.data, self.data)
            ).flatten()

        self.dist_encoded = np.square(
            euclidean_distances(self.data_latent, self.data_latent)
            ).flatten()

    def cal_spearman(self, recalculation=False):
        if self.spearman is None or recalculation:
            self.spearman = spearmanr(self.dist_ori, self.dist_encoded)

        return self.spearman

    def cal_pearson(self, recalculation=False):
        if self.pearson is None or recalculation:
            self.pearson = pearsonr(self.dist_ori, self.dist_encoded)

        return self.pearson

    def cal_rmsd(self):
        rmsd = np.sqrt(np.sum(np.square(
            self.decoded_structure * 10 - self.data * 10), axis=1
            ) / (self.data.shape[1] // 3))

        return np.mean(rmsd), np.std(rmsd)

    def cal_dope(self, template_file, recalculation=False, stride=1):
        if self.dope is None or recalculation:
            # protein template
            protein = Protein()
            protein.extract_template(template_file)

            # store all dope scores
            self.dope = []

            combined_data = np.vstack(
                (self.data[::stride], self.decoded_structure[::stride])
                )

            for frame in combined_data:
                protein.load_coor(frame * 10)
                protein.write_file(self.temp_file)
                self.dope.append(DOPE(self.temp_file))

            # remove this temporary file
            os.system("rm %s" % self.temp_file)

        # calculate differences
        dopes_diff = []
        size = len(self.dope) // 2

        for i in range(size):
            # decoded - real
            cur_diff = self.dope[i + size] - self.dope[i]
            dopes_diff.append([cur_diff, cur_diff / self.dope[i] * 100])

        mean_vals = np.mean(dopes_diff, axis=0)

        return mean_vals