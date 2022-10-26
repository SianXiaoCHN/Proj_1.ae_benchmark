# Sian Xiao, SMU
# 2022 Mar 5

"""model training
Usage: python 10.ubq_training.py {LATENT_DIM}
Example: python 10.ubq_training.py 4
Output: 
"""

import sys

import numpy as np
from sklearn.preprocessing import MinMaxScaler

from models import vae
from utils.evaluation import Evaluation

# define detail
LATENT_DIM = int(sys.argv[1])

# define parameters
LAYERS = 4
EPOCHS = 200
BATCH_SIZE = 256
template_file = '../data/pdbs/ubq.pdb'

######################################## Data ########################################
# data loading
data = np.load("../data/features/ubq.npz")

coor_ori = data['coor']
input_dim = (coor_ori.shape[1],)

# scale data
scaler = MinMaxScaler()
feature_scale = scaler.fit_transform(coor_ori[::4])


######################################## Model ########################################
# hidden layers
hidden_layers = [LATENT_DIM]
while len(hidden_layers) <= LAYERS:
    hidden_layers.append(hidden_layers[-1] * 4 if hidden_layers[-1]
                         * 4 <= coor_ori.shape[1] else coor_ori.shape[1])

hidden_layers = hidden_layers[::-1]
hidden_layers.pop()
assert len(hidden_layers) == LAYERS, "Wrong layers"


# fit and evaluate VAE model
encoder, decoder, model = vae.build_vae(
    input_dim, hidden_layers, latent_dim=LATENT_DIM)
print(encoder.summary())

model.fit(
    x=feature_scale, y=feature_scale,
    shuffle=True,
    epochs=EPOCHS,
    batch_size=BATCH_SIZE)

eva = Evaluation(encoder, decoder, feature_scale, scaler)
spearmanr = eva.cal_spearman()
pearson = eva.cal_pearson()
rmsd = eva.cal_rmsd()
dope = eva.cal_dope(template_file, stride=100)

with open('ubq.txt', 'a') as f:
    f.write(f'{LATENT_DIM}: {spearmanr}, {pearson}, {rmsd}, {dope}\n')
