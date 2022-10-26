import matplotlib.pyplot as plt
import numpy as np

fileName = '13.ubq.txt'

with open(fileName, 'r') as f:
    lines = f.readlines()

dims = [1,2,3,4,5,6,7,8, 9, 10]

spearmanrs = [[] for _ in range(len(dims))]
pearsonrs = [[] for _ in range(len(dims))]
rmsds = [[] for _ in range(len(dims))]
dopes = [[] for _ in range(len(dims))]


for line in lines:
    try:
        idx = dims.index(int(line.split(':')[0]))
    except:
        continue

    if 'nan' in line:
        continue

    spearmanr = float(line.split('=')[1].split(',')[0])
    spearmanrs[idx].append(spearmanr)

    pearsonr = float(line.split('(')[2].split(',')[0])
    pearsonrs[idx].append(pearsonr)

    rmsd = float(line.split(',')[-3].split('(')[-1])
    rmsds[idx].append(rmsd)

    dope = abs(float(line.split(']')[0].split(' ')[-1])) # change end 0
    dopes[idx].append(dope)


# plt.boxplot(spearmanrs, showmeans=True, showfliers=False)
# plt.xticks(range(1, len(spearmanrs)+1), dims)
# plt.ylabel('Spearman R')
# plt.xlabel('Latent dimension')
# plt.title('Spearman R for different dimensions')
# plt.savefig('11.ubq_spearmanr.png')
# plt.cla()

# plt.boxplot(pearsonrs, showmeans=True, showfliers=False)
# plt.xticks(range(1, len(pearsonrs)+1), dims)
# plt.ylabel('Pearsonr R')
# plt.xlabel('Latent dimension')
# plt.title('Pearson R for different dimensions')
# plt.savefig('11.ubq_pearsonr.png')
# plt.cla()

# plt.boxplot(rmsds, showmeans=True, showfliers=False)
# plt.xticks(range(1, len(spearmanrs)+1), dims)
# plt.ylabel('RMSD ($\AA$)')
# plt.xlabel('Latent dimension')
# plt.title('RMSD for different dimensions')
# plt.savefig('11.ubq_rmsd.png')

# plt.cla()
# plt.boxplot(dopes, showmeans=True, showfliers=False)
# plt.xticks(range(1, len(dopes)+1), dims)
# plt.ylabel('DOPE percent error')
# plt.xlabel('Latent dimension')
# plt.title('DOPE for different dimensions')
# plt.savefig('11.ubq_dope.png')


spearmanrs = np.array(spearmanrs, dtype=object)
pearsonrs = np.array(pearsonrs, dtype=object)
rmsds = np.array(rmsds, dtype=object)
dopes = np.array(dopes, dtype=object)

np.savez(
    'ubq_analyze.npz', spearmanrs=spearmanrs, pearsonrs=pearsonrs, rmsds=rmsds, dopes=dopes)
