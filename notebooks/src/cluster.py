import MDAnalysis as mda
from g3 import MixedRadialAngularDistribution as G3
from g3 import generate_similarity_matrix
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

u = mda.Universe('data/spc-oplsaa-cubic-npt-run.gro', 'data/spc-oplsaa-cubic-npt-run.trr')
print(u)
print(len(u.trajectory))
print(set(u.atoms.names))

oxy = u.select_atoms('type O')
hyd = u.select_atoms('type H')
ds = []
for atom in oxy:
    g3 = G3(mda.AtomGroup([atom]), oxy, hyd, 10, (100, 50))
    g3.run()
    dist = g3.distribution.T
    ds.append(dist)
ds = np.array(ds)
dist = ds
vmax = np.percentile(dist[dist!=0], 90)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,3.5), constrained_layout=True)
ax.imshow(dist[0,:,:], extent=[0, 10, -1, 1], vmax=vmax, aspect='auto', interpolation='none')
ax.set_xlabel('r /A')
ax.set_ylabel(r'cos$\Theta$')
fig.savefig('distribution_single.png')
plt.close()

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(7,3.5), constrained_layout=True)

mat = generate_similarity_matrix(ds)
print(mat.min(), mat.mean())
axes[0].imshow(mat)

pca = PCA(2).fit_transform(mat)
axes[1].scatter(pca[:,0], pca[:,1], s=1)
axes[1].set_xlabel('PC-1')
axes[1].set_ylabel('PC-2')
fig.savefig('ssimmat.png')
