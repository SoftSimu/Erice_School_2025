import MDAnalysis as mda
from g3 import MixedRadialAngularDistribution as G3
import matplotlib.pyplot as plt
import numpy as np

u = mda.Universe('data/spc-oplsaa-cubic-npt-run.gro', 'data/spc-oplsaa-cubic-npt-run.trr')
print(u)
print(len(u.trajectory))
print(set(u.atoms.names))

oxy = u.select_atoms('type O')
hyd = u.select_atoms('type H')
g3 = G3(oxy, oxy, hyd, 10)
g3.run()
dist = g3.distribution.T
vmax = np.percentile(dist[dist!=0], 90)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,3.5), constrained_layout=True)
ax.imshow(g3.distribution.T, extent=[0, 10, -1, 1], vmax=vmax, aspect='auto')
ax.set_xlabel('r /A')
ax.set_ylabel(r'cos$\Theta$')
fig.savefig('distribution.png')
