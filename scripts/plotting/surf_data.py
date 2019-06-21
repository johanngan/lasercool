#!/usr/bin/env python3
"""
Generate a surface plot from data in tall ordered format.
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

######## CONFIGURATION ########
# numpy.loadtxt() stuff
fname = 'output/*/*.out'
skiprows = 0
usecols = None  # Can only use three columns at a time

cmap = cm.jet
######## END CONFIGURATION ########

data = np.loadtxt(fname, skiprows=skiprows, usecols=usecols)

# Convert from tall format to meshgrid format
y = np.unique(data[:, 1])
x = data[::len(y), 0]
X, Y = np.meshgrid(x, y)
Z = np.reshape(data[:, 2], X.shape, order='F')

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap=cmap)
fig.colorbar(surf)

plt.show()