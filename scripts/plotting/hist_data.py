#!/usr/bin/env python3
"""
Generate a 1D histogram from data
"""

import numpy as np
import matplotlib.pyplot as plt

fname = 'output/*/*.out'
data = np.loadtxt(fname)
plt.hist(data, bins=100, density=True)
plt.show()
