#!/usr/bin/env python3
"""
Generate a 1D histogram from data
"""

import numpy as np
import matplotlib.pyplot as plt

fname = 'test.txt'
data = np.loadtxt(fname, skiprows=0, usecols=0)
plt.hist(data, bins=100, density=True)
plt.show()
