#!/usr/bin/env python3
"""
Generate multiple 1D plots from data
"""
import numpy as np
import matplotlib.pyplot as plt

######## CONFIGURATION ########
# numpy.loadtxt() stuff
fname = 'test.txt'
skiprows = 0
usecols = None

plot_title = ''
plot_xlabel = ''
plot_ylabel = ''
# list of kwargs dictionaries for each plotted line
plot_options = []
######## END CONFIGURATION ########

data = np.loadtxt(fname, skiprows=skiprows, usecols=usecols)
# assume data points are stratified across the larger dimension
try:
    if data.shape[0] < data.shape[1]:
        data = data.T   # transpose so each column is a variable
    # assume first column is the independent variable
    indep_var = data[:, 0]
    dep_vars = data[:, 1:]
except IndexError:  # only one dimension
    # insert a dummy independent variable
    indep_var = np.arange(len(data))
    dep_vars = np.reshape(data, (data.size, 1)) # force to be a column matrix

fig, ax = plt.subplots()
for i in range(dep_vars.shape[1]):
    # Default to no options if all the options have been used up
    kwargs = plot_options[i] if i < len(plot_options) else {}
    # Recycle plot options if there aren't enough to cover all variables
    line, = ax.plot(indep_var, dep_vars[:, i], **kwargs)
    # Show all labels (strip surrounding underscores)
    line.set_label(line.get_label().strip('_'))
ax.legend()
ax.set_title(plot_title)
ax.set_xlabel(plot_xlabel)
ax.set_ylabel(plot_ylabel)

plt.show()
