#! /usr/bin/env python3

from pandas import *
import matplotlib.pyplot as plt
from pandas.tools.plotting import scatter_matrix
import sys

f = sys.stdin
data = read_csv(f, delim_whitespace=True)

# plt.scatter(data['changes'], data['diff'])
# plt.scatter(data['changes_p'], data['diff'])
# plt.scatter(data['area'], data['diff'])
# plt.scatter(data['areaxchange'], data['diff'])

axes = scatter_matrix(data)
corr = data.corr().as_matrix()
for i, j in zip(*plt.np.triu_indices_from(axes, k=1)):
    axes[i, j].annotate("%.3f" %corr[i,j], (0.8, 0.8), xycoords='axes fraction', ha='center', va='center')

plt.show()
