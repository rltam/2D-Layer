import sys
import numpy as np
import netCDF4
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
import math
import scipy
from scipy.interpolate import BSpline, splrep
from Basis_Funcs import Spline_Basis
import shutil
import os
import datetime
from datetime import datetime
import astropy
from astropy.visualization import hist
from scipy.stats import norm

folder = '/home/rltam/Documents/RunOutputFiles/Residual_Analysis/'

# damping_range = [1e-8, 2.154e-8, 4.641e-8, 1e-7, 2.154e-7, 4.641e-7, 1e-6, 2.154e-6, 4.641e-6, 1e-5]
nbins = 20


meanResids = np.zeros(100)
meanResidsSquared = np.zeros(100)
# meanRS_2 = np.zeros(100)
# meanRS_3 = np.zeros(100)
# meanRS_4 = np.zeros(100)
# meanRS_5 = np.zeros(100)
# meanRS_6 = np.zeros(100)
# meanRS_7 = np.zeros(100)
# meanRS_8 = np.zeros(100)
# meanRS_9 = np.zeros(100)
# best = np.zeros((100,2))

# averes = np.loadtxt(f'{folder}/run0/AveResiduals.txt')

for i in range(100):
    resids = np.loadtxt(f'{folder}/{i}/Inversion/Residuals.txt(25x3)x1')
    meanRS_0[i] = resids[0][1]
    meanRS_1[i] = resids[1][1]
    meanRS_2[i] = resids[2][1]
    meanRS_3[i] = resids[3][1]
    meanRS_4[i] = resids[4][1]
    meanRS_5[i] = resids[5][1]
    meanRS_6[i] = resids[6][1]
    meanRS_7[i] = resids[7][1]
    meanRS_8[i] = resids[8][1]
    meanRS_9[i] = resids[9][1]
    best_val = min(resids[:,1])
    index = np.where(resids == best_val)
    pos = index[0][0]
    best[i][0] = resids[pos][0]
    best[i][1] = best_val

lst = [meanRS_0, meanRS_1, meanRS_2, meanRS_3, meanRS_4, meanRS_5, meanRS_6, meanRS_7, meanRS_8, meanRS_9]

# for i in range(10):
#     data = lst[i]
#     damp = damping_range[i]

#     mu, std = norm.fit(data)
#     plt.hist(data, bins = nbins)
#     xmin, xmax = plt.xlim()
#     x = np.linspace(xmin, xmax, 100)
#     p = norm.pdf(x, mu, std)
#     # plt.plot(x, p, c='k', label=f'Fitted Gaussian, mean = {mu}, std = {std}')
#     # plt.legend()
#     plt.xlabel('Mean Squared Residual')
#     plt.ylabel('Number of Occurences')
#     plt.title(f'Mean Squared Residuals Histogram for Damp = {damp}')
#     plt.savefig(f'{folder}/ResidualHist_{damp}.png')
#     plt.close()

# fig, ax = plt.subplots(2, 5, figsize = (25, 15))

# for i in range(10):
#     data = lst[i]
#     damp = damping_range[i]
#     mu, std = norm.fit(data)

#     if i < 5:
#         ax[0][i].hist(data, bins = nbins)
#         xmin, xmax = ax[0][i].set_xlim()
#         x = np.linspace(xmin, xmax, 100)
#         p = norm.pdf(x, mu, std)
#         # ax[0][i].plot(x, p, c='k')
#         ax[0][i].set_title(f'Damp = {damp}')
#         # ax[0][i].set_title(f'Damp = {damp}, mean = {mu}, std = {std}')
#     if i >= 5:
#         ax[1][i-5].hist(data, bins = nbins)
#         xmin, xmax = ax[1][i-5].set_xlim()
#         x = np.linspace(xmin, xmax, 100)
#         p = norm.pdf(x, mu, std)
#         # ax[1][i-5].plot(x, p, c='k')
#         ax[1][i-5].set_title(f'Damp = {damp}')
#         # ax[1][i-5].set_title(f'Damp = {damp}, mean = {mu}, std = {std}')

# fig.supxlabel('Mean Squared Residuals')
# fig.supylabel('Number of Occurences')
# fig.suptitle(f'Mean Squared Residual Histograms for all Damping')
# plt.savefig(f'{folder}/AllResidualHists')
# plt.close()


print(damping_range)
print(best[:,0])
damp_occurences = np.zeros(10)
for i in range(len(damping_range)):
    damp = damping_range[i]
    counter = 0
    for num in best[:,0]:
        if num == damp:
            counter += 1
    damp_occurences[i] = counter


fig, ax = plt.subplots(nrows=1, ncols=2, figsize = (18, 10))
# print(ax)
# print(ax.shape)
# print(ax[0])
# print(ax[1])

print(damp_occurences)

mu, std = norm.fit(best)
ax[0].hist(best[:,1], bins = nbins) 
xmin, xmax = ax[0].set_xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
# ax[0].plot(x, p, c='k', label=f'Fitted Gaussian, mean = {mu}, std = {std}')
# ax[0].legend()
ax[0].set_xlabel('Mean Squared Residual')
ax[0].set_ylabel('Number of Occurences')
ax[0].set_title(f'Best Fit Mean Squared Residuals Histogram')
ax[1].bar([str(damp) for damp in damping_range], damp_occurences)
ax[1].set_xlabel('Damping Factor')
ax[1].set_ylabel('Best Fit Occurences')
ax[1].set_title('Number of Runs Each Damping Factor Gives Best Fit')
ax[1].set_xticklabels([str(damp) for damp in damping_range], rotation = 60)
plt.savefig(f'{folder}/BestResidualsPlot.png')
plt.close()