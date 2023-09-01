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
from scipy import interpolate
from scipy.interpolate import BSpline, splrep
from Basis_Funcs import Spline_Basis
import shutil
import os
import datetime
from datetime import datetime

def make_hist()