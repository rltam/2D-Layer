import numpy as np
from numpy import math
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from scipy.interpolate import BSpline, splrep, spalde, CubicSpline

Edge_Height=np.loadtxt("Edge_Height.txt")
Edge_Length=np.loadtxt("Edge_Length.txt")

left = -Edge_Length/2
right = Edge_Length/2

def Spline_Basis(num_bases, degree):
    num_knots = int(num_bases) - int(degree) + 1
    inner_knots = np.linspace(left, right, num=num_knots)
    full_knots = np.concatenate([[left]*degree, inner_knots, [right]*degree])
    bases = [BSpline.basis_element(full_knots[i:i + degree + 2], extrapolate=True) for i in range(num_bases)]

    return bases

