#code by Artem, Nov. 10 2025
import numpy as np 
import pandas as pd
# import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
# import seaborn as sns

# Read the file
# file = "windae_nobolo_static_all.dat"
file = "Models/windae_grids/windae_nobolo_static_110325.dat"
df = pd.read_csv(file,sep="\t",header=27)

# Dropping Teq=400K because that dimension is not complete and messes up everything
# df = df[df["Teq"] > 500]

# Keys: 'Mp', 'Rp', 'grav', 'Teq', 'a', 'Fxuv', 'mstar', 'lstar', 'Mdot',
#       'MdotOW17', 'ramp_to', 'run_wind', 'rochelim'

# Simple error handling method: the errorcode is the sum of all errors
df["error"] = np.abs(df["ramp_to"])+np.abs(df["run_wind"])+np.abs(df["rochelim"])

# Easier for the interpolator: sort all input in scending order
list_labels = ['mstar', 'Teq', 'Fxuv', 'grav', 'Mp']
df = df.sort_values(by=['mstar', 'Teq', 'Fxuv', 'grav', 'Mp'])

# Extracting axes
axis_all = [np.unique(df[c]) for c in list_labels]

# Getting dimension sizes
dim_all = [len(axis) for axis in axis_all]
dim_product = np.prod(dim_all).item()

# Getting values
mdot_values = df['Mdot'].to_numpy().reshape(dim_all)
error_values = df['error'].to_numpy().reshape(dim_all)

# Transforming dataframe into numpy array
# mi = pd.MultiIndex.from_product(axis_all, names=list_labels)
# mdot_values = (df.set_index(list_labels)['Mdot']
#               .reindex(mi)
#               .to_numpy()
#               .reshape(dim_all))

# error_values = (df.set_index(list_labels)['error']
#               .reindex(mi)
#               .to_numpy()
#               .reshape(dim_all))

# Making the RGI for mdot and error
# fill_value=None means values can be extrapolated
# I'm choosing fill_value=20 in the error to indicate extrapolation
rgi_mdot = RegularGridInterpolator(points=axis_all, values=mdot_values,
                                  method="slinear",
                                  bounds_error=False, fill_value=None)

rgi_error = RegularGridInterpolator(points=axis_all, values=error_values,
                                  method="slinear",
                                  bounds_error=False, fill_value=20)


# Actually, let's interpolate in log-log space. This works better
# for quantities that are related by power-laws, rather than
# linear relationships. For that, we need to convert all values to log 
# and then convert back to normal. For convenience, let's wrap everything
# in a standard function

# Making the rgi
rgi_mdot_log = RegularGridInterpolator(points=[np.log10(axis) for axis in axis_all], values=np.log10(np.clip(mdot_values,a_min=1e-20,a_max=None)),
                                  method="slinear",
                                  bounds_error=False, fill_value=None)

rgi_error_log = RegularGridInterpolator(points=[np.log10(axis) for axis in axis_all], values=error_values,
                                  method="slinear",
                                  bounds_error=False, fill_value=20)

#'mstar', 'Teq', 'Fxuv', 'grav', 'Mp'
def mdot(ms,teq,fxuv,grav,mp):
    ms = (ms).flatten()[0]
    teq = (teq).flatten()[0]
    fxuv = (fxuv).flatten()[0]
    grav = (grav).flatten()[0]
    mp = (mp).flatten()[0]
    # print(mp)
    # print([ms,teq,fxuv,grav,mp], type([ms,teq,fxuv,grav,mp]))
    input_arr = np.log10(np.array([ms,teq,fxuv,grav,mp]))
    log_mdot = rgi_mdot_log(input_arr)
    return 10**log_mdot

def error(ms,teq,fxuv,grav,mp):
    ms = (ms).flatten()[0]
    teq = (teq).flatten()[0]
    fxuv = (fxuv).flatten()[0]
    grav = (grav).flatten()[0]
    mp = (mp).flatten()[0]
    # mp = float(mp)
    input_arr = np.log10([ms,teq,fxuv,grav,mp])
    error = rgi_error_log(input_arr)
    return error

# -----------------------------------------------------------
# Testing (old, use next)
# Recovering values from the grid
# >>> mdot(0.55,1500,2.36354e3,7.24,2);error(0.55,1500,2.36354e3,7.24,2)
# array([5.97219e+24])
# array([3.])
# >>> mdot(0.55,1500,2.36354e3,7.24,4);error(0.55,1500,2.36354e3,7.24,4)
# array([30208100.00000022])
# array([0.])

# Interpolating in loglog-space makes things smoother
# >>> mdot(0.55,1500,2.36354e3,7.24,3);error(0.55,1500,2.36354e3,7.24,3)
# array([4.55634115e+14])
# array([1.2451125])
# >>> rgi_mdot([0.55,1500,2.36354e3,7.24,3]);rgi_error([0.55,1500,2.36354e3,7.24,3])
# array([2.986095e+24])
# array([1.5])

# -----------------------------------------------------------
# Testing
# Recovering values from the grid
# 0.20000	0.52019	7.240000	400	0.01644	7.33145e-02	0.100	0.00115	2.73029e+03	2.96036e+07	0	0	0
# >>> mdot(0.1,400,7.33145e-02,7.24,0.2);error(0.1,400,7.33145e-02,7.24,0.2)
# array([2730.29])
# array([0.])

