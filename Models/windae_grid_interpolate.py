#code by Artem, Nov. 10 2025
#comments pruned by Emerson, Oct. 20 2026
import numpy as np 
import pandas as pd
# import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
# import seaborn as sns

# Read the file
# file = "windae_nobolo_static_all.dat"
file = "Models/windae_grids/windae_nobolo_static_110325.dat"
df = pd.read_csv(file,sep="\t",header=27)

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