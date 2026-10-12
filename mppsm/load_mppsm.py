# Load and filter MPPSM synthetic planets
#   MPPSM: McMaster Planet Population Synthesis Model
#   This filters for sub-Neptune steam candidates
#   Sun analog planets.
#
# https://ui.adsabs.harvard.edu/abs/2022MNRAS.515.2548A/abstract
# //////////////////////////////////////// #

import numpy as np # type: ignore
import sys
sys.path.insert(0, './Models')
import swe, xpml
from astropy.table import Table # type: ignore


# //////////////////////////////////////// #
# Load original MPPSM planets

# File from B. Skinner
path = './evolution-data/mppsm/'
filename = 'winds_varalpha_noatmosphere_planets.txt'
filepath_init = path+filename

# semi-major axis (au), original mass (M_e), core wmf, 
# original radius (R_e) (not Bennett's model), period, teq (K) from stellar irradiation
sep_mppsm_init,mp_mppsm_init,wmf_mppsm_init,rp_mppsm_init,period_mppsm_init,teq_mppsm_init = np.loadtxt(
    filepath_init,skiprows=(1),usecols=(2,3,8,9,12,14),unpack=True)


# Print results
print("INITIAL MPPSM PLANETS")
mppsm_table = Table()
mppsm_table["Mass (M🜨)"] = mp_mppsm_init
mppsm_table["Radius (R🜨)"] = rp_mppsm_init
mppsm_table["Teq (K)"] = teq_mppsm_init
mppsm_table["WMF"] = wmf_mppsm_init
mppsm_table["Sep (au)"] = sep_mppsm_init
# sn_table["Mstar (M☉)"] = 

mppsm_table["Mass (M🜨)"].format = "{:.3f}"
mppsm_table["Radius (R🜨)"].format = "{:.3f}"
mppsm_table["Teq (K)"].format = "{:.3f}"
mppsm_table["WMF"].format = "{:.3f}"
mppsm_table["Sep (au)"].format = "{:.3f}"
# sn_table["Mstar (M☉)"].format = "{:.3f}"
print(mppsm_table)


# //////////////////////////////////////// #
# FILTER FOR STEAM SUB-NEPTUNE CANDIDATES

# Create a mask to filter out unwanted planets
mask = (mp_mppsm_init>1.0) & (mp_mppsm_init<20.0) & (teq_mppsm_init>400) & (teq_mppsm_init<1200)

sep_mppsm = sep_mppsm_init[mask]
mp_mppsm = mp_mppsm_init[mask]
wmf_mppsm = wmf_mppsm_init[mask]
rp_mppsm = rp_mppsm_init[mask]
period_mppsm = period_mppsm_init[mask]
teq_mppsm = teq_mppsm_init[mask]

# # Worse version of the whole mask thing
# for i in range(len(mp_mppsm)):
#     # Mp [4,12] and Teq >= 750
#     if mp_mppsm[i] >= 1.0 and mp_mppsm[i] <= 20.0 and teq_mppsm[i] >= 400.0 \
#         and teq_mppsm[i] <= 1200 and rp_mppsm[i] < 5 and rp_mppsm[i] >= 1.8:
#         # Is there an easier way to do this?
#         mp_sn_mppsm = np.append(mp_mppsm, mp_mppsm[i])
#         rp_sn_mppsm = np.append(rp_mppsm, rp_mppsm[i])
#         wmf_sn_mppsm = np.append(wmf_mppsm, wmf_mppsm[i])
#         teq_sn_mppsm = np.append(teq_mppsm, teq_mppsm[i])
#         sep_sn_mppsm = np.append(sep_mppsm, sep_mppsm[i])
#         period_sn_mppsm = np.append(period_mppsm, period_mppsm[i])

# Print results
mppsm_table = Table()
mppsm_table["Mass (M🜨)"] = mp_mppsm
mppsm_table["Radius (R🜨)"] = rp_mppsm
mppsm_table["Teq (K)"] = teq_mppsm
mppsm_table["WMF"] = wmf_mppsm
mppsm_table["Sep (au)"] = sep_mppsm
# sn_table["Mstar (M☉)"] = 

mppsm_table["Mass (M🜨)"].format = "{:.3f}"
mppsm_table["Radius (R🜨)"].format = "{:.3f}"
mppsm_table["Teq (K)"].format = "{:.3f}"
mppsm_table["WMF"].format = "{:.3f}"
mppsm_table["Sep (au)"].format = "{:.3f}"
# sn_table["Mstar (M☉)"].format = "{:.3f}"
print("MPPSM PLANETS: PART 2 ELECTRIC BOOGALOO")
print(mppsm_table)


M_sun = 1.9885e30 #Sun mass (kg)
L_sun = 3.828e26 #Total luminosity of the Sun (Watts = Joules/sec)
mstar_mppsm = np.full(len(mp_mppsm),1)
lstar_mppsm = np.full(len(mp_mppsm),1)
# mstar_mppsm = np.full(len(mp_mppsm),M_sun)
# lstar_mppsm = np.full(len(mp_mppsm),L_sun)


# //////////////////////////////////////// #
# Recalculate planet radii based on SWE model (A25)
import swe  # type: ignore
rp_swe = np.array([])
wmf_swe = np.array([])

for i in range(len(rp_mppsm)):
    age = 1e-2 * 3.15e16
    M_e = 5.97237e24 #Earth mass (kg)
    R_e = 6.371e6 #Earth radius (m)
    mp = mp_mppsm[i] * M_e
    rp = rp_mppsm[i] * R_e
    teq = teq_mppsm[i]
    wmf = wmf_mppsm[i]

    # calculate SWE WMF/Rp
    wmf_new = swe.solve_wmf_g(rp,teq,mp,age)
    rp_new = swe.radius_g(wmf,teq,mp,age)

    wmf_swe = np.append(wmf_swe,wmf_new)
    rp_swe = np.append(rp_swe,rp_new/R_e)

# add period and scaleheight as values
period_arr = xpml.period(mstar_mppsm,mp_mppsm*M_e,sep_mppsm)
print(period_arr)
scaleheight_arr = xpml.scaleheight(teq_mppsm,rp_mppsm,mp_mppsm)


# New table of filtered planets and new params
print("FILTERED SUB-NEPTUNES")
sn_table = Table()
sn_table["Mass (M🜨)"] = mp_mppsm
sn_table["Radius (R🜨)"] = rp_mppsm
sn_table["swe Radius (R🜨)"] = rp_swe
sn_table["Teq (K)"] = teq_mppsm
sn_table["WMF"] = wmf_mppsm
sn_table["SWE WMF"] = wmf_swe
sn_table["Sep (au)"] = sep_mppsm
sn_table["Period (days)"] = period_arr
sn_table["Scale Height"] = scaleheight_arr
# sn_table["Mstar (M☉)"] = 

sn_table["Mass (M🜨)"].format = "{:.3f}"
sn_table["Radius (R🜨)"].format = "{:.3f}"
sn_table["swe Radius (R🜨)"].format = "{:.3f}"
sn_table["Teq (K)"].format = "{:.3f}"
sn_table["WMF"].format = "{:.3f}"
sn_table["SWE WMF"].format = "{:.3f}"
sn_table["Sep (au)"].format = "{:.3f}"
sn_table["Period (days)"].format = "{:.3f}"
sn_table["Scale Height"].format = "{:.3f}"
# sn_table["Mstar (M☉)"].format = "{:.3f}"
print(sn_table)


# //////////////////////////////////////// #
# SAVE FILTERED SUB-NEPTUNES
filter_filename = 'filtered_winds_varalpha_noatmosphere_planets.txt'
filepath_filter = path+filter_filename
header = 'Filtered MPPSM Sub-Neptunes # Mp, Rp (M), Rp (SWE), Teq, WMF (M), WMF (SWE), Sep, Period, Scale height, Msun, Lsun'
np.savetxt(filepath_filter,np.c_[mp_mppsm,rp_mppsm,rp_swe,teq_mppsm,wmf_mppsm,
           wmf_swe,sep_mppsm,period_arr,scaleheight_arr,mstar_mppsm,lstar_mppsm], 
           header=header, fmt='%.5f', delimiter='\t')
print('File saved:',filepath_filter)
