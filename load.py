# load.py
import os, sys
sys.path.insert(0, './Models')
import numpy as np

M_e   = 5.972e24  #Earth mass (kg)
planetnum, mp_gp, rp_gp, period_gp = np.loadtxt('paper/5000_xp/5000_xp_gp_data.csv', usecols=(0,1,2),unpack=True)

# mass_mask = mp_gp < 7.6 and mp_gp > 7
# rp_mask = rp_gp > 2.5 and rp_gp < 2.7

# print(mp_gp[mass_mask])
# print(rp_gp[rp_mask])


mdot_all = np.loadtxt("Models/windae_grids/windae_nobolo_static_all.dat", skiprows=(29), usecols=8, unpack=True)

zero_mask = mdot_all == 0
# print(len(mdot_all[zero_mask]))
max_mask = mdot_all == max(mdot_all)
# print(len(mdot_all[max_mask]))

print(max(mdot_all[~max_mask]))
print(min(mdot_all[~zero_mask]))