# Evolution model for atmospheric escape of water from sub-Neptune steam atmospheres
#   Utilizes the McMaster Synthetic Planet Population
#   Optimized for multiprocessing
# 
#   VERSION:        SUMMER 2025
#   LAST UPDATED:   JUNE 30 2025, DESKTOP
#   ATMOSPHERE:     Owen & Wu 2017 [H/He]
#   INTERIOR:       Aguichine et al. 2025
#   UNITS:          EARTH
# //////////////////////////////////////// #

import sys, time
sys.path.insert(0, './Models')
import multiprocessing             # replace all 'multiprocessing' with 'multiprocess' on macos 
from multiprocessing import Pool

import numpy as np
import module_mass_loss as mml
import xpml #, swe
from astropy.table import Table

# Check CPU count
cpu_ct = multiprocessing.cpu_count()
cores = cpu_ct    # Manually set number of cores or use all
print(f'Cores in use: {cores}/{cpu_ct}')



# //////////////////////////////////////// #
# RETRIEVE PLANETS
#   Planets can either be generated or imported.
#   Comment the code to whichever you will not use.
# Output Units: MKS

# ////// LOAD IN planet population /////// #
#   Using: McMaster Planet Population Synthesis Model (MPPSM)
#   WARNING: import in Earth units, need to be converted to MKS
loadin_path = './evolution-data/mppsm/filtered_winds_varalpha_noatmosphere_planets.txt'
# mp_gp,rp_gp,teq_gp,wmf_gp,sep_gp,mstar_gp,lstar_gp = np.loadtxt(loadin_path,
#                                                                 skiprows=(1),usecols=(0,1,3,4,6,9,10), unpack=True)
mp_gp,rp_gp,teq_gp,wmf_gp,sep_gp,mstar_gp,lstar_gp = np.loadtxt(loadin_path,
                                                                skiprows=(1),usecols=(0,2,3,4,6,9,10), unpack=True) 
                                                                #USING ONLY WATER RADIUS...
# Display planet info (Earth units)
print('PLANET INFORMATION')
gp_table = Table()
# gp_table["Number"] = planetnum
gp_table["Mass (M🜨)"] = mp_gp / mml.M_e
gp_table["Radius (R🜨)"] = rp_gp / mml.R_e
# gp_table["Core Radius (R🜨)"] = core_r_gp / mml.R_e
gp_table["Teq (K)"] = teq_gp
gp_table["WMF"] = wmf_gp
gp_table["Sep (au)"] = sep_gp / mml.au
gp_table["Lstar (L☉)"] = lstar_gp / mml.L_sun
gp_table["Mstar (M☉)"] = mstar_gp / mml.M_sun
# gp_table["Scale Height (m)"] = H_gp

gp_table["Mass (M🜨)"].format = "{:.3f}"
gp_table["Radius (R🜨)"].format = "{:.3f}"
# gp_table["Core Radius (R🜨)"].format = "{:.3f}"
gp_table["Teq (K)"].format = "{:.3f}"
gp_table["WMF"].format = "{:.3f}"
gp_table["Sep (au)"].format = "{:.3f}"
gp_table["Lstar (L☉)"].format = "{:.3f}"
gp_table["Mstar (M☉)"].format = "{:.3f}"
# gp_table["Scale Height (m)"].format = "{:.3f}"
print(gp_table)

n_planets = len(mp_gp)
planetnum = np.arange(1,n_planets+1,1,dtype=int)


# ///// GENERATE a Planet Population ///// #
#   Random, uniformly distributed
n_planets = 1 #int(1e5)  #how many planets we want to generate? #10e6
planetnum, mp_gp, rp_gp, core_m_gp, teq_gp, wmf_gp, sep_gp, period_gp, H_gp, lambda_gp, startype_gp, mstar_gp, lstar_gp = xpml.planet_generator(n_planets)

# exit()
# earth = xpml.scaleheight(225,mml.R_e,mml.M_e)
# print(earth)


# //////////////////////////////////////// #
# SIMULATE PLANETS
t_0 = time.time()
print('\n Simulating evolution...')

# Define complete evolution code for generated population
def gp_xpmlevo(planetnum,mp,rp,teq,wmf,sep,mstar,lstar):

    # Run XPML
    time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr = xpml.evolution(mp,rp,teq,wmf,sep,mstar,lstar)

    # SAVE final evolution data for individual planet
    # MPPSM:
    planet_filename = f'planet_{planetnum}-{n_planets}_evolution_data.csv'
    path = f'evolution-data/mppsm/wr/{planet_filename}' #save MPPSM planets

    # planet_filename = f'planet_{planetnum}-{n_planets}_evolution_data.csv'
    # path = f'evolution-data/gp-evo/{n_planets}_xp/{planet_filename}'
    # # gp_evo_header = f'Planet {planetnum}/{n_planets} \ttime [Gyr]\tmass [Me]\twmf\tradius [Re]'
    # np.savetxt(path, np.c_[time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr], fmt='%.5e', delimiter='\t')#, header=gp_evo_header)
    
    # np.savetxt(path, \
    #         np.c_[time_arr/mml.Gyr, mp_arr/mml.M_e, rp_arr/mml.R_e, wmf_arr], fmt='%.5f', delimiter='\t')#, header=gp_evo_header)

    return time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr


# //////////////////////////////////////// #
# OPTIMIZE AND RUN ON MULTIPLE CORES

def xpmlevo_wrapper(param):   #wrapper of function
    return gp_xpmlevo(*param)

#input list of lists per planet
param_xpmlevo = [(planetnum[i], mp_gp[i], rp_gp[i], teq_gp[i], wmf_gp[i], sep_gp[i], mstar_gp[i], lstar_gp[i]) for i in range(n_planets)]

# RUN ON MULTIPLE CORES
# Default number of cores: all
with multiprocessing.Pool(processes=cpu_ct) as pool:
    xpmlevo_results = pool.map(xpmlevo_wrapper, param_xpmlevo)


# print(xpmlevo_results)

# Display runtime
t_f = time.time()
t = t_f - t_0
t_mins = t/60
print(f'\n Total Elapsed: {t:.2f} s  |  {t_mins:.2f} minutes')
print(f'\n Evolution complete. Fish of celebration -> <*))))<< \n')
