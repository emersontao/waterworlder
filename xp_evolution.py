# Evolution model for atmospheric escape of water from sub-Neptune steam atmospheres
#   Generates randomized planet populations
#   Optimized for multiprocessing
# 
#   VERSION:        SUMMER 2026
#   LAST UPDATED:   JUNE 17 2026, KAMA (now runs both atm)
#   ATMOSPHERE:     Broome et al. 2025 [H2O] & Owen & Wu 2017 [H-He]
#   INTERIOR:       Aguichine et al. 2025
#   UNITS:          MKS
# //////////////////////////////////////// #

import sys, time, os 
sys.path.insert(0, './Models')

# Replace all 'multiprocessing' with 'multiprocess' on macos 
import multiprocess
from multiprocess import Pool

import numpy as np
import xpml #type:ignore
import xpml_HHe #type: ignore
import plot_module #type: ignore

#Check CPU count
cpu_ct = multiprocess.cpu_count()
cores = cpu_ct/2    #Manually set number of cores or use all
print(f'\nCores in use: {cores}/{cpu_ct}')

# //////////////////////////////////////// #
# RETRIEVE PLANETS
#   Planets can either be generated or imported.
#   Comment the code to whichever you will not use.
#   Output Units: MKS

#-- IMPORT planet population
#from XPML generated population (gp):
# n_planets = 200
# loadin_path = './evolution-data/paper/5000_xp/5000_xp_gp_data.csv'
# planetnum, mp_gp, rp_gp, mcore_gp, rsurf_gp, teq_gp, wmf_gp, grav_gp, H_gp, lambda_gp, sep_gp, period_gp, startype_gp, mstar_gp, lstar_gp, fxuv_gp = np.loadtxt(loadin_path,skiprows=(1), unpack=True)

# #from NEA: not gp but wtv
# filename = 'PS_2025.07.16_21.12.59.csv'
# planetnum, mp_gp, rp_gp, mcore_gp, rsurf_gp, teq_gp, wmf_gp, grav_gp, H_gp, lambda_gp, sep_gp, period_gp, startype_gp, mstar_gp, lstar_gp = xpml.nea_reformatter(filename)
# n_planets = len(planetnum)
# # print('Hooray!')
# # print(n_planets)
# # exit()

#-- GENERATE planet population
n_planets = 10
savepath="paper/"
mdot_type="H2O"
planetnum, mp_gp, rp_gp, mcore_gp, rsurf_gp, teq_gp, wmf_gp, grav_gp, H_gp, lambda_gp, sep_gp, period_gp, startype_gp, mstar_gp, lstar_gp, fxuv_gp = xpml.planet_generator(n_planets, savepath)

savefolder=f'{n_planets}_xp/'+mdot_type+'/'
try:
    os.mkdir(savepath+savefolder)
    print(f"Directory '{savefolder}' created successfully in {savepath}.")
except FileExistsError:
    print(f"Directory '{savefolder}' already exists in {savepath}. WARNING: Data in this folder may be overwritten.")
except Exception as e:
    print(f"An error occurred: {e}")

# exit()
# //////////////////////////////////////// #
# SIMULATE PLANETS
t_0 = time.time()
print('\n Simulating evolution...')

#Define complete evolution code for generated population
def gp_xpmlevo(planetnum,mp,rp,teq,wmf,sep,mstar,lstar):

    #Run XPML
    # time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr = xpml_HHe.evolution(mp,rp,teq,wmf,sep,mstar,lstar)

    # Save planet evolution data
    # -- Generated Planets
    planet_filename = f'planet_{int(planetnum)}-{int(n_planets)}_evolution_data.csv'
    # path = f'evolution-data/paper/{n_planets}_xp/H-He/{planet_filename}'
    # np.savetxt(path, np.c_[time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr], fmt='%.5e', delimiter='\t')


    #-- run for H2O
    time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr, error_arr = xpml.evolution(mp,rp,teq,wmf,sep,mstar,lstar, mdot_type)

    path = savepath+f'{n_planets}_xp/'+mdot_type+'/'+planet_filename


    # #-- NEA Planets
    # '''planet_filename = f'planet_{planetnum}-{n_planets}_evolution_data.csv'
    # path = f'evolution-data/real-evo/nea/{n_planets}_xp_{filename}/{planet_filename}' '''

    np.savetxt(path, np.c_[time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr, error_arr], fmt='%.5e', delimiter='\t') #,error_arr
    #error_arr
    return time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr, error_arr



# //////////////////////////////////////// #
# OPTIMIZE AND RUN ON MULTIPLE CORES
def xpmlevo_wrapper(param):
    return gp_xpmlevo(*param)

#input list of lists per planet
param_xpmlevo = [(planetnum[i], mp_gp[i], rp_gp[i], teq_gp[i], wmf_gp[i], sep_gp[i], mstar_gp[i], lstar_gp[i]) for i in range(n_planets)]

#-- RUN ON MULTIPLE CORES
# Default number of cores: all
with multiprocess.Pool(processes=cpu_ct) as pool:
    xpmlevo_results = pool.map(xpmlevo_wrapper, param_xpmlevo)

# Display runtime
t_f = time.time()
t = t_f - t_0
t_mins = t/60
print(f'\n Evolution complete. Fish of celebration -> <*))))<<\nPlotting...')


#now let's make some PLOTS, baby!
# plot_module.run(n_planets, mdot_type="H2O", datapath=savepath, plotkde=False, planets_plotted=10)

plot_module.trails_one(n_planets,n_planets,identifier="",datapath=savepath+savefolder)

print(f"\n Plotting complete. oh my god look it's a snail family!!! -> Y_@_ Y@_ v@_\n")
print(f'\n Total Elapsed: {t:.2f} s  |  {t_mins:.2f} minutes\n')
