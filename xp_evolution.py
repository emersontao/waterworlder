#  xp_evo.py
import sys, time, os 
sys.path.insert(0, './Models')
import argparse

# Replace all 'multiprocessing' with 'multiprocess' on macos 
import multiprocess
from multiprocess import Pool
import numpy as np
import pandas as pd
import xpml, plot_module

cpu_ct = multiprocess.cpu_count()

arg = argparse.ArgumentParser()
arg.add_argument('--cores',     default=int(cpu_ct/2), type=int, action='store', help='Number of PCU cores to run for multiprocessing')
arg.add_argument('--savepath',  default='data/',       type=str, action='store', help='Directory to create data folders in')
arg.add_argument('--n_planets', default=int(5000),     type=int, action='store', help='Number of planets to generate')
arg.add_argument('--mdot_type', default='H2O',         type=str, action='store', help='Mass-loss escape type')
arg.add_argument('--generate',                         action='store_true',      help='Call to generate a new population')

args = arg.parse_args()
cores = args.cores
savepath = args.savepath
n_planets = args.n_planets
mdot_type = args.mdot_type
generate = args.generate
print(f'\nCores in use: {cores}/{cpu_ct}')

t_0 = time.time()

savefolder=f'{savepath}/{n_planets}_xp/{mdot_type}/'
try:
    os.mkdir(savefolder)
    print(f"\nDirectory '{savefolder}' created successfully in {os.getcwd()}.")
except FileExistsError:
    print(f"\nDirectory '{savefolder}' already exists in {os.getcwd()}. WARNING: Data in this folder will be overwritten.")
except Exception as e:
    print(f"\nAn error occurred: {e}")


# n_planets,planet_ID,mp,rp,teq,wmf,a,mstar,lstar,savepath='./',H2O=True
if generate:
    print("Generating population of", n_planets, "steam worlds ...")
    xpml.planet_generator(n_planets, savepath)

#Either way, load ts in
# all numbers are in SI units

init_df = pd.read_csv(f'{savepath}/{n_planets}_xp/{n_planets}_xp_gp_data.csv', delimiter='\t', usecols=['# Planet/5000', 'Mp (kg)', 'Rp (m)', 'Teq (K)', 'WMF', 'sep (m)', 'Mstar (kg)', 'Lstar (W)'])
print(init_df.keys())
# exit()

ids = init_df['# Planet/5000'].to_numpy()
mp = init_df['Mp (kg)'].to_numpy()
rp = init_df['Rp (m)'].to_numpy()
teq = init_df['Teq (K)'].to_numpy()
wmf = init_df['WMF'].to_numpy()
a = init_df['sep (m)'].to_numpy()
mstar = init_df['Mstar (kg)'].to_numpy()
lstar = init_df['Lstar (W)'].to_numpy()

print(type(ids), ids[0])
# ids, mp, rp, teq, wmf, a, mstar, lstar = np.loadtxt(f'{savepath}/{n_planets}_xp/{n_planets}_xp_gp_data.csv', 
#                                                     delimiter='\t', skiprows=1, usecols=[0,1,2,5,6,10,13,14], unpack=True)

params = [([n_planets, ids[i], mp[i], rp[i], teq[i], wmf[i], a[i], mstar[i], lstar[i], mdot_type, savepath]) for i in range(len(ids))]
# print(params[0])

def pool_wrapper(param):
    return xpml.gp_xpmlevo(*param)

# exit()
print('\n Simulating evolution...')
# n_planets,planet_ID,mp,rp,teq,wmf,a,mstar,lstar,savepath='./',H2O=True
with Pool(processes=cores) as pool:
    pool.map(pool_wrapper, params)

t_f = time.time()
t = t_f - t_0
t_mins = t/60
print('Evolution duration:', t_mins, "minutes")
print(f'\n Evolution complete. Fish of celebration -> <*))))<<\nPlotting...')

#Plot_module is NOT UPDATED
plot_module.trails_one(250, int(n_planets), datapath=savepath)
plot_module.kde(5000,'H2O', 'PlanetS', datapath=savepath, plotpath='/plots')
print("Plots finished")