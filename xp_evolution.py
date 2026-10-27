#  xp_evo.py
import sys, time, os 
sys.path.insert(0, './Models')
import argparse

# Replace all 'multiprocessing' with 'multiprocess' on macos 
import multiprocess
from multiprocess import Pool
import numpy as np
import xpml, plot_module

cpu_ct = multiprocess.cpu_count()

arg = argparse.ArgumentParser()
arg.add_argument('-c', '--cores',     default=int(cpu_ct/2), type=int, action='store', help='Number of PCU cores to run for multiprocessing')
arg.add_argument('-s', '--savepath',  default='data/',       type=str, action='store', help='Directory to create data folders in')
arg.add_argument('-n', '--n_planets', default=int(5000),     type=int, action='store', help='Number of planets to generate')
arg.add_argument('-m', '--mdot_type', default='H2O',         type=str, action='store', help='Mass-loss escape type')
arg.add_argument('-g', '--generate',  default=True,          type=bool,action='store', help='Generate new population or not')

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
    print(f"Directory '{savefolder}' created successfully in {os.getcwd()}.")
except FileExistsError:
    print(f"Directory '{savefolder}' already exists in {os.getcwd()}. WARNING: Data in this folder will be overwritten.")
except Exception as e:
    print(f"An error occurred: {e}")


# n_planets,planet_ID,mp,rp,teq,wmf,a,mstar,lstar,savepath='./',H2O=True
if generate is True:
    xpml.planet_generator(n_planets, savepath)

#Either way, load ts in
ids, mp, rp, teq, wmf, a, mstar, lstar = np.loadtxt(f'{savepath}/{n_planets}_xp/5000_xp_gp_data.csv', 
                                                    delimiter='\t', skiprows=1, usecols=[0,1,2,5,6,10,13,14], unpack=True)

params = [(n_planets, ids[i], mp[i], rp[i], teq[i], wmf[i], a[i], mstar[i], lstar[i], mdot_type, savepath) for i in range(n_planets)]


print('\n Simulating evolution...')
# n_planets,planet_ID,mp,rp,teq,wmf,a,mstar,lstar,savepath='./',H2O=True
with Pool(processes=cores) as pool:
    pool.map(xpml.gp_xpmlevo(*params))

t_f = time.time()
t = t_f - t_0
t_mins = t/60
print('Evolution duration:', t_mins, "minutes")
print(f'\n Evolution complete. Fish of celebration -> <*))))<<\nPlotting...')

#Plot_module is NOT UPDATED
plot_module.trails_one(n_planets,n_planets,identifier="",datapath=savepath+savefolder)
print("Plots finished")