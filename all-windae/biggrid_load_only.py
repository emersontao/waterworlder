# ------------------ biggrid.py ------------------- #
# LAST UPDATED: NOVEMBER 26, 2025

from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
from wind_ae.wrapper.wrapper_utils import constants as const
from wind_ae.wrapper.wrapper_utils.system import system

import numpy as np
import os, sys, time, signal
from itertools import product
import initialize_weights as init

#Constants
sigma_sb = 5.67037e-5 #erg⋅cm−2⋅s−1⋅K−4
G = 6.67e-8         #whatever that is in cgs
Gyr = 3.153e16      #s

#set up timeout
def signal_handler(signum, frame):
    raise Exception("Timed out!")


# -------------- INITIALIZE THE GRID -------------- #
# mp_axis     = np.array([0.2, 2, 4, 7, 11, 15, 18, 20, 22.5, 25, 27, 30]) #earth mass
mp_axis     = np.array([2, 4, 7, 11, 15, 18, 20, 22.5, 25, 27, 30]) #earth mass
grav_axis   = np.array([515, 600, 724, 1000, 1400, 2000, 3000]) #cm/s^2
teq_axis    = np.array([1500, 1225, 950, 675, 400]) #K
mstar_axis  = np.array([0.1, 0.22, 0.3, 0.55, 1]) #sun mass
Fxuv_axis   = np.array(list(reversed(init.log_Fxuv_axis(mstar_axis,teq_axis,5)))) #erg/s/cm^2

#add grid flipper
num_arg = len(sys.argv) #number of arguments; default is 1 (name of script)
# print(num_arg)
if num_arg > 1:
    if float(sys.argv[1]) == 1:
        teq_axis = np.array(list(reversed(teq_axis)))
        print("Reversing Teq: 400 -> 1500K")
    elif float(sys.argv[1]) == 2:
        mstar_axis = np.array(list(reversed(mstar_axis)))
        print("Reversing Mstar: 1 -> 0.1 Msun")
    elif float(sys.argv[1]) == 3:
        Fxuv_axis = np.array(list(reversed(Fxuv_axis)))
        print("Reversing Fxuv: low -> high")

ms_ext = np.tile(mstar_axis,len(teq_axis)) #extended ms
teq_ext = np.repeat(teq_axis,len(mstar_axis))
# print(teq_ext)
mask_ms_low = ms_ext < 0.43
lstar_ext = np.zeros_like(ms_ext)
lstar_ext[mask_ms_low] = 0.23*(ms_ext[mask_ms_low])**2.3
lstar_ext[~mask_ms_low] = ms_ext[~mask_ms_low]**4
sep_axis = np.sqrt(lstar_ext*const.Lsun / (16*np.pi* sigma_sb * teq_ext**4)) #cm

# Fxuv_axis(list(reversed(Fxuv_axis)))
# print(Fxuv_axis)
# print(type(Fxuv_axis))
# exit()

dim_mp = len(mp_axis)
dim_g = len(grav_axis)
dim_teq = len(teq_axis)
dim_mstar = len(mstar_axis)
dim_Fxuv = len(Fxuv_axis)
dim_sep = dim_teq * dim_mstar

#check for duplicate sep in our parameterspace
_, sep_cts = np.unique(sep_axis, return_counts=True)
if np.any(sep_cts > 1):
    print("WARNING: duplicate semi-major axis values detected in the grid.\nModify your Teq or Mstar to avoid this before running again.")
    exit()

# print(np.sort(sep_axis))
# exit()
#initialize the grid in the shape of given dimensions
input_grid = [mp_axis, grav_axis, teq_axis, mstar_axis, Fxuv_axis]
massloss_grid = np.zeros((dim_mp,dim_g,dim_sep,dim_Fxuv))
ramp_grid = np.zeros((dim_mp,dim_g,dim_sep,dim_Fxuv))
wind_grid  = np.zeros((dim_mp,dim_g,dim_sep,dim_Fxuv))

# with open("grids/rochelim_invalid_list.dat", "w") as file:
#     file.write("#invalid combinations for grid due to Roche limit as of 09/29/25\n")
G = 6.67e-11
Gyr = 3.1536e16 #seconds in a gigayear
#-- Owen & Wu 2017
def efficiency(mp,rp,mstar,sep):
    #efficiency of photoevaporative atmospheric H/He mass loss
    #Owen & Wu 2017
    v_esc = np.sqrt(2*G*mp / rp) #escape velocity at planet surface
    eff = 0.1 * (v_esc / 15e3)**(-2) #mass-loss efficiency, uncorrected
    
    #Roche lobe correction, Salz+ 2016a
    xi = (sep / rp) * (mp / (3*mstar))**(1/3)   #ξ
    k_diff = 1 - 3/(2 * xi) + 1/(2 * xi**3)     #fractional gravitational potential energy difference

    efficiency = eff/k_diff #corrected efficiency
    return efficiency

def L_xuv(age,lstar):
    #high-energy luminosity (X + UV) emitted by the star
    #Sanz-Forçada+2011

    #L_xuv in non-saturation regime
    Lx = 1.89e21 * (age/Gyr)**(-1.55)
    Leuv = 10**(3.82) * Lx**0.86
    Lxuv = Lx + Leuv

    #L_xuv in saturation regime
    Lxsat = lstar * 6.3e-4
    Leuvsat = 10**(3.82) * Lxsat**0.86
    Lxuvsat = Lxsat + Leuvsat

    new_lxuv = np.minimum(Lxuv,Lxuvsat) #final Lxuv
    return new_lxuv

def solve_mdot(mp,rp,mstar,sep,age,lstar):
    L_he = L_xuv(age,lstar)
    eta = efficiency(mp,rp,mstar,sep)
    mdot = eta * (L_he*rp**3) / (4*G*mp*sep**2)
    return mdot




#file details
date = "100226"
gridname = f"grids/windae_nobolo_static_{date}.dat"

header = """##################################################################################
# Grid of water world sub-Neptune atmospheric escape models 
#   wind-ae (Broome+2025)
#
# Atmosphere Composition
#  HI: 1.12e-01
#  OI: 8.88e-01
#
# Parameters:
#   Mp: planet mass (M_earth)
#   Rp: planet radius (R_earth)
#   grav: planet gravity (m/s^2)
#   Teq: equilibrium temperature (K)
#   a: semi-major axis (au)
#   Fxuv: total XUV flux at planet's orbit (Fuv_earth)
#   mstar: stellar mass (Msun)
#   lstar: stellar luminosity (Lsun)
#   Mdot: atmospheric mass-loss rate (kg/s)
#   Mdot_OW17: atmospheric mass-loss rate from Owen & Wu 2017 (kg/s) [IGNORE FOR NOW]
#   ramp_to: 0- ramp success, else- ramp failed
#   run_wind: sim.run_wind() output. 0- wind launched, else- wind failed
#   rochelim: 0- planet outside Roche limit, 1- planet within Roche limit
#   timewall: 0- planet ramped in 20 min., 1- planet/previous solution took longer
#
# Created: October 02, 2026
# Units: SI, astronomical
# Pure H2O atmosphere, no conductive or bolometric heat/cool, static bcs.
##################################################################################\n
Mp\tRp\tgrav\tTeq\ta\tFxuv\tmstar\tlstar\tMdot\tMdotOW17\tramp_to\trun_wind\trochelim\ttimewall
"""

with open(gridname, "w") as file:
    file.write(header)

# iclistname = "grids/incomplete_solution_list.txt"
# with open(iclistname, "w") as file:
#     file.write(f"#Incomplete solution filenames for {gridname}\n")


#file details
timename = f"grids/windae_nobolo_static_info_{date}.dat"
timeheader = """#Output info for grid of windae_nobolo_static_info_100226.dat
#wind-ae inputs and duration
#Mp\tgrav\ta\tFxuv\tmstar\tMdot\tMdotOW17\tramp_to\trun_wind\trochelim\ttimewall\telapsed
"""

with open(timename, "w") as file:
    file.write(timeheader)

# print(mp_axis,rho_axis,teq_axis,mstar_axis,Fxuv_axis)
total_points = dim_mp*dim_g*dim_sep*dim_Fxuv
print('Total combinations to run:',total_points)
print('Maximum run estimate:', (total_points*10)/(60*24), 'days')


# -------------- LOOP THE SIMULATION -------------- #
start_loop = time.time()
count = 0
complete = 0
mp_timewall = False

semimajor_axis = np.array([])
for i,j,k,l in product(range(dim_sep), range(dim_Fxuv), range(dim_g), range(dim_mp)):

    init_time = time.time()
    if l == 0:
        #whenever we begin going through mp again, we reset
        mp_timewall = False

    count += 1
    new_step = time.time()

    # --- set up parameters --- #

    init_params = [mp_axis[l], grav_axis[k], teq_ext[i], ms_ext[i], Fxuv_axis[j]] #for file searching
    mp = init_params[0] *const.Mearth       #g
    grav = init_params[1]                   #cm/s^2
    rp = np.sqrt((mp*G)/(grav))             #cm
    mstar = ms_ext[i] * const.Msun          #g
    lstar = lstar_ext[i] * const.Lsun       #erg/s
    teq = teq_ext[i]                        #K
    semimajor = sep_axis[i]                 #cm
    Fxuv = init_params[4]                   #ergs/cm/s
    T = np.sqrt((4*np.pi**2 * semimajor**3)/(const.G*mstar)) #s
    T_rochelim = np.sqrt(0.593*const.Rearth/rp)*86400 #s
    Mdot_OW17 = solve_mdot(mp*const.Mearth/1000,rp/100,mstar*const.Msun/1000,semimajor/100,3.153e16,lstar*3.828e26/const.Lsun) #kg/s

    #separate array of wind-ae inputs
    target = np.array([mp,rp,mstar,semimajor,Fxuv,lstar]) #can be unpacked with system(*target)    
    savefile = f"{int(grav)}/" + init.make_fname(init_params)
    completefile = init.make_fname(init_params, complete=True)


    # ---- filter solutions --- #
    if T <= T_rochelim:
        rochelim = 1
    else:
        rochelim = 0
    
    #complete solution exists
    if os.path.exists(savefile) is True:
        complete += 1
        target_ramp = 0
        run_wind = 0
        timewall = 0

        if os.path.exists(completefile):
            savefile=completefile
        sim = wind_sim()
        sim.load_planet(savefile, print_atmo=False)
        if rochelim == 0:
            # if NOT inside the roche limit
            Mdot = sim.windsoln.Mdot
        else:
            Mdot = const.Mearth
        print("Solution %s exists." %completefile)

    #if we already hit the time limit
    elif mp_timewall is True:
        target_ramp = 1
        run_wind = -1
        timewall = 1
        Mdot = 0
        print("Skipping Mp=%.2f, g=%.2f m/s^2 due to time limit." %(mp/const.Mearth,grav*0.01))

    #otherwise, run the sim.
    # else:
    #     # print('YAYYYYY')
    #     closest_file, norm_dist = init.initialize(init_params,input_grid)
    #     sim = wind_sim()
        
    #     if norm_dist < 1:
    #         loadfile = closest_file
    #     else:
    #         #Otherwise, load initial guess.
    #         loadfile = f'saves/H2O_neptune_nobolo_static_{grav}.csv'
    #         print("Using initial guess.")


    #     # ---- ramp the planet ---- #
    #     import windae_one as wone #reload it every time to ensure no caching issues

    #     # import subprocess
    #     # yourCommand = "python3 'import windae_one as wone; Mdot, target_ramp, run_wind = wone.run(loadfile, init_params, target)'"
    #     # timeoutSeconds = 2
    #     # subprocess.check_output(yourCommand, shell=False, timeout=timeoutSeconds)

    #     # os.system("timeout 2s python3 -c 'import windae_one as wone; Mdot, target_ramp, run_wind = wone.run(\"%s\", %s, %s)'" %(loadfile, init_params, target))
    #     Mdot, target_ramp, run_wind, timeout = wone.run(loadfile, init_params, target)
    #     # print(Mdot, target_ramp, run_wind)
    #     # break
    #     if timeout:
    #         mp_timewall = True
    #     # if target_ramp == 101:
    #     #     mp_timewall = True
            
    #     if target_ramp + run_wind == 0:
    #         #run is complete and working
    #         complete += 1
    #         savefile = completefile
    #     else:
    #         #log invalid run
    #         with open("grids/invalid_list.dat", 'a') as file:                          
    #             file.write('%.2f %.2i %.2f %.2f %.2f\n' %(init_params[0],init_params[1]/100,init_params[2],init_params[3],init_params[4]))
        

    #     # --- save solution --- #
    #     sim.save_planet(savefile,overwrite=True)
    #     ramp_step = time.time()      
    #     # Mdot = sim.windsoln.Mdot
    #     rochelim = 0
        
    #     #print status
    #     print(" Mass-loss rate: %.4e g/s " %Mdot)
    #     print(" run_wind: ", run_wind)
    #     print(" %i - %.2f min - [%i,%i,%i,%i]\n" %(count,((ramp_step-new_step)/60),l,k,i,j))
    

    # --- save to grids --- #
    massloss_grid[l,k,i,j] = Mdot
    ramp_grid[l,k,i,j]     = target_ramp
    wind_grid[l,k,i,j]     = run_wind

    fin_time = time.time()
    elapsed = (fin_time - init_time)/60 #convert s -> min

    result_line = "%.5f\t%.5f\t%3f\t%i\t%.5f\t%.5e\t%.3f\t%.5f\t%.5e\t%.5e\t%i\t%i\t%i\t%i\n" %(\
        mp/const.Mearth, rp/const.Rearth, grav*0.01, teq, \
        semimajor/const.au, Fxuv/const.FuvEarth, mstar/const.Msun, lstar/const.Lsun, \
        massloss_grid[l,k,i,j]*0.001, Mdot_OW17, ramp_grid[l,k,i,j], wind_grid[l,k,i,j], rochelim, timewall)
    
    result_line_time = "%.5f\t%3f\t%.5f\t%.5e\t%.3f\t%.5e\t%.5e\t%i\t%i\t%i\t%i\t%.5f\n" %(\
        mp/const.Mearth, grav*0.01, semimajor/const.au, Fxuv/const.FuvEarth, mstar/const.Msun, \
        massloss_grid[l,k,i,j]*0.001, Mdot_OW17, ramp_grid[l,k,i,j], wind_grid[l,k,i,j], rochelim, timewall, elapsed)
    
    with open(gridname, 'a') as file:
        file.write(result_line)

    with open(timename, 'a') as file:
        file.write(result_line_time)


end_loop = time.time()
print("Total elapsed: %.2f days" %((end_loop-start_loop)/86400))
print("Average duration: %.2f min" %(((end_loop-start_loop)/60)/count))
# with open("final.txt", "w") as file:
#     file.write("Total elapsed: %.2f days" %((end_loop-start_loop)/86400))
#     file.write("Average duration: %.2f min" %(((end_loop-start_loop)/60)/count))    

print("\nMass-loss grid data saved to", gridname)
print(100 * complete/total_points, "%","solutions complete")
print('\n\nFull grid COMPLETE! Hooray!!!\n')
