# gridmaker.py

from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
# from wind_ae.wrapper.wrapper_utils.plots import energy_plot, six_panel_plot, quick_plot 
from wind_ae.wrapper.wrapper_utils import constants as const
# from wind_ae.wrapper.wrapper_utils.system import system

import os
import numpy as np
import matplotlib.pyplot as plt

from itertools import product
import initialize_weights as init

#Grid setup
mp_axis     = np.array([0.2, 2, 4, 7, 11, 15, 18, 20, 22.5, 25, 27, 30]) #earth mass
grav_axis   = np.array([515, 724, 1400, 3000])
teq_axis    = np.array([400, 675, 950, 1225, 1500]) #K
mstar_axis  = np.array([0.1, 0.22, 0.3, 0.55, 1]) #sun mass
Fxuv_axis = init.log_Fxuv_axis(mstar_axis,teq_axis,5)

# teq_temp = [400, 675, 950]
# dimensions
dim_mp    = len(mp_axis)
dim_grav  = len(grav_axis)
dim_teq   = len(teq_axis)
dim_mstar = len(mstar_axis)
dim_Fxuv  = len(Fxuv_axis)

total_points = dim_mp*dim_grav*dim_teq*dim_mstar*dim_Fxuv
print("total points:", total_points)



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
gridname = f"grids/windae_nobolo_static_102325.dat"

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
#   Mdot_OW17: atmospheric mass-loss rate from Owen & Wu 2017 (kg/s)
#   ramp_to: 0- ramp success, else- ramp failed
#   run_wind: sim.run_wind() output. 0- wind launched, else- wind failed
#   rochelim: 0- planet outside Roche limit, 1- planet within Roche limit
#
# Created: October 23, 2025
# Units: SI, astronomical
# Pure H2O atmosphere, no conductive or bolometric heat/cool, static bcs.
##################################################################################\n
Mp\tRp\tgrav\tTeq\ta\tFxuv\tmstar\tlstar\tMdot\tMdotOW17\tramp_to\trun_wind\trochelim
"""

with open(gridname, "w") as file:
    file.write(header)

# iclistname = "grids/incomplete_solution_list.txt"
# with open(iclistname, "w") as file:
#     file.write(f"#Incomplete solution filenames for {gridname}\n")

#build lists
params = [(mp, grav, teq, mstar, Fxuv) 
          for mp in mp_axis 
          for grav in grav_axis 
          for teq in teq_axis
          for mstar in mstar_axis 
          for Fxuv in Fxuv_axis
          ]

noexist = 0
complete = 0
loops = 0
for x in params:
    loops +=1
    rochelim = 0
    # print('\n',x)
    [mp, grav, teq, mstar, Fxuv] = x

    if mstar < 0.43:
        Lstar = 0.23*const.Lsun*(mstar)**2.3
    else:
        Lstar = const.Lsun*(mstar)**4
    semimajor = np.sqrt(Lstar / (16*np.pi* 5.67037e-5 * teq**4)) #cm
    Fxuv = x[4]                   #ergs/cm/s
    Rp = np.sqrt((mp*const.Mearth*const.G)/(grav)) #cm
    T = np.sqrt((4*np.pi**2 * semimajor**3)/(const.G*mstar*const.Msun)) #s
    T_rochelim = np.sqrt(0.593*const.Rearth/Rp)*86400 #s

    Mdot_OW17 = solve_mdot(mp*const.Mearth/1000,Rp/100,mstar*const.Msun/1000,semimajor/100,3.153e16,Lstar*3.828e26/const.Lsun) #kg/s


    #if period is within roche limit
    if T <= T_rochelim:
        rampto = 1
        rw = 1
        rochelim = 1
        Mdot = const.Mearth
        result_line = "%.5f\t%.5f\t%3f\t%i\t%.5f\t%.5e\t%.3f\t%.5f\t%.5e\t%.5e\t%i\t%i\t%i\n" %(mp,Rp/const.Rearth,grav*0.01,teq,semimajor/const.au,Fxuv/const.FuvEarth,mstar,Lstar/const.Lsun,Mdot*0.001,Mdot_OW17,rampto,rw,rochelim)
        # print("Invalid Planet: period too short (%.2f days)." %(T/(86400)))
        with open(gridname, 'a') as file:                          
            file.write(result_line)

    else:
        file = init.make_fname(x)
        completefile = init.make_fname(x,complete=True)
        sim = wind_sim()

        rp = Rp / const.Rearth
        a = semimajor / const.au
        lstar = Lstar / const.Lsun

        if os.path.exists(completefile) is True:
            complete +=1

            sim.load_planet(completefile,print_atmo=False,print_warnings=False)

            Mdot = sim.windsoln.Mdot
            rampto = 0
            rw = sim.run_wind()

            result_line = "%.5f\t%.5f\t%3f\t%i\t%.5f\t%.5e\t%.3f\t%.5f\t%.5e\t%.5e\t%i\t%i\t%i\n" %(mp,Rp/const.Rearth,grav*0.01,teq,a,Fxuv/const.FuvEarth,mstar,Lstar/const.Lsun,Mdot*0.001,Mdot_OW17,rampto,rw,rochelim)
            with open(gridname, 'a') as file:                          
                file.write(result_line)
        
        else:
            try:
                sim.load_planet(file,print_atmo=False,print_warnings=False)

                Mdot = sim.windsoln.Mdot
                rampto = 1
                rw = sim.run_wind()

                result_line = "%.5f\t%.5f\t%3f\t%i\t%.5f\t%.5e\t%.3f\t%.5f\t%.5e\t%.5e\t%i\t%i\t%i\n" %(mp,Rp/const.Rearth,grav*0.01,teq,a,Fxuv/const.FuvEarth,mstar,Lstar/const.Lsun,Mdot*0.001,Mdot_OW17,rampto,rw,rochelim)
                with open(gridname, 'a') as file:                          
                    file.write(result_line)

            except FileNotFoundError:
                noexist +=1
                # mp,rp,a,Fxuv,mstar,lstar,Mdot,rampto,rw = [-1,-1,-1,-1,-1,-1,-1,-1,-1]
    #

print("\nMass-loss grid data saved to", gridname)
print(100 * complete/loops, "%","solutions complete")
if noexist >0:
    print("WARNING: SOME FILES NOT FOUND.", noexist, "FILES MISSING.")