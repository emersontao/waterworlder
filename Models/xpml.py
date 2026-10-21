# new_xpml.py

#ALL FUNCTIONS INPUT SI. ALL FUNCTIONS OUTPUT SI.
#Functions that don't are WRAPPED and explicitly stated to be so.

import os, sys, json
sys.path.insert(0, './Models')
import numpy as np
# import matplotlib.pyplot as plt
import swe
import windae_grid_interpolate as wgi

#Constants
M_e   = 5.972e24  #Earth mass (kg)
R_e   = 6.371e6   #Earth radius (m)

M_sun = 1.988e30  #Solar mass (kg)
R_sun = 6.955e10  #Solar radius (m)
L_sun = 3.828e26  #luminosity of the Sun (Watts, kg⋅m^2⋅s^-3)


Fuv_e = 4.667e-3     #W/m^2
au    = 1.496e11  #astronomical unit (m)
Gyr   = 3.153e16  #seconds in 1 billion years
G     = 6.674e-11 #gravitational constant (m^3⋅kg^-1⋅s^-2)
sigma_sb = 5.670374419e-8 #W⋅m−2⋅K−4 ; Stefan-Boltzmann constant



def evolution(mp, rp, teq, wmf, a, mstar, lstar, H2O=True, time_init=1e-3*Gyr, time_end=5*Gyr,):
    """
    mp:         planet mass

    rp:         planet radius
    
    teq:        equilibrium temperature of planet (K)
    
    a:          semi-major axis between planet and host star
    
    mstar:      stellar mass
    
    lstar:      stellar luminosity
    
    H2O:        True: Evolve with water mass-loss rates. False: Evolve with OW17 mass-loss rates. 
    
    time_init:  Age of system at simulation start (s). Default: 10 Myr / 1e-3 Gyr
    
    time_end:   Total duration of evolution (s). Default: 5 Gyr

    Takes in SI units.
    The function that runs the evolution for ONE steam world.
    Assumes top of atmoshpere is at 1 microbar radius.
    """
    ### Setup #########################################################
    #Detect stellar type (M or G)
    #create radius wrapper; converts swe's astro units to take in/return SI
    mstar_cutoff = 0.43
    if mstar/M_sun < mstar_cutoff:
        #M-type host
        def radius(wmf,teq,mp,age):
                    swerad = swe.radius_sweet(np.array(mp/M_e).flatten(), 
                                              np.array(0), np.array(0), 
                                              np.array(wmf), np.array(teq),
                                              np.array(age/Gyr))
                    return swerad * R_e
    else:
        #G-type host
        def radius(wmf,teq,mp,age):
            swerad = swe.radius_sweet(np.array(mp/M_e).flatten(), 
                                        np.array(0), np.array(1), 
                                        np.array(wmf), np.array(teq),
                                        np.array(age/Gyr))
            return swerad * R_e

    #Create mass-loss wrapper
    if H2O:
        def solve_mdot(mstar,teq,fxuv,mp,rp,a):
            grav = (mp * G)/rp**2
            return wgi.mdot(mstar,teq,fxuv,grav,mp)

        def error(mstar,teq,fxuv,mp):
            return wgi.error(mstar, teq, fxuv, grav, mp)

    else:
        def solve_mdot(mstar,teq,fxuv,mp,rp,a):
            return mdot_OW17(mstar,fxuv,mp,rp,a)
        def error(mstar,teq,fxuv,mp):
             return np.nan
             
    #Time inputs
    dt = 1e-4 * Gyr         #timestep
    dt_save = 1e-3 * Gyr    #timestep saving age
    time_since_last_save = 0


    #fixed planet parameters
    mcore = mp*(1-wmf)      #extremely simple calculation here
    grav = (mp * G)/rp**2  #m/s

    #starting planet params
    age  = time_init         #set age to initial time
    fxuv = Fxuv(age, lstar, a)
    mdot = solve_mdot(mstar,teq,fxuv,mp,rp,a) #calculate first mass-loss rate, just in case
    err = error(mstar,teq,fxuv,mp)

    #initialize things to track
    time_saved = []
    mp_saved   = []
    rp_saved   = []
    wmf_saved  = []
    fxuv_saved = []
    mdot_saved = []
    err_saved  = []

    #save starting values
    time_saved.append(age)
    mp_saved.append(mp)
    rp_saved.append(rp)
    wmf_saved.append(wmf)
    fxuv_saved.append(fxuv)
    mdot_saved.append(mdot)
    err_saved.append(err)
    

    while age < time_end:
        #Calculate mass-loss rate
        mdot_old = mdot

        #update timestep
        dtmin = (mp / mdot_old) * 1e-16 #<-------------------------------- !!!! CHANGED: USES mdot_old INSTEAD OF CALCULATING NEW.
        if age > time_init:
            np.seterr(divide='ignore') #ignore the divide by 0 errors
            dtmax = ( 2 * mp ) / np.abs(mdot - mdot_old)
            dt = 10**(0.5 * (np.log10(dtmin)) + 0.5*np.log10(dtmax))
            dt = min(dt,0.5*dt_save)     #stops dt > time_step_save         
        else:
            dt = 10**(0.5 * (np.log10(dtmin)))

        # Possible conditions that would end the loop
        if (dt + age) > time_end:
            dt = time_end - age #ensures the last possible step is 5 Gyr

        if (1 - (mcore/(mp - mdot_old*dt))) < 0.001:  #<------------------- !!!! CHANGED: USES mdot_old INSTEAD OF CALCULATING NEW.
            dt_new = mp*wmf/mdot
            age += dt_new
            mp -= mdot*dt_new
            wmf = 0     #1 - (mcore/mp)
            rp = swe.zeng_radius(mcore/M_e, 0.325) #assume an earth-like core
            break

        # Update state
        age += age
        mp -= mdot*dt
        wmf = 1 - (mcore/mp)
        fxuv = Fxuv(age, lstar, a)
        mdot = solve_mdot(mstar, teq, fxuv, mp, rp, a) #<------------------- We solve for mdot down here because mdot depends on fxuv, which depends on time.
        rp = radius(wmf, teq, mp, age) #intakes all updated parameters

        time_since_last_save += dt

        if any(x <= 0 for x in [mp, rp, mdot]):
            #This should NOT happen.
            print("WARNING: Evolution terminated. 0 in mp, rp, or mdot.")
            break

        if time_since_last_save >= dt_save:
            time_saved.append(age)
            mp_saved.append(mp)
            rp_saved.append(rp)
            wmf_saved.append(wmf)
            fxuv_saved.append(fxuv)
            mdot_saved.append(mdot)
            err_saved.append(err)

            time_since_last_save=0

    #save the last points
    time_saved.append(age)
    mp_saved.append(mp)
    rp_saved.append(rp)
    wmf_saved.append(wmf)
    fxuv_saved.append(fxuv)
    mdot_saved.append(mdot)
    err_saved.append(err)

    return np.array(time_saved), np.array(mp_saved), np.array(rp_saved), np.array(wmf_saved), np.array(fxuv_saved), np.array(mdot_saved), np.array(err_saved)




def planet_generator(n_planets, savepath, params_path='/generator_params.json', printtable=False):
    """
    n_planets:   number of planets to generate

    savepath:    directory that a new directory named '{n_planets}_xp' will be created in and filled with planet info
    
    params_path: loads the generator_params.json file that allows for customization of some generated values. default: '/generator_params.json'
    
    printtable:  print table displaying generated values. Default: False

    Generates a .csv file with randomized parameters for n_planets
    """
    #Load in parameters
    with open(params_path, "r") as file:
        js = json.load(file)

    #Consistent with the evolution script
    age_init = 1e-3*Gyr
    mstar_cutoff = 0.43


    #Begin initializing ranges
    ids = np.arange(1, n_planets+1, 1, dtype=int) #Create an integer ID for each planet


    #Initialize randomized parameters
    mstar_gp = np.random.uniform(js["mstar"][0], js["mstar"][1], n_planets) * M_sun     # kg
    mp_gp = np.random.uniform(js["mp"][0], js["mp"][1], n_planets) * M_e     #M_e
    teq_gp = np.random.uniform(js["teq"][0], js["teq"][1], n_planets)      #K, 750 initially
    wmf_gp = np.random.uniform(js["wmf"][0], js["wmf"][1],n_planets)

    G_startypes = mstar_gp/M_sun > mstar_cutoff #Boolean mask

    #empty lists for now
    rp_gp = []#np.array([])                                #m
    mcore_gp = []#np.array([])                             #M_e
    rsurf_gp = []#np.array([])                             #R_e
    lstar_gp = []#np.array([])                             #W


    #Loop to load arrays with info
    for p in range(n_planets):

        #Core mass
        mcore = mp_gp[p] * (1-wmf_gp[p])
        mcore_gp.append(mcore)

        #"Surface" radius
        rsurf = swe.zeng_radius(mcore/M_e, 0.325)
        rsurf_gp.append(rsurf)

        #Boolean for startyle
        G_star = G_startypes[p]

        if G_star:
            lstar = 0.23 * L_sun * (mstar_gp[p] / M_sun)**2.3
            rp = swe.radius_m(wmf_gp[p],teq_gp[p],mp_gp[p]/M_e,age_init)*R_e #m
        else:
            # Host star is G-TYPE
            lstar = L_sun * (mstar_gp[p] / M_sun)**4
            rp = swe.radius_g(wmf_gp[p],teq_gp[p],mp_gp[p]/M_e,age_init)*R_e #m

        rp_gp.append(rp)
        lstar_gp.append(lstar)

    #exit loop, convert lists to arrays
    rp_gp = np.array(rp_gp)
    mcore_gp  = np.array(mcore_gp)
    rsurf_gp  = np.array(rsurf_gp)
    period_gp = np.array(period_gp)

    #additional arrays
    H_gp = scaleheight(teq_gp, rp_gp, mp_gp)    #initial scaleheight of atmosphere
    lambda_gp = rp_gp/H_gp #Jean's escape parameter
    grav_gp   = mp_gp * G / rp_gp**2
    sep_gp    = a_from_teq(lstar_gp, teq_gp)
    fxuv_gp   = Fxuv(age_init, lstar_gp, sep_gp)
    period_gp = period_solver(mstar_gp, lstar_gp, teq_gp) #days

    print(f'\n{n_planets} planets generated.')

    ### Save data
    #create folder to save evolution data
    savefolder=f'{n_planets}_xp'
    location = savepath+savefolder

    try:
        os.mkdir(location)
        print(f"Directory '{savefolder}' created successfully in {savepath}.")
    except FileExistsError:
        print(f"Directory '{savefolder}' already exists in {savepath}. WARNING: Data in this folder will be overwritten.")
    except Exception as e:
        print(f"An error occurred: {e}")



    #BE EXTRA CAREFUL MODIFYING THESE
    header     = f'Planet/{n_planets}\tMp (kg)\tRp (m)\tCore M (kg)\tSurface Radius (m)\tTeq (K)\tWMF\tgravity(m/s^2)\tscale height (m)\tlambda\tsep (m)\tperiod (days)\tG-Star Type\tMstar (kg)\tLstar (W)\tFxuv (W/m^2)'
    header_e   = f'Planet/{n_planets}\tMp (Me)\tRp (Re)\tCore M (Me)\tSurface Radius (Re)\tTeq (K)\tWMF\tgravity(m/s^2)\tscale height (m)\tlambda\tsep (au)\tperiod (days)\tG-Star Type\tMstar (Msun)\tLstar (Lsun)\tFxuv (Fuv_e)'
    formats    = ['%d','%.5e','%.5e','%.5e','%.5e','%.3f','%.5e','%.5e','%.3f','%.1f','%.5e','%.5f','%i','%.5e','%.5e','%.5e']
    formats_e  = ['%d','%.5f','%.5f','%.5f','%.5f','%.3f','%.5f','%.5f','%.3f','%.1f','%.5f','%.3f','%i','%.5f','%.5f','%.5f']
    filename   = f'{n_planets}_xp_gp_data.csv'
    filename_e = f'{n_planets}_xp_gp_data_unitEarth.csv'

    #save in SI units
    np.savetxt(f'{location}/{filename}', \
               np.c_[ids, mp_gp, rp_gp, mcore_gp, rsurf_gp, teq_gp, wmf_gp, grav_gp, H_gp, lambda_gp, \
                     sep_gp, period_gp, G_startypes, mstar_gp, lstar_gp, fxuv_gp], \
               fmt=formats, delimiter='\t', header=header) 
    print(f'{filename} written to {location}.')
    
    #and in Earth/astronomical units:    
    np.savetxt(f'{location}/{filename_e}', \
               np.c_[ids, mp_gp/M_e, rp_gp/R_e, mcore_gp/M_e, rsurf_gp/R_e, teq_gp, wmf_gp, grav_gp, H_gp, lambda_gp, \
                     sep_gp/au, period_gp, G_startypes, mstar_gp/M_sun, lstar_gp/L_sun, fxuv_gp/Fuv_e], \
               fmt=formats_e, delimiter='\t', header=header_e)
    print(f'{filename_e} written to {location}.')

    if printtable == True:
        # Display planet info (using most intuitive units)
        print('\nPLANET POPULATION INFORMATION')
        from astropy.table import Table
        gp_table = Table()
        gp_table["Number"] = ids
        gp_table["Mass (M𛲜)"] = mp_gp/M_e
        gp_table["Radius (R𛲜)"] = rp_gp/R_e
        gp_table["Surface Radius (R𛲜)"] = rsurf_gp/R_e
        gp_table["Teq (K)"] = teq_gp
        gp_table["WMF"] = wmf_gp
        gp_table["Gravity (m/s^2)"] = grav_gp
        gp_table["Scale Height (km)"] = H_gp/1000
        gp_table["Lambda"] = lambda_gp
        gp_table["Period (days)"] = period_gp
        gp_table["Separation (au)"] = sep_gp/au
        gp_table["G-star Type"] = G_startypes
        gp_table["Lstar (L☉)"] = lstar_gp/L_sun
        gp_table["Mstar (M☉)"] = mstar_gp/M_sun
        gp_table["Fxuv (W/m^2)"] = fxuv_gp

        gp_table["Mass (M𛲜)"].format = "{:.3f}"
        gp_table["Radius (R𛲜)"].format = "{:.3f}"
        gp_table["Surface Radius (R𛲜)"].format = "{:.3f}"
        gp_table["Teq (K)"].format = "{:.3f}"
        gp_table["WMF"].format = "{:.3f}"
        gp_table["Gravity (m/s^2)"].format = "{:.3f}"
        gp_table["Scale Height (km)"].format = "{:.1f}"
        gp_table["Lambda"].format = "{:.3f}"
        gp_table["Period (days)"].format = "{:.1f}"
        gp_table["Separation (au)"].format = "{:.3f}"
        gp_table["Lstar (L☉)"].format = "{:.3f}"
        gp_table["Mstar (M☉)"].format = "{:.3f}"
        gp_table["Fxuv (W/m^2)"].format = "{:.3f}"
        print(gp_table)

    print(f'{n_planets} generated.\n')



### Functions for running in parallel

#Define complete evolution code for generated population
def parallel_evo(n_planets,planet_ID,mp,rp,teq,wmf,a,mstar,lstar,savepath='./',H2O=True):
    '''
    n_planets: Number of planets, int

    planet_ID: planet (out of n_planets), int

    mp:        planet mass (kg)

    rp:        planet radius (m)
    
    teq:       equilibrium temperature (K)
    
    wmf:       water-mass fraction [0,1]
    
    a:         semi-major axis (m)
    
    mstar:     stellar mass (kg)
    
    lstar:     stellar luminosity (W)
    
    savepath:  directory that will have folder 'H2O' or 'H_He' created to store evolution tracks. str
    
    H2O:       True: evolve with water mass-loss rates. False: Evolve with H-He mass-loss rates. boolean
    '''
    ### Save data
    # <--------------------------- !!! Likely will end up putting the folder-maker in xp_evolution. Here for now, it should work
    if H2O:

        #create folder to save evolution data
        savefolder='H2O'
        location = savepath+f'{n_planets}_xp/'+savefolder
        try:
            os.mkdir(location)
            print(f"Directory '{savefolder}' created successfully in {savepath}.")
        except FileExistsError:
            print(f"Directory '{savefolder}' already exists in {savepath}. WARNING: Data in this folder will be overwritten.")
        except Exception as e:
            print(f"An error occurred: {e}")

        time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr, error_arr = evolution(mp,rp,teq,wmf,a,mstar,lstar,True)

        #Save planet evolution data
        planet_filename = f'planet_{int(planet_ID)}-{int(n_planets)}_evolution_data.csv'
        path = savepath + '/' + planet_filename
        np.savetxt(path, np.c_[time_arr, mp_arr, rp_arr, wmf_arr, mdot_arr, error_arr], fmt='%.5e', delimiter='\t')

    else:
        #H2-He evolution
        savefolder='H_He'
        location = savepath+f'{n_planets}_xp/'+savefolder
        try:
            os.mkdir(location)
            print(f"Directory '{savefolder}' created successfully in {savepath}.")
        except FileExistsError:
            print(f"Directory '{savefolder}' already exists in {savepath}. WARNING: Data in this folder will be overwritten.")
        except Exception as e:
            print(f"An error occurred: {e}")

        #run evolution
        time_arr, mp_arr, rp_arr, wmf_arr, fxuv_arr, mdot_arr, error_arr = evolution(mp,rp,teq,wmf,a,mstar,lstar,False)

        #Save planet evolution data
        planet_filename = f'planet_{int(planet_ID)}-{int(n_planets)}_evolution_data.csv'
        path = savepath + planet_filename
        header = f"Planet #{planet_ID} Time [s]\tM_p [kg]\tR_p [m]\tWMF\tF_xuv [W/m^2]\tMdot [kg/s]\tError"
        np.savetxt(path, np.c_[time_arr, mp_arr, rp_arr, wmf_arr, fxuv_arr, mdot_arr, error_arr], fmt='%.5e', delimiter='\t')

        
    return time_arr, mp_arr, rp_arr, wmf_arr, fxuv_arr, mdot_arr, error_arr



def parallel_evo_wrapper(param):
    return parallel_evo(*param)





### Important OW17 Mass-loss functions
def efficiency(mp,rp,mstar,sep):
    """
    Efficiency of photoevaporative atmospheriv H-He mass-loss.
    Owen & Wu 2017, Equation #31
    """
    v_esc = np.sqrt(2*G*mp / rp) #escape velocity at planet surface
    eff = 0.1 * (v_esc / 15e3)**(-2) #mass-loss efficiency, uncorrected
    
    #Roche lobe correction, Salz+ 2016a
    xi = (sep / rp) * (mp / (3*mstar))**(1/3)   #ξ
    k_diff = 1 - 3/(2 * xi) + 1/(2 * xi**3)     #fractional gravitational potential energy difference

    return eff/k_diff #corrected efficiency


def Lxuv(age,lstar):
    """
    high-energy luminosity (X + UV) emitted by the star
    Sanz-Forçada+2011
    SI units
    """
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


def Fxuv(age, lstar, a):
    # SI units
    return Lxuv(age, lstar)*1000/(4*np.pi*a**2)


def mdot_OW17(mstar,fxuv,mp,rp,a):
    '''
    Owen & Wu 2017, Equation #19 equivalent
    Adapted for wind-ae inputs
    SI inputs
    '''
    eta = efficiency(mp,rp,mstar,a)
    mdot = eta * np.pi * (fxuv*rp) / (G*mp) #kg/s
    return mdot






###useful for generator
#planet period
def period_solver(mstar,lstar,teq):
    a = (lstar / (16 * np.pi * sigma_sb * teq**4))**(3/2)
    
    T = 365.25 * np.sqrt(M_sun * a / (mstar * au**3)) #days
    return T


#calculate semi-major axis/equilibrium temperature
def a_from_teq(lstar,teq):
    a = np.sqrt (lstar / (16*np.pi*sigma_sb*teq**4))
    return a

def teq_from_a(lstar,a):
    teq = (lstar / (16*np.pi*sigma_sb*a**2))**0.25
    return teq


#calculate scale height
def scaleheight(teq,rp,mp):
    mmolm = 18.01528 * 0.001 # kg/mol, mean molecular mass of the atmosphere (g/mol*0.001kg/g)
    R = 8.31446 #J/Kmol, molar gas constant
    g = (mp*G)/rp**2 # m/s, gravity of planet
    
    H = (teq*R)/(g*mmolm) # m
    return H