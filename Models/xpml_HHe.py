# eXoPlanet MassLoss Module
#   VERSION:        SPRING 2026
#   LAST UPDATED:   March 04 2026, KAMA
#   INTERIOR:       Aguichine et al. 2025
#   ATMOSPHERE:     Owen and Wu 2017
#   UNITS:          MKS, also exports EARTH
# //////////////////////////////////////// #

import os, sys
sys.path.insert(0, './Models')
import numpy as np
# import matplotlib.pyplot as plt
import swe
# import windae_grid_interpolate as wgi

#constants and astronomical units
M_e   = 5.972e24  #Earth mass (kg)
R_e   = 6.371e6   #Earth radius (m)

M_sun = 1.988e30  #Solar mass (kg)
R_sun = 6.955e10  #Solar radius (m)
L_sun = 3.828e26  #luminosity of the Sun (Watts, kg⋅m^2⋅s^-3)
T_sun = 5778      #Solar temperature (K)


Fuv_e = 4.667     #W/m^2
au    = 1.496e11  #astronomical unit (m)
Gyr   = 3.153e16  #seconds in 1 billion years
G     = 6.674e-11 #gravitational constant (m^3⋅kg^-1⋅s^-2)
sigma_sb = 5.670374419e-8 #W⋅m−2⋅K−4 ; Stefan-Boltzmann constant


#Evolution for a single exoplanet (MKS)
#  Returns evolution file tracking age, mass, radius, wmf, mass-loss rate
def evolution(mp,rp,teq,wmf,sep,mstar,lstar):
    """
    mp: planet mass [kg]
    rp: planet radius [m]
    teq: equilibrium temperature [K]
    wmf: water-mass fraction
    sep: planet-star separation [m]
    mstar: stellar mass [kg]
    lstar: stellar luminosity [W]
    """
    #Input units: MKS

    #time values
    time_init = 1e-3 * Gyr      #1 Myr
    time_end = 5 * Gyr
    delta_time = 1e-4 * Gyr
    time_step_save = 1e-3 * Gyr #default: 1e-3 Gyr
    time_since_last_save = 0
    age = time_init #s

    #L*-dependent wrappers
    #https://en.wikipedia.org/wiki/Stellar_classification
    #0.43 comes from luminosity-relations
    mstar_cutoff = 0.43     #0.5585 -> midpoint of K-type min/max bolometric luminosities, solar masses. 
    if mstar/M_sun < mstar_cutoff:
        def radius(wmf,teq,mp,age):
            swerad = swe.radius_sweet(np.array(mp/M_e).flatten(),
                                      np.array(0), np.array(0), 
                                      np.array(wmf), np.array(teq),
                                      np.array(age/Gyr))
            return swerad * R_e
    else:
        #Host star is G-TYPE
        def radius(wmf,teq,mp,age):
            swerad = swe.radius_sweet(np.array(mp/M_e).flatten(), 
                                      np.array(0), np.array(1), 
                                      np.array(wmf), np.array(teq),
                                      np.array(age/Gyr))
            return swerad * R_e
        
    mcore = mp * (1 - wmf)
    rsurf = swe.zeng_radius(mcore/M_e,0.325)
    grav = gravity(mp,rp)
    F_xuv = find_Fxuv(age,lstar,sep)
    # print(mcore)
    mdot = solve_mdot(mp, rp, mstar, sep, age, lstar)

    #create arrays
    time_saved = np.array([])
    mp_saved = np.array([])
    wmf_saved = np.array([])
    rp_saved = np.array([])
    mdot_saved = np.array([])

    #add the starting values
    time_saved = np.append(time_saved, age)
    mp_saved = np.append(mp_saved, mp)
    wmf_saved = np.append(wmf_saved, wmf)
    rp_saved = np.append(rp_saved, rp)
    mdot_saved = np.append(mdot_saved, mdot)
    iterations = 0

    # t_start = time.time()

    #simulation
    while age < time_end :
        #calculate timestep here
        mdot_old = mdot

        mdot = solve_mdot(mp, rp, mstar, sep, age, lstar)
        # print('mstar:', mstar/M_sun,'T_eq:',teq,'Fxuv:', F_xuv/Fuv_e, 'Mdot:', mdot, 'Errcode:', err)
        # break
        
        dtmin = (mp / mdot) * 1e-16
        if iterations > 0:
            np.seterr(divide='ignore') #ignore the divide by 0 errors
            dtmax = ( 2 * mp ) / np.abs(mdot - mdot_old)
            delta_time = 10**(0.5 * (np.log10(dtmin)) + 0.5*np.log10(dtmax))
            delta_time = min(delta_time,0.5*time_step_save)     #stops dt > time_step_save            
        else:
            delta_time = 10**(0.5 * (np.log10(dtmin)))

        if (delta_time + age) > time_end:
            #timestep failsafe
            delta_time = time_end - age


        #WMF Safety condition
        expected_wmf = 1 - (mcore/(mp - mdot*delta_time))
        if expected_wmf < 0.001:
            #Will fully deplete remaining 0.1% of water in ONE timestep.
            dt_new = (mp*wmf)/mdot
            # print(dt_new)
            # mdot_original = mdot

            age = age + dt_new
            mp = float(mp - mdot*dt_new)
            wmf = 1 - (mcore/mp) #which should be 0 now
            rp = rsurf

            time_saved = np.append(time_saved, age)
            mp_saved = np.append(mp_saved, mp)
            rp_saved = np.append(rp_saved, rp)
            wmf_saved = np.append(wmf_saved, wmf)
            mdot_saved = np.append(mdot_saved, mdot)
            # diditwork = 'Yes'
            # # print('wmf < 0 final rp:',rp)
            break
        
        #Update state
        age = age+delta_time
        mp = float(mp - mdot*delta_time)
        wmf = 1 - (mcore/mp)
        rp = radius(wmf,teq,mp,age)

        grav = gravity(mp,rp)
        F_xuv = find_Fxuv(age,lstar,sep)
        time_since_last_save += delta_time
        iterations += 1

        #Break conditions
        if mp < 0 :
            break
        if rp < 0 :
            break
        if mdot < 0:
            break
        if mdot == 0:
            break
        #save to arrays every time step
        if time_since_last_save >= time_step_save :
            time_saved = np.append(time_saved, age)
            mp_saved = np.append(mp_saved, mp)
            rp_saved = np.append(rp_saved, rp)
            wmf_saved = np.append(wmf_saved, wmf)
            mdot_saved = np.append(mdot_saved, mdot)
            time_since_last_save = 0     #reset time since last save
        

    #add final points
    time_saved = np.append(time_saved, age)
    mp_saved = np.append(mp_saved, mp)
    wmf_saved = np.append(wmf_saved, wmf)
    rp_saved = np.append(rp_saved, rp)
    mdot_saved = np.append(mdot_saved, mdot)
    # print(reason)
    # end_t = time.time()
    #print(f'Runtime: {end_t - t_start:.2f}s')

    #save data
    # header = f'time [s]\tmass [kg]\twmf\tradius [m]\tmdot [kg]'
    # np.savetxt(f'./evolution-data/gp-evo/5000_xp/planet_1-1_evolution_data_wmf2.csv', \
    #         np.c_[time_saved, mp_saved, rp_saved, wmf_saved, mdot_saved], \
    #         fmt='%.9f', delimiter='\t', header=header)
    # print('File saved')

    return time_saved, mp_saved, rp_saved, wmf_saved, mdot_saved


#calculate planet period
def period_solver_m(mstar,teq):
    Ls = L_sun*0.23*(mstar/M_sun)**2.3
    a = (Ls / (16 * np.pi * sigma_sb * teq**4))**(3/2)
    
    T = 365.25 * np.sqrt(M_sun * a / (mstar * au**3)) #days
    return T

def period_solver_g(mstar,teq):
    Ls = L_sun*(mstar/M_sun)**4
    a = (Ls/(16*np.pi*sigma_sb*teq**4))**(3/2)
    T = 365.25 * np.sqrt(M_sun * a / (mstar * au**3)) #days
    return T


#calculate planet gravity
def gravity(mp,rp):
    grav = (mp * G)/(rp**2)  #m/s
    return grav


#calculate scale height
def scaleheight(teq,rp,mp):
    mmolm = 18.01528 * 0.001 # kg/mol, mean molecular mass of the atmosphere (g/mol*0.001kg/g)
    R = 8.31446 #J/Kmol, molar gas constant
    g = gravity(mp,rp) # m/s, gravity of planet
    
    H = (teq*R)/(g*mmolm) # m
    return H


#calculate semi-major axis/equilibrium temperature
def sep_from_teq(lstar,teq):
    a = np.sqrt (lstar / (16*np.pi*sigma_sb * teq**4))
    return a

def teq_from_sep(lstar,sep):
    teq = (lstar / (16*np.pi*sigma_sb * sep**2))**0.25
    return teq


#find F_xuv at planet orbit
def find_Fxuv(age,lstar,sep):
    '''
    age:    system age [s]
    lstar:  stellar luminosity [W]
    sep:    semi-major axis [m]
    '''
    #SI units
    Fxuv =  find_Lxuv(age,lstar) / (4*np.pi*sep**2)

    return Fxuv


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

def find_Lxuv(age,lstar):
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
    L_he = find_Lxuv(age,lstar)
    eta = efficiency(mp,rp,mstar,sep)
    mdot = eta * (L_he*rp**3) / (4*G*mp*sep**2)

    return mdot
