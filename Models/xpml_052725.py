# eXoPlanet MassLoss Module
#   VERSION: SPRING 2025
#   LAST UPDATED: MAY 27 2025, KAMA
#   INTERIOR:   Aguichine et al. 2025

import sys, time
sys.path.insert(0, './Models')
import numpy as np
import matplotlib.pyplot as plt
import swe  #new steam world evolution model
import module_mass_loss as mml # type: ignore


def evolution(mp,rp,teq,wmf,sep,mstar,lstar,planetnum):
    # Input units: mkgs

    # time values
    time_init = 1e-2 * mml.Gyr # 10 Myr ... 1 Gyr = 3.15e16 s
    time_end = 20 * mml.Gyr
    delta_time = 1e-4 * mml.Gyr
    time_step_save = 1e-3 * mml.Gyr     #default: 1e-3 Gyr
    time_since_last_save = 0
    age = time_init #s
    
    # Stellar luminosity-dependent function wrappers
    # rp (swe_radius), wmf (solve_wmf)
    # TEMPORARY reference: https://en.wikipedia.org/wiki/Stellar_classification
    mstar_cutoff = 0.56 * mml.M_sun #midpoint of K-type min/max bolometric luminosities
    
    if mstar > mstar_cutoff:
        # Host star is G-TYPE
        def swe_radius(wmf,teq,mp,age):
            return swe.radius_g(wmf,teq,mp,age) #m
        
        # def solve_wmf(rp,teq,mp,age):
        #     return swe.solve_wmf_g(rp,teq,mp,age)

    else:
        # Host star is M-TYPE
        def swe_radius(wmf,teq,mp,age):
            return swe.radius_m(wmf,teq,mp,age) #m
        
        # def solve_wmf(rp,teq,mp,age):
        #     return swe.solve_wmf_m(rp,teq,mp,age)

    # # Only calculate wmf if it is NOT given
    # wmf = solve_wmf(rp,teq,mp,age)
    mcore = mp * (1 - wmf) #kg
    mdot = mml.mdot(mp,rp,mstar,sep,age,lstar)

    #create arrays
    time_saved = np.array([])
    mp_saved = np.array([])
    wmf_saved = np.array([])
    rp_saved = np.array([])
    
    #add the starting values
    time_saved = np.append(time_saved, age)
    mp_saved = np.append(mp_saved, mp)
    wmf_saved = np.append(wmf_saved, wmf)
    rp_saved = np.append(rp_saved, rp)
    iterations = 0

    t_start = time.time()

    #simulation
    while age < time_end :
        
        #calculate timestep here
        mdot_old = mdot
        mdot = mml.mdot(mp,rp,mstar,sep,age,lstar)
    
        dtmin = (mp / mdot) * 1e-16
        #Change this, if within if is not good
        if iterations > 0:
            np.seterr(divide='ignore') #ignore the divide by 0 errors
            dtmax = ( 2 * mp ) / ( np.abs(mdot - mdot_old) )
            delta_time = 10**((0.5*(np.log10(dtmin))) + 0.5*(np.log10(dtmax)))
            delta_time = min(delta_time,0.5*time_step_save)     # stops dt > time_step_save            
        else: 
            delta_time = 10**(0.5*(np.log10(dtmin)))

        #simulation
        age = age + delta_time
        mp = mp - mdot*delta_time
        wmf = 1 - (mcore/mp)
        rp = swe_radius(wmf,teq,mp,age)

        time_since_last_save += delta_time
        iterations += 1
        
        #breaking conditions
        if mp < 0 :
            # print('mp break')
            break
        if rp < 0 :
            # print('rp break')
            break
        if wmf < 0 :
            # print('wmf break')
            break
        if mdot < 0 :
            # print('mdot < 0 break')
            break
        if mdot == 0 :
            # print('mdot = 0 break')
            break

        #save to arrays every time step
        if time_since_last_save >= time_step_save :
            time_saved = np.append(time_saved, age)
            mp_saved = np.append(mp_saved, mp)
            wmf_saved = np.append(wmf_saved, wmf)
            rp_saved = np.append(rp_saved, rp)
            time_since_last_save = 0     #reset time since last save

    #add final points
    time_saved = np.append(time_saved, age)
    mp_saved = np.append(mp_saved, mp)
    wmf_saved = np.append(wmf_saved, wmf)
    rp_saved = np.append(rp_saved, rp)
    #print('Number of iterations:',iterations)
    
    # for i, rp in enumerate (rp_saved):      #calls out all values
    #     print(i, rp)
    
    end_t = time.time()
    #print(f'Runtime: {end_t - t_start:.2f}s')

    # #save data
    # header = f'{planet}\ttime [Gyr]\tmass [Me]\twmf\tradius [Re]'
    # np.savetxt(f'./evolution-data/10Myr_{planet}_swerp_mass_loss_data.csv', \
    #         np.c_[time_saved, mp_saved, wmf_saved, rp_saved], \
    #         fmt='%.12f', delimiter='   ', header=header)
    # print('File saved')

    return time_saved, mp_saved, rp_saved, wmf_saved
    


# Uniformly generated steam sub-Neptune planets for n_planets
def planet_generator(n_planets):
    
    planetnum = np.arange(1,n_planets+1,1,dtype=int) #Assign a number to each random planet

    # Mstar Parameters drawn from Cifuentes et al. (2007) https://arxiv.org/pdf/2007.15077 (Table 6)
    mstar_gp = np.random.uniform(0.077,1.040,n_planets) * mml.M_sun     # kg
    lstar_gp = np.array([])
    for k in range(n_planets):
        #using https://en.wikipedia.org/wiki/Mass%E2%80%93luminosity_relation
        
        if mstar_gp[k]/mml.M_sun < 0.43: 
            lstar = 0.23 * mml.L_sun * ( mstar_gp[k] / mml.M_sun )**2.3
            lstar_gp = np.append(lstar_gp,lstar)
            
        elif mstar_gp[k]/mml.M_sun < 2:
            lstar = mml.L_sun * ( mstar_gp[k] / mml.M_sun )**4
            lstar_gp = np.append(lstar_gp,lstar)

    # PLANET PARAMETERS
    age_init = 1e-2 * mml.Gyr # 10 Myr; consistent age
    mp_gp = np.random.uniform(7,12,n_planets) * mml.M_e    # kg
    teq_gp = np.random.uniform(750,1200,n_planets)         # K
    wmf_gp = np.random.uniform(0.001,0.6,n_planets)
    sep_gp = mml.get_sep_from_teq(lstar_gp,teq_gp)         # m

    rp_gp = np.array([])

    # TEMPORARY reference: https://en.wikipedia.org/wiki/Stellar_classification
    # Eventually use https://articles.adsabs.harvard.edu/pdf/1981A%26AS...46..193H
    mstar_cutoff = 0.5585 * mml.M_sun #midpoint of K-type min/max bolometric luminosities
    for j in range(n_planets):
        # Identify stellar type and relevant rp function.
        
        if mstar_gp[j] > mstar_cutoff:
        # Host star is G-TYPE
            def swe_radius(wmf,teq,mp,age):
                return swe.radius_g(wmf,teq,mp,age) #m    
        else:
            # Host star is M-TYPE
            def swe_radius(wmf,teq,mp,age):
                return swe.radius_m(wmf,teq,mp,age) #m

        # calculate radii
        rp = swe_radius(wmf_gp[j],teq_gp[j],mp_gp[j],age_init)
        print(j, rp)
        rp_gp = np.append(rp_gp, rp)

    period_gp = period(mstar_gp,mp_gp,sep_gp)
    # print(period_gp)
    H_gp = scaleheight(teq_gp,rp_gp,mp_gp)
    # print(H_gp)

    print(f'\n{n_planets} planets generated.\n')

    # //////////////////////////////////////// #
    # SAVE STARTING DATA

    # header = f'{n_planets}\tMass \tRadius \tTeq [K]\tWMF [1]\tSep [au]\tMstar [kg]\tLstar [W]'
    header = f'{n_planets} Planets \tMp (Me), Rp (Re), Teq (K), WMF, sep (au), period (days), scale height (m), Mstar (Msun), Lstar (Lsun)'
    savepath_init = f'evolution-data/gp_evo/{n_planets}_xp'
    filename = f'{n_planets}_xp_gp_data.csv'
    np.savetxt(f'{savepath_init}/{filename}', np.c_[planetnum, 
            mp_gp/mml.M_e, rp_gp/mml.R_e, teq_gp, wmf_gp, sep_gp/mml.au, period_gp, H_gp,
            mstar_gp/mml.M_sun, lstar_gp/mml.L_sun], 
            fmt='%.5f', delimiter='\t', header=header)
    # Saved units: Earth
    
    # Return units: mkgs
    return mp_gp, rp_gp, teq_gp, wmf_gp, sep_gp, period_gp, H_gp, mstar_gp, lstar_gp



# calculate period
def period(mstar,mp,sep):
    day_s = 60*60*24 # s/day

    C = (4 * np.pi**2 )/(mml.Ggrav * (mstar + mp)) #constant, m^3/s^2
    T = np.sqrt((sep**3 * C)) / day_s

    return T
    


def scaleheight(Teq,rp,mp):
    G = mml.Ggrav
    rp = rp*mml.R_e #convert to m
    mp = mp*mml.M_e #convert to kg

    # kb = 1.381e-23 # J/K, Boltzmann constant
    mmolm = 18.01528 # kg/mol, mean molecular mass of the atmosphere; PURE STEAM here
    R = 8.31446 #J/Kmol, molar gas constant
    
    gravp = mp*G/rp**2 #gravity of planet
    print(gravp)

    # calculate scale height, km?
    H = (Teq*R)/(gravp*mmolm)
    
    return H