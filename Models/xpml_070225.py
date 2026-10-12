# eXoPlanet MassLoss Module
#   VERSION:        SPRING 2025
#   LAST UPDATED:   JUNE 30 2025, DESKTOP
#   INTERIOR:       Aguichine et al. 2025
#   UNITS:          EARTH

import os, sys
sys.path.insert(0, './Models')
import numpy as np
import matplotlib.pyplot as plt
import swe # type: ignore
import module_mass_loss as mml # type: ignore


def evolution(mp,rp,teq,wmf,sep,mstar,lstar):
    # Input units: Earth, needs to use MKS
    mp = mp*mml.M_e
    rp = rp*mml.R_e
    sep = sep*mml.au
    mstar = mstar*mml.M_sun
    lstar = lstar*mml.L_sun

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
            return swe.radius_g(wmf,teq,mp/mml.M_e,age/mml.Gyr)*mml.R_e
        
        # def solve_wmf(rp,teq,mp,age):
        #     return swe.solve_wmf_g(rp,teq,mp,age)

    else:
        # Host star is M-TYPE
        def swe_radius(wmf,teq,mp,age):
            return swe.radius_m(wmf,teq,mp/mml.M_e,age/mml.Gyr)*mml.R_e
        
        # def solve_wmf(rp,teq,mp,age):
        #     return swe.solve_wmf_m(rp,teq,mp,age)

    # # Only calculate wmf if it is NOT given
    # wmf = solve_wmf(rp,teq,mp,age)
    mcore = mp * (1 - wmf)  #non-water mass
    mdot = mml.mdot(mp,rp,mstar,sep,age,lstar)
        
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
        mdot = mml.mdot(mp,rp,mstar,sep,age,lstar)

        if mdot < 0 :
            # print('mdot<0 break')
            break
        if mdot == 0 :
            # print('mdot==0 break')
            break
        
        # if mdot <= 0:
        #     # print('mdot<=0 break')
        #     time_saved = np.append(time_saved, age)
        #     mp_saved = np.append(mp_saved, mp)
        #     wmf_saved = np.append(wmf_saved, wmf)
        #     rp_saved = np.append(rp_saved, rp)
        #     mdot_saved = np.append(mdot_saved, mdot)
        #     break
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
            
        #save to arrays every time step
        if time_since_last_save >= time_step_save :
            time_saved = np.append(time_saved, age)
            mp_saved = np.append(mp_saved, mp)
            rp_saved = np.append(rp_saved, rp)
            wmf_saved = np.append(wmf_saved, wmf)
            mdot_saved = np.append(mdot_saved, mdot)
            time_since_last_save = 0     #reset time since last save

    
    # Commented because SWE can't handle wmf=0
    #add final points
    time_saved = np.append(time_saved, age)
    mp_saved = np.append(mp_saved, mp)
    wmf_saved = np.append(wmf_saved, wmf)
    rp_saved = np.append(rp_saved, rp)
    mdot_saved = np.append(mdot_saved, mdot)
    #print('Number of iterations:',iterations)
    
    # for i, rp in enumerate (rp_saved):      #calls out all values
    #     print(i, rp)
    
    # end_t = time.time()
    #print(f'Runtime: {end_t - t_start:.2f}s')

    # #save data
    # header = f'{planet}\ttime [Gyr]\tmass [Me]\twmf\tradius [Re]'
    # np.savetxt(f'./evolution-data/10Myr_{planet}_swerp_mass_loss_data.csv', \
    #         np.c_[time_saved, mp_saved, wmf_saved, rp_saved], \
    #         fmt='%.12f', delimiter='   ', header=header)
    # print('File saved')

    # time_saved = time_saved/mml.Gyr
    # mp_saved = mp_saved/mml.M_e
    # rp_saved = rp_saved/mml.R_e
    # mdot_saved = mdot_saved/mml.M_e

    return time_saved, mp_saved, rp_saved, wmf_saved, mdot_saved



# Uniformly generated steam sub-Neptune planets for n_planets
# Makes file including:planetnum, mp_gp, rp_gp, core_m_gp, teq_gp, wmf_gp, H_gp, lambda_gp, sep_gp, period_gp, startype_gp, mstar_gp, lstar_gp
#   planet index, mass, radius, core mass, equilibrium temperature, water-mass fraction, scale height, jean's escape rate,
#   distance from star, orbital period, star type (Mtype:0, Gtype:1), host star mass, host star luminosity
def planet_generator(n_planets):
    
    planetnum = np.arange(1,n_planets+1,1,dtype=int) #Assign a number to each random planet
    
    age_init = 1e-2 #* mml.Gyr # 10 Myr; consistent age
    
    # STELLAR PARAMETERS
    # Mstar Parameters drawn from Cifuentes et al. (2007) https://arxiv.org/pdf/2007.15077 (Table 6)
    mstar_gp = np.random.uniform(0.077,1.040,n_planets) #* mml.M_sun     # kg
    lstar_gp = np.array([])
    startype_gp = np.zeros(n_planets)

    for k in range(n_planets):
        #using https://en.wikipedia.org/wiki/Mass%E2%80%93luminosity_relation
        # Luminosity of M and G
        if mstar_gp[k] < 0.43: 
            lstar = 0.23 * (mstar_gp[k])**2.3
            lstar_gp = np.append(lstar_gp,lstar)
            
        elif mstar_gp[k] < 2:
            lstar = mstar_gp[k] **4
            lstar_gp = np.append(lstar_gp,lstar)

    # PLANET PARAMETERS
    mp_gp = np.random.uniform(1,12,n_planets)   # Earth Masses
    core_m_gp = mp_gp*0.325
    teq_gp = np.random.uniform(750,1200,n_planets)  # K
    wmf_gp = np.random.uniform(0.001,0.6,n_planets)
    sep_gp = mml.get_sep_from_teq(lstar_gp*mml.L_sun,teq_gp)/mml.au # au
    rp_gp = np.array([])
    # core_r_gp = np.array([])
    # TEMPORARY reference: https://en.wikipedia.org/wiki/Stellar_classification, eventually use https://articles.adsabs.harvard.edu/pdf/1981A%26AS...46..193H
    mstar_cutoff = 0.5585 #midpoint of K-type min/max bolometric luminosities, solar masses
    for j in range(n_planets):
        # Identify stellar type and relevant rp function.
        
        if mstar_gp[j] > mstar_cutoff:
            # Host star is G-TYPE
            def swe_radius(wmf,teq,mp,age):
                return swe.radius_g(wmf,teq,mp,age) #m
            
            # Mark star type. [1]
            startype_gp[j] = 1

        else:
            # Host star is M-TYPE
            def swe_radius(wmf,teq,mp,age):
                return swe.radius_m(wmf,teq,mp,age) #m

        # calculate radii
        rp = swe_radius(wmf_gp[j],teq_gp[j],mp_gp[j],age_init)
        # print(j, rp)
        rp_gp = np.append(rp_gp, rp)

        '''#core radius; NOT RELIABLE.
        core_r = swe_radius(0.001,teq_gp[j],core_m_gp[j],age_init) #smallest wmf allowed for swe model
        # print(core_r)
        core_r_gp = np.append(core_r_gp,core_r)
        # print(core_r_gp)'''

    period_gp = period(mstar_gp,mp_gp,sep_gp)   #days
    H_gp = scaleheight(teq_gp,rp_gp,mp_gp) #m
    lambda_gp = (rp_gp*mml.R_e)/H_gp
    grav_gp = grav(mp_gp,rp_gp)     #m/s^2

    print(f'\n{n_planets} planets generated.\n')


    # //////////////////////////////////////// #
    # SAVE GENERATED DATA

    # Location details
    savepath = 'evolution-data/gp-evo/'
    datafolder = f'{n_planets}_xp'
    location = savepath + datafolder
    
    # Create folder to save evolution data
    try:
        os.mkdir(location)
        print(f"Directory '{datafolder}' created successfully in {savepath}.")
    except FileExistsError:
        print(f"Directory '{datafolder}' already exists in {savepath}.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # header = f'{n_planets}\tMass \tRadius \tTeq [K]\tWMF [1]\tSep [au]\tMstar [kg]\tLstar [W]'
    # info = f'{n_planets} Planets - Earth units - 1-12 Mp, 0.001-0.6 WMF, 750-1200K, M/G stars'
    header = f'Planet/{n_planets}\tMp (Me)\tRp (Re)\tCore M (Me)\tTeq (K)\tWMF\tscale height (m)\tlambda\tsep (au)\tperiod (days)\tStar Type\tMstar (Msun)\tLstar (Lsun)'
    # header = info + labels
    filename = f'{n_planets}_xp_gp_data.csv'
     #path for each individual set of data
    
    #use MKS units
    np.savetxt(f'{location}/{filename}', 
               np.c_[planetnum, mp_gp, rp_gp, core_m_gp, teq_gp, wmf_gp, H_gp, lambda_gp, sep_gp, period_gp, startype_gp, mstar_gp, lstar_gp], 
               fmt='%.5f', delimiter='\t', header=header)
    
    # #Save in Earth/Sol units:
    # np.savetxt(f'{savepath_init}/{filename}', np.c_[planetnum, mp_gp/mml.M_e, rp_gp/mml.R_e, teq_gp, wmf_gp, sep_gp/mml.au, period_gp, H_gp, mstar_gp/mml.M_sun, lstar_gp/mml.L_sun], 
    #         fmt='%.5f', delimiter='\t', header=header)
    
    # Display planet info (Earth units)
    print('PLANET POPULATION INFORMATION')
    from astropy.table import Table
    gp_table = Table()
    gp_table["Number"] = planetnum
    gp_table["Mass (M🜨)"] = mp_gp
    gp_table["Radius (R🜨)"] = rp_gp
    gp_table["Teq (K)"] = teq_gp
    gp_table["WMF"] = wmf_gp
    gp_table["Gravity (m/s^2)"] = grav_gp
    gp_table["Scale Height (km)"] = H_gp/1000
    gp_table["Lambda"] = lambda_gp
    gp_table["Period (days)"] = period_gp
    gp_table["Sep (au)"] = sep_gp
    gp_table["Star Type"] = startype_gp
    gp_table["Lstar (L☉)"] = lstar_gp
    gp_table["Mstar (M☉)"] = mstar_gp
    
    gp_table["Mass (M🜨)"].format = "{:.3f}"
    gp_table["Radius (R🜨)"].format = "{:.3f}"
    gp_table["Teq (K)"].format = "{:.3f}"
    gp_table["WMF"].format = "{:.3f}"
    gp_table["Gravity (m/s^2)"].format = "{:.3f}"
    gp_table["Scale Height (km)"].format = "{:.1f}"
    gp_table["Lambda"].format = "{:.3f}"
    gp_table["Period (days)"].format = "{:.1f}"
    gp_table["Sep (au)"].format = "{:.3f}"
    gp_table["Lstar (L☉)"].format = "{:.3f}"
    gp_table["Mstar (M☉)"].format = "{:.3f}"
    print(gp_table)

    # Return units: Earth
    return planetnum, mp_gp, rp_gp, core_m_gp, teq_gp, wmf_gp, sep_gp, period_gp, H_gp, lambda_gp, startype_gp, mstar_gp, lstar_gp



# calculate planet period
def period(mstar,mp,sep):
    mstar = mstar*mml.M_sun
    mp = mp*mml.M_e
    sep = sep*mml.au
    day_s = 60*60*24    # s/day

    C = (4 * np.pi**2 )/(mml.Ggrav * (mstar + mp)) #constant, m^3/s^2
    T = np.sqrt((sep**3 * C))

    return T / day_s
    

# calculate planet gravity
def grav(mp,rp):
    rp = rp*mml.R_e #m
    mp = mp*mml.M_e #kg

    G = 6.67e-11    # mml.Ggrav
    grav = (mp * G)/(rp**2)  #m/s
    
    return grav

# calculate scale height
def scaleheight(Teq,rp,mp):

    mmolm = 18.01528 * 0.001 # kg/mol, mean molecular mass of the atmosphere (g/mol*0.001kg/g)
    R = 8.31446 #J/Kmol, molar gas constant
    
    g = grav(mp,rp) # m/s, gravity of planet

    # calculate scale height
    # H = (Teq*R)/(gravp*mmolm)  # m
    
    H = (Teq*R)/(g*1e-3) # m
    
    return H
    