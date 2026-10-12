# exoplanets for Kama. solid planets
import module_mass_loss as mml
M_e = mml.M_e  # kg (mass of Earth)
R_e = mml.R_e  # m (radius of Earth)
M_sun = mml.M_sun  # kg (mass of the Sun)
R_sun = mml.R_sun  # m (radius of the Sun)
L_sun = mml.L_sun  # or 3.828e26 W (luminosity of the Sun)
T_sun = mml.T_sun  # K (temperature of the Sun)
au = mml.au
Gyr = mml.Gyr

K2_138b = {
    'name': 'K2_138b',
    'mp': 3.10 * M_e,
    'rp': 1.51 * R_e,
    'teq': 1157, #K (temperature of the planet)
    'sep': 0.03385 * au,
    'period': 2.35309, #days
    'mstar': 0.9160 * M_sun, #kg (mass of the star)
    'rstar': 0.839 * R_sun, #m (radius of the star)
    'lstar': L_sun * 10**(-0.2877621), #W (luminosity of the star)
    'tstar': 5283, #K (temperature of the star)
    #'age': 0,
}

HD_23472f = {   #GOOD
    'name': 'HD_23472f',
    'mp': 0.770 * M_e,
    'rp': 1.137 * R_e,
    'teq': 630,
    'sep': 0.0906 * au,
    'period': 12.1621,
    'mstar': 0.67 * M_sun,
    'rstar': 0.71 * R_sun,
    'lstar': L_sun * 10**(-0.625),
    'tstar': 4684,
    #'age': 
}

HD_86226c = {
    'name': 'HD_86226c',
    'mp': 7.25 * M_e,
    'rp': 2.16 * R_e,
    'teq': 1311,
    'sep': 0.049 * au,
    'period': 3.98442,
    'mstar': 1.019 * M_sun,
    'rstar': 1.053 * R_sun,
    'lstar': L_sun * 10**(0.07188),
    'tstar': 5863,
    #'age': 
}   #https://exoplanetarchive.ipac.caltech.edu/overview/HD%2086226

Kepler_19b = {              # very similar to HD_110067b. Pick one to run
    'name': 'Kepler_19b',
    'mp': 6.1 * M_e, 
    'rp': 2.209 * R_e,
    'teq': 851,
    'sep': 0.0846 * au,
    'period': 9.286990,
    'mstar': 0.9670 * M_sun,
    'rstar': 0.859 * R_sun,
    'lstar': L_sun * 10**(-0.1605629),
    'tstar': 5541,
    #'age': 
} 

Kepler_51b = {  #GOOD
    'name': 'Kepler_51b',
    'mp': 2.1 * M_e, 
    'rp': 7.1 * R_e,
    'teq': 543,
    'sep': 0.2514 * au,
    'period': 45.1540,
    'mstar': 1.04 * M_sun,
    'rstar': 0.940 * R_sun,
    'lstar': L_sun * 10**(-0.1802617),
    'tstar': 6018,
    'age': 0.3 * Gyr
}

EPIC_229004835b = {
    'name': 'EPIC_229004835b',
    'mp': 10.4 * M_e, 
    'rp': 2.332 * R_e,
    'teq': 804,
    'sep': 0.1237 * au,
    'period': 16.141132,
    'mstar': 0.970 * M_sun,
    'rstar': 0.999 * R_sun,
    'lstar': L_sun * 10**(0.0298928),
    'tstar': 5868,
    #'age': 
}