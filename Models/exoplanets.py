# module of exoplanet data

import module_mass_loss as mml
M_e = mml.M_e  # kg (mass of Earth)
R_e = mml.R_e  # m (radius of Earth)
M_sun = mml.M_sun  # kg (mass of the Sun)
R_sun = mml.R_sun  # m (radius of the Sun)
L_sun = mml.L_sun  # or 3.828e26 W (luminosity of the Sun)
T_sun = mml.T_sun  # K (temperature of the Sun)
au = mml.au
Gyr = mml.Gyr


Wolf_327b = {
    'name': 'Wolf 327b',
    'mp': 2.53 * M_e, 
    'rp': 1.24 * R_e,
    'teq': 996,
    'sep': 0.01 * au,
    'period': 0.57347,
    'mstar': 0.405 * M_sun,
    'rstar': 0.406 * R_sun,
    'lstar': L_sun * 10**(-1.630),
    'tstar': 3542,
    #'age': 
}
#https://exoplanetarchive.ipac.caltech.edu/overview/Wolf%20327%20b

#not complete
TOI_1450A_b = {
    'name': 'K2_38c',
    'mp': 1.25800000 * M_e, 
    'rp': 1.13000000 * R_e,
    'teq': 722,
    'sep': 0.0959 * au,
    'period': 2.04392740000,
    'mstar': 1.054 * M_sun,
    'rstar': 1.141 * R_sun,
    'lstar': L_sun * 10**(-0.0374),
    'tstar': 3437,
    #'age': 
}


K2_38c = {
    'name': 'K2 38c',
    'mp': 7.4 * M_e, 
    'rp': 2.49 * R_e,
    'teq': 949,
    'sep': 0.0959 * au,
    'period': 10.56103,
    'mstar': 1.054 * M_sun,
    'rstar': 1.141 * R_sun,
    'lstar': L_sun * 10**(-0.0374),
    'tstar': 5705,
    #'age': 
} #https://exoplanetarchive.ipac.caltech.edu/overview/K2-38%20c

K2_138b = {     #good
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
}   #https://exoplanetarchive.ipac.caltech.edu/overview/k2-138

K2_180b = {
    'name': 'K2_180b',
    'mp': 11.4 * M_e,
    'rp': 2.466 * R_e,
    'teq': 797.5,
    'sep': 0.0756 * au,
    'period': 8.8656, #days
    'mstar': 0.735 * M_sun,
    'rstar': 0.719 * R_sun,
    'lstar': L_sun * 10**(-0.413),
    'tstar': 5365,
    'age': 7.5 * Gyr,
}   #https://exoplanetarchive.ipac.caltech.edu/overview/K2-180%20b

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

GJ_3473b = {
    'name': 'GJ 3473b',
    'mp': 1.86 * M_e,
    'rp': 1.264 * R_e,
    'teq': 773,
    'sep': 0.01589 * au,
    'period': 1.198,
    'mstar': 0.360 * M_sun,
    'rstar': 0.364 * R_sun,
    'lstar': L_sun * 10**(4.81),
    'tstar': 3347,
    #'age': 
}   #https://exoplanetarchive.ipac.caltech.edu/overview/HD%2086226

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

HD_110067b = {              # very similar to Kepler-19b. Pick one to run
    'name': 'HD_110067b',
    'mp': 5.69 * M_e,
    'rp': 2.20 * R_e,
    'teq': 800,
    'sep': 0.07930 * au,
    'period': 9.113678,
    'mstar': 0.798 * M_sun,
    'rstar': 0.788 * R_sun,
    'lstar': L_sun * 10**(-0.3622405),
    'tstar': 5266,
    #'age': 
}   #https://exoplanetarchive.ipac.caltech.edu/overview/HD%20110067%20b




HD_106315b = {
    'name': 'HD_106315b',
    'mp': 12.6 * M_e, 
    'rp': 2.44 * R_e,
    'teq': 1153,
    'sep': 0.0907 * au,
    'period': 9.55237,
    'mstar': 1.079 * M_sun,
    'rstar': 1.31 * R_sun,
    'lstar': L_sun * 10**(0.394543),
    'tstar': 6256,
    #'age': 
} #https://exoplanetarchive.ipac.caltech.edu/overview/HD%20106315%20b

TOI_561b = {            # VERY SMALL WMF not good
    'name': 'TOI_561b',
    'mp': 2.00 * M_e,
    'rp': 1.425 * R_e,
    'teq': 2310,
    'sep': 0.0106 * au,
    'period': 0.44656,
    'mstar': 0.806 * M_sun,
    'rstar': 0.843 * R_sun,
    'lstar': L_sun * 10**(-0.273),
    'tstar': 5372,
    'age': 11.0 * Gyr,
}

TOI_1416b = {           #short evolution, not good 
    'name': 'TOI_1416b',
    'mp': 3.48 * M_e,
    'rp': 1.62 * R_e,
    'teq': 1517,
    'sep': 0.0190 * au,
    'period': 1.069756,
    'mstar': 0.798 * M_sun,
    'rstar': 0.793 * R_sun,
    'lstar': L_sun * 10**(-0.4600775),
    'tstar': 4884,
    #'age': 
} #https://exoplanetarchive.ipac.caltech.edu/overview/TOI-1416

TOI_1836 = {
    'name': 'TOI_1836b',
    'mp': 8.0 * M_e,
    'rp': 2.59 * R_e, 
    'teq': 2007.59,
    'sep': 0.03089 * au,
    'period': 1.7727,
    'mstar': 1.2518 * M_sun,
    'rstar': 1.6454 * R_sun,
    'lstar': L_sun * 10**(0.5698),
    'tstar': 6237,
}

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
} #https://exoplanetarchive.ipac.caltech.edu/overview/Kepler-19%20b

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

Kepler_68b = {
    'name': 'Kepler_68b',
    'mp': 8.03 * M_e, 
    'rp': 2.357 * R_e,
    'teq': 1275,
    'sep': 0.06135 * au,
    'period': 5.39875259,
    'mstar': 1.0420 * M_sun,
    'rstar': 1.2564 * R_sun,
    'lstar': L_sun * 10**(0.23106034),
    'tstar': 5847,
    #'age': 
} #https://exoplanetarchive.ipac.caltech.edu/overview/Kepler-68%20b


TOI_544b = {
    'name': 'TOI 544b',
    'mp': 2.89 * M_e, 
    'rp': 2.018 * R_e,
    'teq': 999,
    'sep': 0.0225 * au,
    'period': 1.54,
    'mstar': 0.631 * M_sun,
    'rstar': 0.623 * R_sun,
    'lstar': L_sun * 10**(-0.9004854),
    'tstar': 4169,
    #'age': 
} #https://exoplanetarchive.ipac.caltech.edu/overview/TOI-544%20b#planet_TOI-544-b_collapsible

Kepler_109b = {
    "name": 'Kepler 109b',
    "mp": 4.9 * M_e,
    "rp": 2.49 * R_e,
    "teq": 1276,
    "sep": 0.0701 * au,
    "period": 6.4816307,
    "mstar": 1.094 * M_sun,
    "rstar": 1.387 * R_sun,
    "lstar": L_sun * 10**(0.3465),
    "tstar": 5950,
    "age": 6.2 * Gyr,
}

KOI_94b = {
    'name': 'KOI_94b',
    'mp': 10.5 * M_e,
    'rp': 1.71 * R_e,
    'teq': 1486,
    'sep': 0.05119 * au,
    'period': 3.7432,
    'mstar': 1.277 * M_sun,
    'rstar': 1.52 * R_sun,
    'lstar': L_sun * 10**(4.181),
    'tstar': 6182,
    'age': 3.16 * Gyr,
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
} #https://exoplanetarchive.ipac.caltech.edu/overview/EPIC%20229004835%20b

L_98_59b = {
    'name': 'L_98_59b',
    'mp': 0.40 * M_e,
    'rp': 0.850 * R_e,
    'teq': 627,
    'sep': 0.02191 * au,
    'period': 2.25311,
    'mstar': 0.273 * M_sun,
    'rstar': 0.303 * R_sun,
    'lstar': L_sun * 10**(-1.948),
    'tstar': 3415,
}


# Kepler_87c = {
#     'name': 'Kepler_87c',
#     'mp':  * M_e,
#     'rp': * R_e,
#     'teq': ,
#     'sep':  * au,
#     'period': ,
#     'mstar':  * M_sun,
#     'rstar': * R_sun,
#     'lstar': L_sun * 10**(),
#     'tstar':  ,
# }

# --- ---

Neptune = {
    "name" : "Neptune",
    "mp": 17.2 * M_e,
    "rp": 3.88 * R_e,
    "teq": 50, #K
    "sep": 30.0611 * mml.au * 0.2, # m
    "period": 59800,
    "mstar": M_sun, # kg (mass of the Sun)
    "rstar": R_sun, #m (radius of the star)
    "lstar": 3.828e26, #W (luminosity of the sun)
    "tstar": T_sun,       # Sun temperature
    "age": 4.571 * mml.Gyr # Planet/star age in s
}

# Sample = {
#     'name': 'Sample',
#     'mp': * M_e,
#     'rp': * R_e,
#     'teq': ,
#     'sep':  * au,
#     'period': ,
#     'mstar':  * M_sun,
#     'rstar': * R_sun,
#     'lstar': L_sun * 10**(),
#     'tstar':  ,
# }