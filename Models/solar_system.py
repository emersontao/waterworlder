import module_mass_loss as mml

M_e = mml.M_e  # kg (mass of Earth)
R_e = mml.R_e  # m (radius of Earth)
M_sun = mml.M_sun  # kg (mass of the Sun)
R_sun = mml.R_sun  # m (radius of the Sun)
L_sun = mml.L_sun  # W (luminosity of the Sun)
#L_sun = 3.828e26    #Total luminosity of the Sun (Watts = Joules/sec)
T_sun = mml.T_sun  # K (temperature of the Sun)

Mercury = {
    "name": "Mercury",
    "mp": 0.0553 * M_e,  # Earth masses
    "rp": 0.3829 * R_e,   # Earth radii
    "teq": 340,       # K
    "sep": 0.3871 * mml.au, # m
    "period": 87.97,      # days
    "mstar": M_sun,       # Sun masses
    "rstar": R_sun,       # Sun radii
    "lstar": L_sun,       # Sun luminosity
    "tstar": T_sun,       # Sun temperature
    "age": 4.571 * mml.Gyr # Planet/star age in s
}

Venus = {
    "name" : "Venus",
    "mp": 0.815 * M_e,
    "rp": 0.9499 * R_e,
    "teq": 735,
    "sep": 0.7233 * mml.au, # m
    "period": 224.7,
    "mstar": M_sun,
    "rstar": R_sun,
    "lstar": L_sun,
    "tstar": T_sun,       # Sun temperature
    "age": 4.571 * mml.Gyr # Planet/star age in s
}

Earth = {
    "name" : "Earth",
    "mp": M_e,
    "rp": R_e,
    "teq": 255,
    "sep": 1.0 * mml.au, # m
    "period": 365.25,
    "mstar": M_sun,
    "rstar": R_sun,
    "lstar": L_sun,
    "tstar": T_sun,       # Sun temperature
    "age": 4.571 * mml.Gyr # Planet/star age in s
}

Mars = {
    "name" : "Mars",
    "mp": 0.107 * M_e,
    "rp": 0.532 * R_e,
    "teq": 215,
    "sep": 1.5237 * mml.au, # m
    "period": 687,
    "mstar": M_sun,
    "rstar": R_sun,
    "lstar": L_sun,
    "tstar": T_sun,        # Sun temperature
    "age": 4.571 * mml.Gyr # Planet/star age in s
}

Jupiter = {
    "name" : "Jupiter",
    "mp": 317.8 * M_e,
    "rp": 11.21 * R_e,
    "teq": 110,
    "sep": 5.2028 * mml.au, # m
    "period": 4331,
    "mstar": M_sun,
    "rstar": R_sun,
    "lstar": L_sun,
    "tstar": T_sun,       # Sun temperature
    "age": 4.571 * mml.Gyr # Planet/star age in s
}

Saturn = {
    "name" : "Saturn",
    "mp": 95.2 * M_e,
    "rp": 9.45 * R_e,
    "teq": 85,
    "sep": 9.5388 * mml.au, # m
    "period": 10747,
    "mstar": M_sun,
    "rstar": R_sun,
    "lstar": L_sun,
    "tstar": T_sun,       # Sun temperature
    "age": 4.571 * mml.Gyr # Planet/star age in s
}

Uranus = {
    "name" : "Uranus",
    "mp": 14.6 * M_e,
    "rp": 4.01 * R_e,
    "teq": 60,
    "sep": 19.1914 * mml.au, # m
    "period": 30589,
    "mstar": M_sun,
    "rstar": R_sun,
    "lstar": L_sun,
    "tstar": T_sun,       # Sun temperature
    "age": 4.571 * mml.Gyr # Planet/star age in s
}

Neptune = {
    "name" : "Neptune",
    "mp": 17.2 * M_e,
    "rp": 3.88 * R_e,
    "teq": 50,
    "sep": 30.0611 * mml.au, # m
    "period": 59800,
    "mstar": M_sun,
    "rstar": R_sun,
    "lstar": L_sun,
    "tstar": T_sun,       # Sun temperature
    "age": 4.571 * mml.Gyr # Planet/star age in s
}
