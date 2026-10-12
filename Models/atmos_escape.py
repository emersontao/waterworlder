import numpy as np

# fundamental constants

# constants all in cgs units

kb = 1.38065e-16 # boltzmann's constant
mh = 1.673558e-24 # mass of hydrogen atom in grams
sigma = 5.6704e-5 # stefan-Boltzmann constant

G=6.67408e-8 # Gravitational constant

earth_mass_to_g = 5.972e27
earth_radius_to_cm = 6.371e8
rho_earth_cgs = 5.51 # density of Earth
year_to_sec = 31557600. # number of seconds in a year
gyr_to_sec = 1e9*year_to_sec
day_to_sec = 86400. # number of secconds in a day

Tsun = 5778. #K
Rsun = 6.9551e10 # cm
Msun = 1.9885e33
au_to_cm = 1.496e13
Period_Earth = 365.25 #days

#My units in SI
me = 5.972e24
re = 6.371e6

Ggrav=6.67408e-11
sigmasb=5.670374419e-8
gyr=gyr_to_sec
Rmol=8.31446261815324
Nav=6.02214076e23
Teq_e=5772.0*np.sqrt(1.362e9/(4.0*1.5e11))


#from EvapMass.mass_loss import efficiency

get_shape = np.load("../eff_shape_file.npy", allow_pickle=True).item()

def get_scaled_hydro_eff(R_planet,M_planet):

    #return 0.05

    ## convert Rp and Mp from Re and Me to m and kg
    Rp = R_planet * re
    Mp = M_planet * me

    ## Returns efficiency according to Owen & Wu 2017

    vesc = np.sqrt ( 2. * Ggrav * Mp / Rp )

    eff = 0.1 * (vesc / 15e3)**(-2.)

    return eff

    # convert Rp and Mp from Re and Me to cm and grams
    Rp = R_planet * earth_radius_to_cm
    Mp = M_planet * earth_mass_to_g

    ## returns the scaled efficiency factor of the Owen & Jackson (2012) evaporation rates

    # the shape is scaled via the escape velocity to a mass of 7.603262769401823e+27 g

    Mp_scale = 7.603262769401823e+27

    ## Radius scale first 

    scaled_eff = 10.**get_shape(np.log10(Rp*Mp_scale/Mp))

    ## this scaled efficiency returns either a scaled value (scaled to max value in table) for the
    ## efficiency or extropolates the efficiency at a constant value of large planets
    ## that would be undergoing Roche Lobe overflow in the Owen & Jackson (2012) tables
    ## E.g. grey region of Figure 5 in Owen & Jackson (2012)


    ### scale mass

    mass_scale = (1. + (np.sqrt(Mp/1e29))**10.)**(1./10.)

    return scaled_eff * mass_scale

# Jeans escape parameter (no units)
def Lambda(rpe,mpe,teq,mu):
    rp = rpe*re
    mp = mpe*me
    return (Ggrav*mp*mu)/(rp*Rmol*teq)

# mass loss rate by Jeans escape (kg/s)
def mdot_Jeans(rpe,mpe,teq,mu):
    n_exo=1e6
    rp = rpe*re
    mp = mpe*me

    mdot = 2.0*np.sqrt(np.pi)*(rp)**2*mu/Nav*n_exo \
        *np.sqrt(2.0*Ggrav*mp/(rp)) \
            *(1.0+Lambda(rpe,mpe,teq,mu))/np.sqrt(Lambda(rpe,mpe,teq,mu)) \
                *np.exp(-Lambda(rpe,mpe,teq,mu))

    return mdot


# mass loss rate by photoevaporation (kg/s)
def mdot_photoevap(rpe,mpe,teq,ltot,lstar,mstar):
    epsilon = get_scaled_hydro_eff(rpe,mpe)
    rp = rpe*re
    mp = mpe*me
    sep = np.sqrt(lstar/(4.0*4.0*np.pi*sigmasb*teq**4))
    xi = (sep/rp)*(mpe/mstar/3)**(1/3)
    ktide = 1 - 3/(2*xi) + 1/(2*(xi)**3) # Correction for the Roche lobe
    return epsilon/ktide*(ltot)*(rp)**3/(Ggrav*mp*(2.0*sep)**2)

# total mass lost by xuv blowoff over planet lifetime (kg)
def deltam_photoevap(rpe,mpe,teq,exuv,lstar,mstar):
    epsilon = get_scaled_hydro_eff(rpe,mpe)
    rp = rpe*re
    mp = mpe*me
    sep = np.sqrt(lstar/(4.0*4.0*np.pi*sigmasb*teq**4))
    xi = (sep/rp)*(mpe/mstar/3)**(1/3)
    ktide = 1 - 3/(2*xi) + 1/(2*(xi)**3) # Correction for the Roche lobe
    return epsilon/ktide*(exuv)*(rp)**3/(Ggrav*mp*(2.0*sep)**2)

