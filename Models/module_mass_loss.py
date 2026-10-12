# UNITS: MKS
import numpy as np

#physical constants
Ggrav = 6.67e-11 #gravitational constant m^3/kg/sec^2
Rmol = 8.314 #ideal gas constant (J/mol/K)

#astronomical constants
M_e = 5.97237e24 #Earth mass (kg)
R_e = 6.371e6 #Earth radius (m)
M_sun = 1.9885e30 #Sun mass (kg)
R_sun = 6.9551e10 #Sun radius (m)
L_sun = 3.828e26 #Total luminosity of the Sun (Watts = Joules/sec)
T_sun = 5778. #Temperature of the Sun (K)
au = 1.496e11 #Astronomical Unit (m)
Gyr = 3.15e16 #Seconds in 1 Gyr

#other constants
L_he = 4.5e-6*L_sun + 4.6e-7*L_sun #EUV + X-Ray luminosities of the Sun (W)
eta = 0.1 #photoevaporation efficiency constant
sigmasb = 5.670374419e-8 # Stefan-Boltzmann constant  (W.m-2.K-4)


#efficiency of atmospheric mass loss by photoevaporation
#from Owen & Wu 2017
def efficiency(mp,rp,mstar,sep):
    #escape velocity at the surface of the planet (m/s)
    vesc = np.sqrt(2 * Ggrav * mp / rp)

    #mass loss efficiency
    eff = 0.1 * (vesc / 15e3)**(-2.)

    #Correction for the Roche lobe
    #from Salz+2016
    xi = (sep/rp)*(mp/mstar/3)**(1/3)
    ktide = 1 - 3/(2*xi) + 1/(2*(xi)**3) 
    
    #total efficiency
    efficiency = eff/ktide
    return efficiency

#high-energy luminosity (X + UV) emitted by the star
#from Sanz-Forçada+2011
def L_xuv(age,lstar):

    # duration of the saturation regime in seconds
    tsat=5.72e15*lstar**(-0.65)

    # lxuv in non-saturation regime
    lx = 1.89e21*(age/Gyr)**(-1.55)
    leuv = 10**(3.82)*lx**0.86
    lxuv = lx + leuv

    # lxuv in saturation regime
    lxsat = 6.3e-4*lstar
    leuvsat = 10**(3.82)*lxsat**0.86
    lxuvsat = lxsat + leuvsat

    new_lxuv = np.minimum(lxuv,lxuvsat)
    return new_lxuv

    
#mass-loss rate mdot in kg/sec
def mdot(mp,rp,mstar,sep,age,lstar):
    Lhe = L_xuv(age,lstar)      #luminosity of high-energy photons
    eta = efficiency(mp,rp,mstar,sep)
    mdot = eta * Lhe * rp**3 / (4.0 * Ggrav * mp * sep**2)
    return mdot


#get separation (planet-star distance) from equilibrium temperature
def get_sep_from_teq(lstar,teq):
    sep = np.sqrt(lstar / (16 * np.pi * sigmasb * teq**4))
    return sep

#get equilibrium temperature from separation (planet-star distance)
def get_teq_from_sep(lstar,sep):
    teq = (lstar / (16 * np.pi * sigmasb * sep**2))**0.25
    return teq
