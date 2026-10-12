# this module creates an interpolator in the grid of mr_all_2024
# Last updated: 04/20/25 
import numpy as np 
from scipy.interpolate import RegularGridInterpolator
import iop_fint as iop
import module_mass_loss as mml

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
#L_he = 4.5e-6*L_sun + 4.6e-7*L_sun #EUV + X-Ray luminosities of the Sun (W)
eta = 0.1 #photoevaporation efficiency constant
sigmasb = 5.670374419e-8 # Stefan-Boltzmann constant  (W.m-2.K-4)

# Load the data
path_models = "./Models/"

path_swe = path_models + "A25_SWE_all.dat"

swe_top = np.array([0,1])
swe_labels_top = ["20 mbar", "1 µbar"]
swe_labels_host = ["Type M", "Type G"]
swe_teqs = [400, 500, 700, 900, 1100, 1300, 1500]
# swe_teqs = [500, 600, 700]
swe_wmfs = [0.1,1,10,20,30,40,50,60,70,80,90,100]  # 12 water mass fractions from 0.05 to 0.60
swe_masses = [0.2       ,  0.254855  ,  0.32475535,  0.41382762,  0.52733018,
        0.67196366,  0.85626648,  1.09111896,  1.39038559,  1.77173358,
        2.25767578,  2.87689978,  3.66596142,  4.67144294,  5.95270288,
        7.58538038,  9.66586048, 12.31696422, 15.69519941, 20.        ]  # 20 points in mass from 0.1 to 2.0

swe_ages = np.array([0.001,0.0015,0.002,0.003,0.005,0.01,
                        0.02,0.03,0.05,
                        0.1,0.2,0.5,
                        1.0,2.0,5.0,
                        10,20])

listrpfull1 = np.loadtxt(path_swe,skiprows=36,unpack=True,usecols=(6))
listrpfull2 = np.loadtxt(path_swe,skiprows=36,unpack=True,usecols=(7))
listrpfull = np.array((listrpfull1,listrpfull2))
listrpfull = np.delete(listrpfull, np.arange(17, listrpfull.size, 18))


# Make SWE interpolator
swe_dim_wmf = np.array(swe_wmfs)/100
swe_dim_teq = np.array(swe_teqs)
swe_dim_mass = np.array(swe_masses)
swe_dim_age = swe_ages
swe_dim_star = np.array([0,1])
swe_dim_top = np.array([0,1])

swe_data_radius = np.reshape(listrpfull,(2,2,12,3,20,17))

fill_value = np.nan
interp_swe = RegularGridInterpolator((swe_dim_top,
                                      swe_dim_star,
                                      swe_dim_wmf,
                                      swe_dim_teq,
                                      swe_dim_mass,
                                      swe_dim_age), 
                                     swe_data_radius, method='slinear', bounds_error=False, fill_value=fill_value)


#ADD WHAT 'TOP' IS
def radius_swe_m(top,wmf,teq,mp,age):
    mp = mp / M_e
    age = age / Gyr
    return interp_swe([top,0,wmf,teq,mp,age]).item(0) * R_e #m, for M dwarf

def radius_swe_g(top,wmf,teq,mp,age):
    mp = mp / M_e
    age = age / Gyr
    return interp_swe([top,1,wmf,teq,mp,age]).item(0)  * R_e #m, for G-star

'''
host_stars = ['M', 'G']
equilibrium_temperatures = [500, 600, 700]
water_mass_fractions = [0.1,1,10,20,30,40,50,60,70,80,90,100]  # 12 water mass fractions from 0.05 to 0.60
planet_masses = [0.2       ,  0.254855  ,  0.32475535,  0.41382762,  0.52733018,
        0.67196366,  0.85626648,  1.09111896,  1.39038559,  1.77173358,
        2.25767578,  2.87689978,  3.66596142,  4.67144294,  5.95270288,
        7.58538038,  9.66586048, 12.31696422, 15.69519941, 20.        ]  # 20 points in mass from 0.1 to 2.0

age_ref = np.array([0.01,0.0101,0.0105,0.011,0.012,0.015,
                        0.02,0.03,0.05,
                        0.1,0.2,0.5,
                        1.0,2.0,5.0,
                        10,20])

listrpfull = np.loadtxt(path_swe,unpack=True,usecols=(5))
listrpfull = np.delete(listrpfull, np.arange(17, listrpfull.size, 18))

listrpfull_m = listrpfull[0:int(len(listrpfull)/2)]
listrpfull_g = listrpfull[int(len(listrpfull)/2):]


# Make IOP interpolator

dimwmf = np.array(water_mass_fractions)/100
dimteq = np.array(equilibrium_temperatures)
dimmass = np.array(planet_masses)
dimage = age_ref

data_radius_m = np.reshape(listrpfull_m,(12,3,20,17))
data_radius_g = np.reshape(listrpfull_g,(12,3,20,17))

interp_radius_swe_m = RegularGridInterpolator((dimwmf, dimteq, dimmass, dimage), data_radius_m, method='linear', bounds_error=False, fill_value=None)
interp_radius_swe_g = RegularGridInterpolator((dimwmf, dimteq, dimmass, dimage), data_radius_g, method='linear', bounds_error=False, fill_value=None)


def radius_swe_m(wmf,teq,mp,age):
    mp = mp / M_e
    age = age / Gyr
    return interp_radius_swe_m([wmf,teq,mp,age]).item(0) * R_e #m, for M dwarf

def radius_swe_g(wmf,teq,mp,age):
    mp = mp / M_e
    age = age / Gyr
    return interp_radius_swe_g([wmf,teq,mp,age]).item(0)  * R_e #m, for G-star

# print(interp_radius_swe_m([0.001,500,0.2,0.01]))
# print(interp_radius_swe_m([0.001,500,0.2,0.02]))
# print(interp_radius_swe_m([1,700,20,20]))

'''


#dichotomy/bisection root-finding method
def bisection(func,a,b,precision, max_steps):
    steps = 0
    for _ in range(max_steps):
        c = (a + b) / 2
        steps += 1

        if abs(a - b) < precision:
            #print(f'[BISECTION: SUCCESS]\nSteps taken = {steps} \nInterval size = { 1/(2**steps) }')
            return c
        if abs(func(c)) < precision:
            return c
        if func(c) < 0:
            a = c
        else:
            b = c
    
    #print('ERROR: Did not converge?')
    #print(f'[BISECTION: FAILED]\nDid not converge\nSteps taken = {steps} \nInterval size = { 1/(2**steps) }')
    return c

#invert the radius_iop equation to solve for wmf
def solve_wmf_m(rp,teq,mp,age,precision=1e-8,max_steps=100):
    a, b = 0, 1 #define interval
    rp = rp / R_e      # in Earth radii
    mp = mp / M_e      # in Earth masses
    age = age / Gyr    # in Gyr
        
    #rpmax = radius_swe(1,teq,mp,age)
    rpmax = interp_radius_swe_m([1,teq,mp,age])
    rpmin = interp_radius_swe_m([0,teq,mp,age])
    #rpmin = radius_swe(0,teq,mp,age)

    if rp > rpmax:
        rp = rpmax
        print('RADIUS:',rp,'R_e \tRP(wmf=1.0):',rpmax,'R_e \nRadius > 100% Water')
        return 1
    
    if rp < rpmin:
        rp = rpmin
        print('RADIUS:',rp,'R_e \tRP(wmf=1.0):',rpmax,'R_e \nRadius < 0% Water')
        return 0

    def objective(wmf):
        return interp_radius_swe_m([wmf,teq,mp,age]) - rp
    
    #make sure the sign changes between endpoints
    wmf_0 = objective(0)
    wmf_1 = objective(1)
    #print(f'[SOLVE WMF]  WMF(0.0) = {wmf_0}\tWMF(1.0) = {wmf_1}')
    if objective(0) * objective(1) > 0:
        raise ValueError('The function does not change sign in the interval [0, 1]')
    
    wmf = bisection(objective,0,1,precision,max_steps)
    
    # #IF USING EXISTING EXOPLANET: 
    # # check for radius
    # calculated_rp = interp_radius_swe_m([wmf,teq,mp,age])
    # if rp > calculated_rp :
    #     rp = calculated_rp
    #     print('WARNING: Calculated radius override given radius')
    #     wmf = bisection(objective,0,1,precision,max_steps)      #recalculate wmf

    #print(f'WMF: {wmf:.4f}')
    return wmf
