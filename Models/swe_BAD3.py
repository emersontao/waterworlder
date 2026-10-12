# this module creates an interpolator in the grid of mr_all_2024
# I tried to fix this one and it doesn't matter. It's fucked.
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

path_swe = path_models + "A24_SWE_all_v3_fixed.dat"

labels_host = ["Type M", "Type G"]
teqs = [400, 500, 700, 900, 1100, 1300, 1500]
# teqs = [500, 600, 700]
wmfs = [0.1,1,10,20,30,40,50,60,70,80,90,100]  # 12 water mass fractions from 0.05 to 0.60
masses = [0.2       ,  0.254855  ,  0.32475535,  0.41382762,  0.52733018,
        0.67196366,  0.85626648,  1.09111896,  1.39038559,  1.77173358,
        2.25767578,  2.87689978,  3.66596142,  4.67144294,  5.95270288,
        7.58538038,  9.66586048, 12.31696422, 15.69519941, 20.        ]  # 20 points in mass from 0.1 to 2.0

ages = np.array([0.001,0.0015,0.002,0.003,0.005,0.01,
                        0.02,0.03,0.05,
                        0.1,0.2,0.5,
                        1.0,2.0,5.0,
                        10,20])

listrpfull = np.loadtxt(path_swe,skiprows=36,unpack=True,usecols=(5))
listrpfull = np.delete(listrpfull, np.arange(17, listrpfull.size, 18))

mask = listrpfull == -1.0
listrpfull[mask] == np.inf

listrpfull_m = listrpfull[0:int(len(listrpfull)/2)]
listrpfull_g = listrpfull[int(len(listrpfull)/2):]


# Make SWE interpolator

dim_wmf = np.array(wmfs)/100
dim_teq = np.array(teqs)
dim_mass = np.array(masses)
dim_age = ages

data_radius_m = np.reshape(listrpfull_m,(12,7,20,17))
data_radius_g = np.reshape(listrpfull_g,(12,7,20,17))

fill_value = None
interp_m = RegularGridInterpolator((dim_wmf, dim_teq, dim_mass, dim_age), data_radius_m, method='slinear', bounds_error=False, fill_value=fill_value)
interp_g = RegularGridInterpolator((dim_wmf, dim_teq, dim_mass, dim_age), data_radius_g, method='slinear', bounds_error=False, fill_value=fill_value)


# INTERPOLATION RADIUS WRAPPERS
# Dependent on star luminosity
def radius_m(wmf,teq,mp,age):
    #M-type
    mp = mp / M_e
    age = age / Gyr
    return interp_m([wmf,teq,mp,age]).item(0) * R_e #m

def radius_g(wmf,teq,mp,age):
    # G-type
    mp = mp / M_e
    age = age / Gyr
    return interp_g([wmf,teq,mp,age]).item(0)  * R_e #m


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
    
    #print(f'[BISECTION: FAILED]\nDid not converge\nSteps taken = {steps} \nInterval size = { 1/(2**steps) }')
    return c


# WMF SOLVER
# invert the radius_iop equation to solve for wmf
def solve_wmf_m(rp,teq,mp,age,precision=1e-8,max_steps=100):
    #M-type stars
    a, b = 0, 1 #define interval
    rp = rp / R_e      # in Earth radii
    mp = mp / M_e      # in Earth masses
    age = age / Gyr    # in Gyr
        
    rpmax = radius_m(1,teq,mp,age)
    rpmin = radius_m(0,teq,mp,age)

    if rp > rpmax:
        rp = rpmax
        print('RADIUS:',rp,'R_e \tRP(wmf=1.0):',rpmax,'R_e \nRadius > 100% Water')
        return 1
    
    if rp < rpmin:
        rp = rpmin
        print('RADIUS:',rp,'R_e \tRP(wmf=1.0):',rpmax,'R_e \nRadius < 0% Water')
        return 0

    def objective(wmf):
        return radius_m(wmf,teq,mp,age) - rp
    
    #make sure the sign changes between endpoints
    wmf_0 = objective(0)
    wmf_1 = objective(1)
    #print(f'[SOLVE WMF]  WMF(0.0) = {wmf_0}\tWMF(1.0) = {wmf_1}')
    if objective(0) * objective(1) > 0:
        raise ValueError('The function does not change sign in the interval [',a,',',b,']')
    
    wmf = bisection(objective,0,1,precision,max_steps)
    
    # #IF USING EXISTING EXOPLANET:
    # # check for radius
    # calculated_rp = interp_radius_m([wmf,teq,mp,age])
    # if rp > calculated_rp :
    #     rp = calculated_rp
    #     print('WARNING: Calculated radius override given radius')
    #     wmf = bisection(objective,0,1,precision,max_steps)      #recalculate wmf

    #print(f'WMF: {wmf:.4f}')
    return wmf


def solve_wmf_g(rp,teq,mp,age,precision=1e-8,max_steps=100):
    # G-type stars
    a, b = 0, 1 #define interval
    # rp = rp / R_e      # in Earth radii
    # mp = mp / M_e      # in Earth masses
    # age = age / Gyr    # in Gyr
        
    rpmax = radius_g(0,teq,mp,age)
    rpmin = radius_g(1,teq,mp,age)

    if rp > rpmax:
        rp = rpmax
        print('RADIUS:',rp,'R_e \tRP(wmf=1.0):',rpmax,'R_e \nRadius > 100% Water')
        return 1
    
    if rp < rpmin:
        rp = rpmin
        print('RADIUS:',rp,'R_e \tRP(wmf=1.0):',rpmax,'R_e \nRadius < 0% Water')
        return 0

    def objective(wmf):
        return radius_g(wmf,teq,mp,age) - rp
    
    #make sure the sign changes between endpoints
    wmf_0 = objective(0)
    wmf_1 = objective(1)
    #print(f'[SOLVE WMF]  WMF(0.0) = {wmf_0}\tWMF(1.0) = {wmf_1}')
    if objective(0) * objective(1) > 0:
        raise ValueError('The function does not change sign in the interval [0, 1]')
    
    wmf = bisection(objective,0,1,precision,max_steps)
    
    # #IF USING EXISTING EXOPLANET: 
    # # check for radius
    # calculated_rp = interp_radius_m([wmf,teq,mp,age])
    # if rp > calculated_rp :
    #     rp = calculated_rp
    #     print('WARNING: Calculated radius override given radius')
    #     wmf = bisection(objective,0,1,precision,max_steps)      #recalculate wmf

    #print(f'WMF: {wmf:.4f}')
    return wmf
