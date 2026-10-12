# this module creates an interpolator in the grid of mr_all_2024
#it's a version that avoids using dictionaries i think? -- Emerson

import module_mass_loss as mml
import numpy as np 
from scipy.interpolate import RegularGridInterpolator


# Load the data
path_models = "./"

path_aguichine = path_models + "Aguichine2021_mr_all_2024.dat"

listmbfull,listmafull,listrbfull,listrafull,listerrfull = np.loadtxt(path_aguichine,unpack=True,usecols=(4,5,6,7,8))

# Make IOP interpolator

dimcmf = np.linspace(0.0,0.9,10)
dimwmf = np.linspace(0.0,1.0,11)
dimteq = np.linspace(400.0,1300.0,10)
dimmass = np.logspace(np.log10(0.2),np.log10(20.0),20)

data_radius = np.transpose( np.reshape(listrbfull+listrafull,(10,10,11,20)) , (0,2,1,3))
data_errcode = np.transpose( np.reshape(listerrfull,(10,10,11,20)) , (0,2,1,3))

interp_radius_iop = RegularGridInterpolator((dimcmf, dimwmf, dimteq, dimmass), data_radius, method='linear', bounds_error=False, fill_value=None)
interp_errcode_iop = RegularGridInterpolator((dimcmf, dimwmf, dimteq, dimmass), data_errcode, method='linear', bounds_error=False, fill_value=None)

#computes the radius of a planet, using the IOP model from Aguichine+2021
#the planet input should contain:
# - planet mass
# - planet water mass fraction (wmf)
# - planet equilibrium temperature (teq)
def radius_iop(cmf,wmf,teq,mp):
    mp = mp / mml.M_e
    # Earth-like value of the core mass fraction
    #cmf = 0.325

    # radius is computed by linear interpolation in the data grid of the IOP model
    rp = interp_radius_iop([cmf,wmf,teq,mp]).item(0)
    return rp * mml.R_e

#determines if the radius from the IOP model is valid or not
# returns 0 if all data points are valid
# returns >0 if one or more data points are not valid
def error_iop(cmf,wmf,teq,mp):
    mp = mp / mml.M_e
    # Earth-like value of the core mass fraction
    #cmf = 0.325

    # radius is computed by linear interpolation in the data grid of the IOP model
    errcode = interp_errcode_iop([cmf,wmf,teq,mp]).item(0)
    return errcode

# ---------------------------------------------------------------------

#dichotomy/bisection root-finding method
def bisection(func,a,b,precision, max_steps):
    steps = 0
    for _ in range(max_steps):
        c = (a + b) / 2
        steps += 1

        if abs(a - b) < precision:
            print(f'[BISECTION] Steps taken = {steps} \t Interval size = { 1/(2**steps) }')
            return c
        if abs(func(c)) < precision:
            return c
        if func(c) < 0:
            a = c
        else:
            b = c
    
    print('ERROR: Did not converge?')
    print(f'Steps taken = {steps} \t Interval size = { 1/(2**steps) }')
    return c

#inverts the radius_iop equation to solve for wmf

def solve_wmf(cmf,teq,rp,mp,precision=1e-8,max_steps=100):
    a, b = 0, 1 #define interval
    #rp = ss.Neptune['rp']
    mp = mp / mml.M_e
    
    rpmax = radius_iop(cmf,1,teq,mp)
    rpmin = radius_iop(cmf,0,teq,mp)

    if rp > rpmax:
        rp = rpmax
        print(f'RP: {rp}\nRPMAX: {rpmax}\nRadius > 100% Water')
        return 1
    
    if rp < rpmin:
        rp = rpmin
        print(f'Radius < 0% Water')
        return 0

    def objective(wmf):
        return radius_iop(cmf,wmf,teq,mp) - rp   #this equation = 0
    
    #make sure the sign changes between endpoints
    print(f'[SOLVE WMF] wmf(0) = {objective(0)} \t wmf(1) = {objective(1)}')
    if objective(0) * objective(1) > 0:
        raise ValueError('The function does not change sign in the interval [0, 1]')
    
    wmf = bisection(objective,0,1,precision,max_steps)

    #add a check for radius
    calculated_rp = radius_iop(cmf,wmf,teq,mp)
    if rp > calculated_rp :
        rp = calculated_rp        #if the given rp (from ss.py) is GREATER than possible given params
        print('WARNING: Calculated radius override given radius')
        wmf = bisection(objective,0,1,precision,max_steps)      #recalculate wmf

    return wmf