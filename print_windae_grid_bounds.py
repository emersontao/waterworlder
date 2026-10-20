# test.py
import numpy as np

Lsun = 3.828e26 #Total luminosity of the Sun (Watts = Joules/sec)
sigma_sb = 5.670374419e-8 # Stefan-Boltzmann constant  (W.m-2.K-4)
Gyr = 3.15e16 #Seconds in 1 Gyr

#xuv flux at the planet's orbit
def find_Fxuv(Lxuv,a):
    '''
    Lxuv:   XUV luminosity [W]
    a:      semi-major axis [m]
    '''
    Fxuv =  Lxuv*1000 / (4*np.pi*a**2) 
    #take in SI, convert to cgs
    return Fxuv

def find_Lxuv(age,lstar):
    #high-energy luminosity (X + UV) emitted by the star
    #Sanz-Forçada+2011

    #L_xuv in non-saturation regime
    Lx = 1.89e21 * (age/Gyr)**(-1.55)
    Leuv = 10**(3.82) * Lx**0.86
    Lxuv = Lx + Leuv

    #L_xuv in saturation regime
    Lxsat = lstar * 6.3e-4
    Leuvsat = 10**(3.82) * Lxsat**0.86
    Lxuvsat = Lxsat + Leuvsat

    new_lxuv = np.minimum(Lxuv,Lxuvsat) #final Lxuv
    return new_lxuv

def log_Fxuv_axis(mstar_axis,teq_axis,ct):
    '''
    mstar_axis: array-like shape (N,-1) of stellar masses [solar mass]
    teq_axis:   array-like shape (M,-1) of planet equilibrium temperatures [K]
    ct:         int, number of points in Fxuv axis
    -----------
    returns:    array-like shape (ct,-1) of Fxuv values [erg/s
    '''
    #figure out fxuv constraints
    lowest_lstar = 0.23*Lsun*(np.min(mstar_axis))**2.3
    lowest_Lxuv = find_Lxuv(20*Gyr, lowest_lstar)
    lowest_semimajor = np.sqrt(lowest_lstar / (16*np.pi* sigma_sb * np.min(teq_axis)**4))
    lowest_Fxuv = find_Fxuv(lowest_Lxuv,lowest_semimajor)

    highest_lstar = 0.23*Lsun*(np.max(mstar_axis))**2.3
    highest_Lxuv = find_Lxuv(1e-2*Gyr, highest_lstar)
    highest_semimajor = np.sqrt(highest_lstar / (16*np.pi* sigma_sb * np.max(teq_axis)**4))
    highest_Fxuv = find_Fxuv(highest_Lxuv, highest_semimajor)

    Fxuv_axis = np.array(np.logspace(np.log10(lowest_Fxuv), np.log10(highest_Fxuv), ct))
    return Fxuv_axis

# -------------- INITIALIZE THE GRID -------------- #
# mp_axis     = np.array([0.2, 2, 4, 7, 11, 15, 18, 20, 22.5, 25, 27, 30]) #earth mass
mp_axis     = np.array([2, 4, 7, 11, 15, 18, 20, 22.5, 25, 27, 30]) #earth mass
grav_axis   = np.array([515, 600, 724, 1000, 1400, 2000, 3000]) #cm/s^2
teq_axis    = np.array([1500, 1225, 950, 675, 400]) #K
mstar_axis  = np.array([0.1, 0.22, 0.3, 0.55, 1]) #sun mass
Fxuv_axis   = np.array(list(reversed(log_Fxuv_axis(mstar_axis,teq_axis,5)))) #erg/s/cm^2


print("mp_axis: ", mp_axis[0],  mp_axis[-1], len(mp_axis))
print("grav_axis: ", grav_axis[0], grav_axis[-1], len(grav_axis))
print("teq_axis: ", teq_axis[0], teq_axis[-1],   len(teq_axis))
print("mstar_axis: ", mstar_axis[0], mstar_axis[-1], len(mstar_axis))
print("Fxuv_axis: ", Fxuv_axis[0],Fxuv_axis[-1], len(Fxuv_axis)) #FuvEarth = 14/3          # Earth's recieved uv [10—91.2 nm] flux (erg/s/cm**2)

print(Fxuv_axis*3/14)