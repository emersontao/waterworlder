#wind-ae initialization module
# ... but with weights!
import numpy as np
import os, sys
from scipy.spatial import cKDTree
from itertools import zip_longest, product
from wind_ae.wrapper.wrapper_utils import constants as const

#Constants
sigma_sb = 5.67037e-5 #erg⋅cm−2⋅s−1⋅K−4
G = 6.67e-8         #whatever that is in cgs
Gyr = 3.153e16      #s

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

#xuv flux at the planet's orbit
def find_Fxuv(Lxuv,a):
    '''
    Lxuv:   XUV luminosity [W]
    a:      semi-major axis [m]
    '''
    Fxuv =  Lxuv*1000 / (4*np.pi*a**2) 
    #take in SI, convert to cgs
    return Fxuv

def log_Fxuv_axis(mstar_axis,teq_axis,ct):
    '''
    mstar_axis: array-like shape (N,-1) of stellar masses [solar mass]
    teq_axis:   array-like shape (M,-1) of planet equilibrium temperatures [K]
    ct:         int, number of points in Fxuv axis
    -----------
    returns:    array-like shape (ct,-1) of Fxuv values [erg/s
    '''
    #figure out fxuv constraints
    lowest_lstar = 0.23*const.Lsun*(np.min(mstar_axis))**2.3
    lowest_Lxuv = find_Lxuv(20*Gyr, lowest_lstar)
    lowest_semimajor = np.sqrt(lowest_lstar / (16*np.pi* sigma_sb * np.min(teq_axis)**4))
    lowest_Fxuv = find_Fxuv(lowest_Lxuv,lowest_semimajor)

    highest_lstar = 0.23*const.Lsun*(np.max(mstar_axis))**2.3
    highest_Lxuv = find_Lxuv(1e-2*Gyr, highest_lstar)
    highest_semimajor = np.sqrt(highest_lstar / (16*np.pi* sigma_sb * np.max(teq_axis)**4))
    highest_Fxuv = find_Fxuv(highest_Lxuv, highest_semimajor)

    Fxuv_axis = np.array(np.logspace(np.log10(lowest_Fxuv), np.log10(highest_Fxuv), ct))
    return Fxuv_axis

def normalize(x,input_grid):
    '''
    x: array-like shape (4,) of given coordinate
    input_grid: list of mp, g, sep, Fxuv arrays ([],[],[],[])
    '''
    mp_ax, g_ax, sep_ext, Fxuv_ax = input_grid
    #where x is a coordinate for parameters
    #logspace or linear, use log for gravity eventually?

    x_norm = np.array([(x[0]-mp_ax[0])/(mp_ax[-1]-mp_ax[0]),
                       (x[1]-g_ax[0])/(g_ax[-1]-g_ax[0]),
                       (np.log10(x[2]) - np.log10(sep_ext[0])) / (np.log10(sep_ext[-1]) - np.log10(sep_ext[0])),
                       (np.log10(x[3]) - np.log10(Fxuv_ax[0])) / (np.log10(Fxuv_ax[-1]) - np.log10(Fxuv_ax[0]))
                       ])
    return x_norm

# def normalize(x,input_grid):
#     '''
#     x: array-like shape (5,) of given coordinate
#     input_grid: list of mp, g, teq, mstar, Fxuv arrays ([],[],[],[],[])
#     '''
#     mp_ax, g_ax, teq_ax, ms_ax, Fxuv_ax = input_grid
#     #where x is a coordinate for parameters
#     #logspace or linear, use log for gravity eventually?

#     x_norm = np.array([(x[0]-mp_ax[0])/(mp_ax[-1]-mp_ax[0]),
#                        (x[1]-g_ax[0])/(g_ax[-1]-g_ax[0]),
#                        (x[2]-teq_ax[0])/(teq_ax[-1]-teq_ax[0]),
#                        (np.log10(x[3]) - np.log10(ms_ax[0])) / (np.log10(ms_ax[-1]) - np.log10(ms_ax[0])),
#                        (np.log10(x[4]) - np.log10(Fxuv_ax[0])) / (np.log10(Fxuv_ax[-1]) - np.log10(Fxuv_ax[0]))
#                        ])
#     return x_norm

#find the nearest point to the given coordinate in a cloud
def nearest_in_cloud(coords, x, weights=[0.5, 5, 2, 1]):
    """
    coords: array shape (N,4)
    x:      array-like shape (4,)
    weights: array-like shape (4,). Weights: Mp, g, sep, Fxuv
    """
    weights = np.array(weights)/np.linalg.norm(weights)

    points = np.asarray(coords, float)
    query = np.asarray(x, float)
    # if weights is not None:
    weights = np.asarray(weights, float)
    points = points * weights
    query = query * weights
    
    tree = cKDTree(points)
    distance, index = tree.query(query, k=1)
    return index, points[index], distance

def make_fname(x,complete=False):
    """
    x:          array-like shape (5,)
    complete:   Defaults to False. Adds '_complete' string to filename.
    """
    #savefile = f'H2O_{mp/const.Mearth:.3f}Me_{teq}K_{mstar/const.Msun:.3f}Ms_{Fxuv:.3f}Fxuv'
    
    str = f'saves/coarse_H2O_nobolo/{int(x[1])}/H2O_{x[0]:.3f}Me_{int(x[2])}K_{x[3].item():.3f}Ms_{x[4]:.3f}Fxuv'
    if complete:
        str = str+'_complete'
    
    return str+'.csv'

def initialize(init_params, input_grid):
    '''
    init_params: array-like shape (5,) of planet coordinates to initialize
    input_grid: list of mp, g, teq, mstar, Fxuv arrays ([],[],[],[],[])
    '''
    #load in inputs
    mp_ax =     input_grid[0]
    g_ax =      input_grid[1]
    teq_ax =    input_grid[2]
    ms_ax =     input_grid[3]
    Fxuv_ax =   input_grid[4]

    ms_ext = np.tile(ms_ax,len(teq_ax)) #extended ms
    teq_ext = np.repeat(teq_ax,len(ms_ax))
    mask_ms_low = ms_ext < 0.43
    lstar_ext = np.zeros_like(ms_ext)
    lstar_ext[mask_ms_low] = 0.23*(ms_ext[mask_ms_low])**2.3
    lstar_ext[~mask_ms_low] = ms_ext[~mask_ms_low]**4

    sep_ax = np.sqrt(lstar_ext*const.Lsun / (16*np.pi* sigma_sb * teq_ext**4)) #cm

    #check for duplicate sep in our parameterspace
    _, sep_cts = np.unique(sep_ax, return_counts=True)
    if np.any(sep_cts > 1):
        print("WARNING: duplicate semi-major axis values detected in the grid.\nModify your Teq or Mstar to avoid this.")
        exit()

    dim_mp =    len(mp_ax)
    dim_g =     len(g_ax)
    # dim_teq =   len(teq_ax)
    # dim_mstar = len(ms_ax)
    dim_Fxuv =  len(Fxuv_ax)
    dim_sep =   len(sep_ax)

    #Check for existing files
    exist_file = np.array([])
    all_index = np.array([])

    # for i,j,k,l,m in product(range(dim_mp),range(dim_g),range(dim_teq),range(dim_mstar),range(dim_Fxuv)):
    for i,j,k,l in product(range(dim_mp), range(dim_g), range(dim_sep), range(dim_Fxuv)):
        mp =    mp_ax[i]
        g =     g_ax[j]
        sep =   sep_ax[k]
        teq =   teq_ext[k]
        mstar = ms_ext[k]
        Fxuv =  Fxuv_ax[l] #m

        filename = make_fname([mp, g, teq, mstar, Fxuv], complete=True)
        # print(filename)

        coordinate = np.array([mp, g, sep, Fxuv])
        new_grid = [mp_ax,g_ax,sep_ax,Fxuv_ax]
        if os.path.isfile(filename):
            exist_file = np.append(exist_file, normalize(coordinate, new_grid),axis=0)
            all_index = np.append(all_index, coordinate, axis=0)

    #reshape arrays
    n_exists = int(len(exist_file)/4)
    exist_file = exist_file.reshape((n_exists,4))
    all_index = all_index.reshape((n_exists,4))

    init_idx,_,init_dist = nearest_in_cloud(
        exist_file,
        normalize(init_params, new_grid))
    closest = all_index[init_idx]
    sep_index = np.where(sep_ax == closest[2])
    closest_5d = [closest[0], closest[1], teq_ext[sep_index], ms_ext[sep_index], closest[3]]
    closest_file = make_fname(closest_5d,complete=True)

    # print("closest indexes :", init_idx)
    # print("closest coords :", closest)
    # print("closest coords norm :", init_pts)
    # print("nromalized distance :", init_dist)
    # print("filename :", make_fname(closest))

    return closest_file, init_dist

# Example weights: [1, 0.5, 1, 2, 1] (adjust as needed)
# closest_file, dist = initialize(init_params, input_grid, weights=[1, 0.5, 1, 2, 1])