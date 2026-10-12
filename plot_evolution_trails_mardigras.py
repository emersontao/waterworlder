# Scatter & Evolution trails
import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, './Models')
from scipy.interpolate import RegularGridInterpolator
# import seaborn as sns
import swe, xpml # type: ignore
# sns.set_theme(style="darkgrid")

# Earth params
M_e = 5.97237e24 #Earth mass (kg)
R_e = 6.371e6 #Earth radius (m)

# Number of samples to generate for each dataset
n_planets = 148


## Plot parameters
xmin, xmax, nx = 0.05, 19.5, 50 # 0.5, 20
ymin, ymax  = 0.75, 4.25
# xmin, xmax     = 2, 13
# ymin, ymax     = 1, 3.25

x = np.logspace(np.log10(xmin), np.log10(xmax), nx)

##############################################
#
#   MAKE SWE INTERPOLATOR (Aguichine+2024)
#
##############################################

# Load the data
path_models = "./Models/"

path_swe = path_models+"A24_SWE_all_v3_fixed.dat"

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

listrpfull = np.loadtxt(path_swe,skiprows=36,unpack=True,usecols=(5))
listrpfull = np.delete(listrpfull, np.arange(17, listrpfull.size, 18))

mask = listrpfull == -1.0
listrpfull[mask] == np.inf

listrpfull_m = listrpfull[0:int(len(listrpfull)/2)]
listrpfull_g = listrpfull[int(len(listrpfull)/2):]

# Make SWE interpolator

swe_dim_wmf = np.array(swe_wmfs)/100
swe_dim_teq = np.array(swe_teqs)
swe_dim_mass = np.array(swe_masses)
swe_dim_age = swe_ages

swe_data_radius_m = np.reshape(listrpfull_m,(12,7,20,17))
swe_data_radius_g = np.reshape(listrpfull_g,(12,7,20,17))

fill_value = None
interp_swe_m = RegularGridInterpolator((swe_dim_wmf, swe_dim_teq, swe_dim_mass, swe_dim_age), swe_data_radius_m,
                                       method='slinear', bounds_error=False, fill_value=fill_value)
interp_swe_g = RegularGridInterpolator((swe_dim_wmf, swe_dim_teq, swe_dim_mass, swe_dim_age), swe_data_radius_g,
                                       method='slinear', bounds_error=False, fill_value=fill_value)



# The parametrized function to be plotted
def rad_swe_m(x,wmf_swe,teq_swe,age_swe):
    input0 = np.stack((np.full(len(x),wmf_swe),np.full(len(x),teq_swe),x,np.full(len(x),age_swe)), axis=-1)
    return interp_swe_m(input0)

def rad_swe_g(x,wmf_swe,teq_swe,age_swe):
    input0 = np.stack((np.full(len(x),wmf_swe),np.full(len(x),teq_swe),x,np.full(len(x),age_swe)), axis=-1)
    return interp_swe_g(input0)


# Define initial parameters
init_host_swe = 0
init_wmf_swe = 0.5
init_teq_swe = 600.0
init_age_swe = 1.0

init_wmf_swe100 = 1.0

# line_swe_m, = plt.plot(
#                        x, rad_swe_m(x,init_wmf_swe,init_teq_swe,init_age_swe),color='blue',zorder=19)
line_swe_g, = plt.plot(
    x, rad_swe_g(x,init_wmf_swe,init_teq_swe,init_age_swe),lw=0.75,linestyle=':',color='teal',zorder=19)

# line_swe_m100, = plt.plot(
#                           x, rad_swe_m(x,init_wmf_swe100,init_teq_swe,init_age_swe),color='blue',zorder=19)
line_swe_g100, = plt.plot(
    x, rad_swe_g(x,init_wmf_swe100,init_teq_swe,init_age_swe),lw=0.75,color='teal',zorder=19)

#swe stuff
swe_func = [rad_swe_m, rad_swe_g]
swe_labels = ["Type M", "Type G"]  # Optional labels for each 


##############################################
# Load Zeng data
# path_zeng = path_models + "./Zeng2016.dat"

# Load curves from Zeng et al. 2016
# zeng_mass,zeng_purefe,zeng_rock,zeng_50wat,zeng_100wat,zeng_earth = np.loadtxt(path_zeng,delimiter="\t",skiprows=1,unpack=True)

# open files
list_zeng_fe_m,list_zeng_fe_r = np.loadtxt(path_models+"zeng2016-iron.dat",unpack=True,usecols=(0,1))
list_zeng_ea_m,list_zeng_ea_r = np.loadtxt(path_models+"zeng2016-earth.dat",unpack=True,usecols=(0,1))
list_zeng_mg_m,list_zeng_mg_r = np.loadtxt(path_models+"zeng2016-rock.dat",unpack=True,usecols=(0,1))

list_zeng_masses = np.logspace(np.log10(0.01), np.log10(100.0), 30)

# reshape radius data to the same vector
list_zeng_fe_radii = np.interp(list_zeng_masses, list_zeng_fe_m, list_zeng_fe_r)
list_zeng_ea_radii = np.interp(list_zeng_masses, list_zeng_ea_m, list_zeng_ea_r)
list_zeng_mg_radii = np.interp(list_zeng_masses, list_zeng_mg_m, list_zeng_mg_r)

# create interpolator
dimcmf_zeng = np.array([0.0,0.325,1.0])
data_zeng = np.vstack((list_zeng_mg_radii, list_zeng_ea_radii,list_zeng_fe_radii)).T
interp_zeng = RegularGridInterpolator((list_zeng_masses,dimcmf_zeng), data_zeng, method='linear', bounds_error=False, fill_value=None)

# parameterized func to plot
def rad_zeng(x,cmf):
    input0 = np.stack((x,np.full(len(x),cmf)), axis=-1)
    rp = interp_zeng(input0)
    return rp

init_cmf_zeng = 0.325
# line_zeng, = plt.plot(x, rad_zeng(x,init_cmf_zeng),lw=2,color='brown',alpha=0.5,zorder=-10)

# Zeng fixed
zeng_width = 0.75
# plt.plot(list_zeng_masses,list_zeng_fe_radii,linewidth=zeng_width,color='black',zorder=5)
plt.plot(list_zeng_masses,list_zeng_ea_radii,linewidth=zeng_width,color='brown',alpha=0.75,zorder=6)
plt.plot(list_zeng_masses,list_zeng_mg_radii,linewidth=zeng_width,color='grey',alpha=0.75,zorder=7)
# plt.plot(zeng_mass,zeng_50wat,linewidth=zeng_width,color='darkturquoise',alpha=0.75,zorder=-8)
# plt.plot(zeng_mass,zeng_100wat,linewidth=zeng_width,color='darkturquoise',alpha=0.75,zorder=-9)

# Zeng fixed labels
# plt.text(0.55, 0.60, '100% Core',fontsize=5,
#         color='black',rotation=10)
        
# plt.text(2.1, 1.26, 'Earth-like',fontsize=5,
#         color='brown',rotation=26)

# plt.text(2.1, 1.36, '100% Mantle',fontsize=5,
#         color='dimgrey',rotation=26)

# plt.text(2.1, 1.6, '50% Liquid H2O',fontsize=5,
#         color='teal',rotation=28)

# plt.text(2.1, 1.77, '100% Liquid H2O',fontsize=5,
#         color='teal',rotation=32)




plt.grid(visible=True, color='slategrey', alpha=0.2, which='major', axis='both')

# //////////////////////////////////////// #
# Evolution Track Scatterplot
# Creates a scatterplot of exoplanet starting params and evolutionary trails.

#import starting data
# gp_file = f'./evolution-data/gp-evo/100-5000_xp/{n_planets}_xp_gp_data.csv'
# mp_init, rp_init, wmf_init = np.loadtxt(gp_file,skiprows=1,usecols=(1,2,3),unpack=True)

#-- NEA data
gp_file = './evolution-data/real-evo/nea/PS_2025.07.16_21.12.59.csv'
period_gp,sep_gp,rp_gp,mp_gp,teq_gp,mstar_gp,lstar_gp = np.loadtxt(gp_file,dtype='float',delimiter=',',skiprows=(41),usecols=(4,5,6,7,10,15,16),unpack=True)


wmf_gp = np.array([])
rsurf_gp = np.array([])
age_init = 1e-2 * xpml.Gyr
#loop to load and plot data
for i in range(145):
    planetnum = i + 1
    ''''''
    #because we're using NEA
    if mstar_gp[i]/xpml.M_sun > 0.43:
        # Host star is M-TYPE
        wmf = swe.solve_wmf_g(rp_gp[i],teq_gp[i],mp_gp[i],age_init)
    else:
        # Host star is G-TYPE
        wmf = swe.solve_wmf_m(rp_gp[i],teq_gp[i],mp_gp[i],age_init)
    wmf_gp = np.append(wmf_gp, wmf)
    rsurf = swe.zeng_radius(mp_gp/xpml.M_e,0.325) * xpml.R_e
    rsurf_gp = np.append(rsurf_gp, rsurf)
    ''''''

    #Load planet data
    filename = f'./evolution-data/real-evo/nea/planet_{planetnum}-{n_planets}_evolution_data.csv' #load nea
    # filename = f'./evolution-data/gp-evo/100-5000_xp/planet_{planetnum}-{n_planets}_evolution_data.csv' #load gp
    time_arr,mp_arr,rp_arr,wmf_arr = np.loadtxt(filename,usecols=(0,1,2,3),unpack=True)
    
    zorder = 50+planetnum
    plt.plot(mp_arr[0]/M_e, rp_arr[0]/R_e, marker='o', c='#F38457', markersize=4, markeredgecolor='none', alpha=0.5, zorder=zorder) #, alpha=0.75
    # plt.annotate(planetnum, (mp_gp[0], rp_gp[0]), textcoords="offset points", xytext=(5,-10),
    #             ha='left',bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white', alpha=0.7),zorder=50+i)
    plt.plot(mp_arr[:-2]/M_e, rp_arr[:-2]/R_e, c='orange', alpha=0.25, zorder=zorder)
    plt.plot(mp_arr[-2]/M_e, rp_arr[-2]/R_e, marker='o', c='#ADBD53', markersize=3.5, markeredgecolor='none', alpha=0.5, zorder=zorder) #, alpha=0.75
# plt.legend()
# # alternatively
# sns.scatterplot(data=df,x='Mp', y='Rp', size=2) #, hue='WMF', palette="crest")

#load final plot
time_arr,mp_arr,rp_arr,wmf_arr = np.loadtxt(f'./evolution-data/real-evo/nea/planet_146-{n_planets}_evolution_data.csv',usecols=(0,1,2,3),unpack=True)
# time_gp,mp_gp,rp_gp,wmf_gp = np.loadtxt(f'./evolution-data/gp-evo/100-5000_xp/planet_100-{n_planets}_evolution_data.csv',skiprows=1,usecols=(0,1,2,3),unpack=True)
plt.plot(mp_arr[0]/M_e, rp_arr[0]/R_e, marker='o', c='#F38457', markersize=4, markeredgecolor='none', label='Initial State', alpha=0.5, zorder=50) #, alpha=0.75
plt.plot(mp_arr[:-2]/M_e, rp_arr[:-2]/R_e, c='orange', label='Evolution track', alpha=0.25, zorder=50)
plt.plot(mp_arr[-2]/M_e, rp_arr[-2]/R_e, marker='o', c='#ADBD53', markersize=3.5, markeredgecolor='none', label='Final State', alpha=0.5, zorder=50) #, alpha=0.75


#plot details
plt.title(f'146 NEA Steam Worlds Mass-Radius Evolution')
#plt.set_figheight(4)
#plt.set_figwidth(6)

plt.xlabel('Mass (M_e)')
plt.ylabel('Radius (R_e)')
plt.xlim(0.25,12.25)
plt.ylim(0.75,4.25)
plt.legend(loc='upper right', fontsize=7)

# //////////////////////////////////////// #
# SAVE ALL GENERATED DATA & CORRESPONDING PLOT
# identifier= '100-5000xp_1-12Me_500-1200K'
# savefile = f'{n_planets}-{identifier}_evolution_trails_lines'
savefile = f'{n_planets}_evolution_trails_lines'
plt.savefig(f'./{savefile}.png', dpi=600)
# print(f'File saved: {os.getcwd()}\{filename}.png \n')

plt.show()


plt.figure()
plt.scatter(mp_gp[:len(wmf_gp)],rp_gp[:len(wmf_gp)],c=wmf_gp, alpha=0.25, cmap='viridis')
plt.xlim(0.25,12.25)
plt.ylim(0.75,4.25)
plt.colorbar(label='WMF')
plt.show()