# KDE Plot for gpevo
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import seaborn as sns
# sns.set_theme(style="darkgrid")

# Number of samples to generate for each dataset
n_planets = 5000

M_e = 5.97237e24 #Earth mass (kg)
R_e = 6.371e6 #Earth radius (m)

path_models=  './Models/'

# des_file = "Models/desert_boundaries.txt"
# cont_data = np.loadtxt(des_file, skiprows=1)
# px_cont = cont_data[:, 0]
# ry_cont = cont_data[:, 1]
# lvl_cont = cont_data[:, 2]
# seg_cont = cont_data[:, 3]

# mask = (
#     (px_cont >= -0.30) & (px_cont <= 0.60)) #&
#     # (ry_cont >= -0.05) & (ry_cont <= min(RY_MAX, np.log10(R_DESERT_MAX)))
# # )

# px_cont = px_cont[mask]
# ry_cont = ry_cont[mask]
# lvl_cont = lvl_cont[mask]
# seg_cont = seg_cont[mask]

# for lev in np.unique(lvl_cont):
#     mask_lev = lvl_cont == lev
#     seg_ids = np.unique(seg_cont[mask_lev])
#     for sid in seg_ids:
#         m = mask_lev & (seg_cont == sid)
#         px_seg = px_cont[m]
#         ry_seg = ry_cont[m]
#         if len(px_seg) < 2:
#             continue
#         plt.plot(px_seg, ry_seg, color="k", linewidth=4.2, alpha=0.95, zorder=4)
# plt.xticks(
#         [-0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4],
#         [
#             round(10 ** -0.4, 1),
#             "1.0",
#             round(10 ** 0.4, 1),
#             round(10 ** 0.8, 1),
#             round(10 ** 1.2, 1),
#             round(10 ** 1.6, 1),
#             round(10 ** 2.0, 1),
#             round(10 ** 2.4, 1),
#         ],
#         fontsize=13.5,
#     )

# plt.yticks(
#     [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
#     [
#         "1.0",
#         round(10 ** 0.2, 1),
#         round(10 ** 0.4, 1),
#         round(10 ** 0.6, 1),
#         round(10 ** 0.8, 1),
#         round(10 ** 1.0, 1),
#         round(10 ** 1.2, 1),
#         round(10 ** 1.4, 1),
#     ],
#     fontsize=13.5,
# )

# plt.xlim(-0.5, 2.2)
# plt.ylim(0.07, 1.36)
# plt.xlabel("Orbital period (days)", fontsize=16)
# plt.ylabel(r"Planet radius ($\rm R_{\oplus}$)", fontsize=16)

# plt.tick_params(
#     which="both",
#     direction="in",
#     length=6,
#     width=1.3,
#     color="k",
#     top=True,
#     right=True,
#     bottom=True,
#     left=True,
# )


# # des_x, des_y, = np.loadtxt(des_file,usecols=(0,1), skiprows=1, unpack=True)
# plt.plot(10**des_x, 10**des_y)
# plt.xscale("log")
plt.xlim(None,30)
plt.ylim(1,4)
# loop import planets starting and ending 
mp_start = np.array([])
rp_start = np.array([])
mp_end = np.array([])
rp_end = np.array([])

datafolder = './paper/5000_xp/H2O'
period_gp = np.loadtxt(datafolder+'/5000_xp_gp_data.csv',skiprows=(1),usecols=(11),unpack=True)

for i in range(n_planets):
    planetnum = i+1
    # import a planet's evolutionary track (just mp_gp and rp_gp)
    filename = f'{datafolder}/planet_{planetnum}-{n_planets}_evolution_data.csv'
    mp_gp,rp_gp = np.loadtxt(filename,usecols=(1,2),unpack=True)
    
    # mp_s = mp_gp[0]
    # mp_e = mp_gp[-1]
    # rp_s = rp_gp[0]
    # rp_e = rp_gp[-1]

    mp_start = np.append(mp_start, mp_gp[0])
    rp_start = np.append(rp_start, rp_gp[0])
    mp_end = np.append(mp_end, mp_gp[-2])
    rp_end = np.append(rp_end, rp_gp[-2])



##############################################
#
#   LOAD DATA
#
##############################################


import requests
import os
from datetime import datetime
import argparse

# Define command-line arguments
parser = argparse.ArgumentParser(description="Run the mardigras tool with optional features.")

# Flag to update the NEA catalog
parser.add_argument(
    "--update-nea-catalog",
    action="store_true",
    help="Update the NASA Exoplanet Archive catalog before starting the tool.",
)

# Flag to update the PlanetS catalog
parser.add_argument(
    "--update-planets-catalog",
    action="store_true",
    help="Update the PlanetS catalog before starting the tool.",
)

# Flag to choose another type of catalog
parser.add_argument(
    "--catalog",
    choices=["NEA", "PlanetS"],
    default="PlanetS",
    help="Choose the exoplanet catalog to use. Default is PlanetS."
)

args = parser.parse_args()



# NASA Exoplanet catalog

def check_internet_connection():
    """Check if there is internet access by pinging a known URL."""
    try:
        requests.get("https://www.google.com", timeout=5)
        return True
    except requests.ConnectionError:
        return False

def update_nea_exoplanet_catalog(catalog_url, output_file):
    """
    Updates the NEA exoplanet catalog from the NASA Exoplanet Archive TAP.
    Parameters:
        catalog_url (str): The TAP URL with the SQL query for the catalog.
        output_file (str): Path to save the downloaded catalog.
    """
    if not check_internet_connection():
        print("No internet connection. Unable to update the exoplanet catalog.")
        return

    try:
        response = requests.get(catalog_url, timeout=10)
        response.raise_for_status()  # Raise an error for HTTP issues
        
        # Prepare the header
        header = [
            "# NASA Exoplanet Catalog",
            f"# Source: {catalog_url}",
            f"# Catalog last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        ]
        
        # Write the header and data to the file
        with open(output_file, "w") as f:
            f.write("\n".join(header) + "\n")
            # Add '#' to the parameter names
            data_lines = response.text.splitlines()
            f.write("# " + data_lines[0] + "\n")  # Add # to parameter names
            f.write("\n".join(data_lines[1:]) + "\n")
        
        print(f"Catalog updated successfully and saved to {output_file}")
    except requests.RequestException as e:
        print(f"Error fetching the catalog: {e}")

def read_nea_last_update(output_file):
    """
    Reads the date of the last update from the catalog file and prints it.
    Parameters:
        output_file (str): Path to the catalog file.
    """
    if not os.path.exists(output_file):
        print("Catalog file does not exist.")
        return

    try:
        with open(output_file, "r") as f:
            for line in f:
                text_lookup = "# Catalog last updated:"
                if line.startswith(text_lookup):
                    print("# NEA catalog last updated: ",line.replace(text_lookup,'').replace('\n',''))
                    break
    except Exception as e:
        print(f"Error reading the catalog: {e}")

# File paths and URL of the NEA Catalog
nea_catalog_url = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
               "query=select+pl_name,pl_rade,pl_radeerr1,pl_radeerr2,"
               "pl_masse,pl_masseerr1,pl_masseerr2,pl_eqt+from+ps+where+"
               "default_flag=1+and+pl_controv_flag=0+and+pl_rade+is+not+null+"
               "and+pl_masse+is+not+null+and+pl_bmassprov='Mass'&format=tsv")
nea_output_file = "./resources/MARDIGRAS/data/catalog_exoplanets.dat"

# Update NEA catalog
if args.update_nea_catalog:
    # Update the catalog
    update_nea_exoplanet_catalog(nea_catalog_url, nea_output_file)
else:
    # Check for existing catalog and print the last update date
    read_nea_last_update(nea_output_file)


# PlanetS catalog

def update_planets_exoplanet_catalog(catalog_url, output_file):
    """
    Updates the PlanetS exoplanet catalog from the DACE website.
    Warning: Work in progress, currently not working
    """
    if True:
        print("PlanetS auto update currently not implemented, skipping updating attempt.")
        return

    if not check_internet_connection():
        print("No internet connection. Unable to update the exoplanet catalog.")
        return

    try:
        response = requests.get(catalog_url, timeout=10)
        response.raise_for_status()  # Raise an error for HTTP issues
        
        # Prepare the header
        header = [
            "# NASA Exoplanet Catalog",
            f"# Source: {catalog_url}",
            f"# Catalog last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        ]
        
        # Write the header and data to the file
        with open(output_file, "w") as f:
            f.write("\n".join(header) + "\n")
            # Add '#' to the parameter names
            data_lines = response.text.splitlines()
            f.write("# " + data_lines[0] + "\n")  # Add # to parameter names
            f.write("\n".join(data_lines[1:]) + "\n")
        
        print(f"Catalog updated successfully and saved to {output_file}")
    except requests.RequestException as e:
        print(f"Error fetching the catalog: {e}")

def read_planets_last_update(output_file):
    """
    Reads the date of the last update from the catalog file and prints it.
    Parameters:
        output_file (str): Path to the catalog file.
    """
    if not os.path.exists(output_file):
        print("Catalog file does not exist.")
        return

    try:
        with open(output_file, "r") as f:
            for line in f:
                text_lookup = "# Catalog last updated:"
                if line.startswith(text_lookup):
                    print("# PlanetS catalog last updated: ",line.replace(text_lookup,'').replace('\n',''))
                    break
    except Exception as e:
        print(f"Error reading the catalog: {e}")

# File paths and URL of the PlanetS Catalog
planets_catalog_url = ("None")
planets_output_file = "./Models/catalog_exoplanets_planets.dat"

# Update PlanetS catalog
if args.update_planets_catalog:
    # Update the catalog
    update_planets_exoplanet_catalog(planets_catalog_url, planets_output_file)
else:
    # Check for existing catalog and print the last update date
    read_planets_last_update(planets_output_file)

# Load exoplanet catalog
if args.catalog == "NEA":
    # Load NEA data
    data_nea = np.genfromtxt(nea_output_file,delimiter="\t",unpack=True,usecols=(1,2,3,4,5,6),filling_values=0.0)

    list_catalog_rp = data_nea[0]
    list_catalog_rpe1 = data_nea[1]
    list_catalog_rpe2 = data_nea[2]
    list_catalog_mp = data_nea[3]
    list_catalog_mpe1 = data_nea[4]
    list_catalog_mpe2 = data_nea[5]

    # Load the exoplanet names
    list_catalog_names = np.genfromtxt(
        nea_output_file, delimiter="\t", dtype=str, usecols=0
    )

    # Remove lines with empty Mp, Rp, or mass precision less than 50%
    exo_mass_prec = 0.5
    filter_nea = (list_catalog_rp!=0) & (list_catalog_mp!=0.0) \
    & ( abs(list_catalog_mpe1)/list_catalog_mp < exo_mass_prec) \
    & ( abs(list_catalog_mpe2)/list_catalog_mp < exo_mass_prec)

    list_exo_rp = list_catalog_rp[filter_nea]
    list_exo_rpe1 = list_catalog_rpe1[filter_nea]
    list_exo_rpe2 = list_catalog_rpe2[filter_nea]
    list_exo_mp = list_catalog_mp[filter_nea]
    list_exo_mpe1 = list_catalog_mpe1[filter_nea]
    list_exo_mpe2 = list_catalog_mpe2[filter_nea]
    list_exo_names = list_catalog_names[filter_nea]
elif args.catalog == "PlanetS":
    # Load PlanetS data
    data_planets = np.genfromtxt(planets_output_file,delimiter="\t",unpack=True,usecols=(1,2,3,4,5,6),filling_values=0.0)

    list_catalog_rp = data_planets[0]
    list_catalog_rpe1 = data_planets[1]
    list_catalog_rpe2 = data_planets[2]
    list_catalog_mp = data_planets[3]
    list_catalog_mpe1 = data_planets[4]
    list_catalog_mpe2 = data_planets[5]

    # Load the exoplanet names
    list_catalog_names = np.genfromtxt(
        planets_output_file, delimiter="\t", dtype=str, usecols=0
    )

    # Remove lines with empty Mp and Rp
    exo_mass_prec = 0.5
    filter_planets = (list_catalog_rp!=0) & (list_catalog_mp!=0.0)

    list_exo_rp = list_catalog_rp[filter_planets]
    list_exo_rpe1 = list_catalog_rpe1[filter_planets]
    list_exo_rpe2 = list_catalog_rpe2[filter_planets]
    list_exo_mp = list_catalog_mp[filter_planets]
    list_exo_mpe1 = list_catalog_mpe1[filter_planets]
    list_exo_mpe2 = list_catalog_mpe2[filter_planets]
    list_exo_names = list_catalog_names[filter_planets]



## Plot parameters
xmin, xmax, nx = 0.5, 19.5, 50 # 0.5, 20
ymin, ymax  = 0.75, 4.25
# xmin, xmax     = 2, 13
# ymin, ymax     = 1, 3.25

x = np.logspace(np.log10(xmin), np.log10(xmax), nx)


#MARDIGRAS stuff

xmin, xmax, nx = 0.1, 19.99999, 50 # 0.5, 20
ymin, ymax  = .75, 4.25

##############################################
#
#   MAKE SWEET INTERPOLATOR (Aguichine+2025)
#
##############################################

path_sweet = "Models/Aguichine2025_SWEET_all.dat"

# top of the atmosphere pressure (index)
sweet_top = np.array([0,1])
sweet_labels_top = ["20 mbar", "1 µbar"]
# host star type
sweet_labels_host = ["Type M", "Type G"]
# equilibrium temperature at zero albedo (K)
sweet_teqs = np.array([400, 500, 700, 900, 1100, 1300, 1500])
# bulk water mass fraction (%)
sweet_wmfs = np.array([0.1,1,10,20,30,40,50,60,70,80,90,100])
# total planet mass (Me)
sweet_masses = np.array([0.2       ,  0.254855  ,  0.32475535,  0.41382762,  0.52733018,
        0.67196366,  0.85626648,  1.09111896,  1.39038559,  1.77173358,
        2.25767578,  2.87689978,  3.66596142,  4.67144294,  5.95270288,
        7.58538038,  9.66586048, 12.31696422, 15.69519941, 20.        ])
# planet age (Gyr)
sweet_ages = np.array([0.001,0.0015,0.002,0.003,0.005,0.01,
                        0.02,0.03,0.05,
                        0.1,0.2,0.5,
                        1.0,2.0,5.0,
                        10,20])

# make dimension axes for interpolator
sweet_dim_wmf = sweet_wmfs.copy()/100
sweet_dim_teq = sweet_teqs.copy()
sweet_dim_mass = sweet_masses.copy()
sweet_dim_age = sweet_ages.copy()
sweet_dim_star = np.array([0,1])
sweet_dim_top = np.array([0,1])

# read radius tracks data and reshape
listrpfull1 = np.loadtxt(path_sweet,skiprows=36,unpack=True,usecols=(6))
listrpfull2 = np.loadtxt(path_sweet,skiprows=36,unpack=True,usecols=(7))
listrpfull = np.array((listrpfull1,listrpfull2))
listrpfull = np.delete(listrpfull, np.arange(17, listrpfull.size, 18))

sweet_data_radius = np.reshape(listrpfull,(2,2,12,7,20,17))
mask = np.isnan(sweet_data_radius)
sweet_data_radius[mask] = -1

# fill_value: value to return when extrapolating beyond grid validity range
fill_value = None
interp_sweet = RegularGridInterpolator((sweet_dim_top,
                                      sweet_dim_star,
                                      sweet_dim_wmf,
                                      sweet_dim_teq,
                                      sweet_dim_mass,
                                      sweet_dim_age), 
                                     sweet_data_radius, method='slinear', bounds_error=False, fill_value=fill_value)


# The parametrized function to be plotted
def radius_m(x,wmf_swe,teq_swe,age_swe):
    input0 = np.stack([np.full(len(x),0), np.full(len(x),0), np.full(len(x),wmf_swe), np.full(len(x),teq_swe), x, np.full(len(x),age_swe)],axis=-1)
    # input0 = np.stack(0, 0, (np.full(len(x),wmf_swe),np.full(len(x),teq_swe),x,np.full(len(x),age_swe)), axis=-1)
    return interp_sweet(input0)

def radius_g(x,wmf_swe,teq_swe,age_swe):
    input0 = np.stack([np.full(len(x),0), np.full(len(x),1), np.full(len(x),wmf_swe), np.full(len(x),teq_swe), x, np.full(len(x),age_swe)],axis=-1)
    # input0 = np.stack(0, 1, (np.full(len(x),wmf_swe),np.full(len(x),teq_swe),x,np.full(len(x),age_swe)), axis=-1)
    return interp_sweet(input0)

# Define initial parameters
init_host_swe = 0
init_wmf_swe = 0.5
init_teq_swe = 600.0
init_age_swe = 1.0

init_wmf_swe100 = 1.0

# line_swe_m, = plt.plot(
#                        x, rad_swe_m(x,init_wmf_swe,init_teq_swe,init_age_swe),color='blue',zorder=19)
# line_swe_g, = plt.plot(
#     x, radius_g(x,init_wmf_swe,init_teq_swe,init_age_swe),lw=0.75,linestyle=':',color='teal',zorder=19)

# line_swe_m100, = plt.plot(
#                           x, rad_swe_m(x,init_wmf_swe100,init_teq_swe,init_age_swe),color='blue',zorder=19)
# line_swe_g100, = plt.plot(
    # x, radius_m(x,init_wmf_swe100,init_teq_swe,init_age_swe),lw=0.75,color='teal',zorder=19)

plt.plot(x, radius_m(x,init_wmf_swe,init_teq_swe,init_age_swe),lw=0.75,color='teal',linestyle='-.',zorder=19,label='50% Steam')
plt.plot(x, radius_g(x,init_wmf_swe100,init_teq_swe,init_age_swe),lw=0.75,color='teal',zorder=19,label='100% Steam')

#swe stuff
# swe_func = [rad_swe_m, rad_swe_g]
swe_labels = ["Type M", "Type G"]  # Optional labels for each 


##############################################
# Load Zeng data
path_zeng = path_models + "./Zeng2016.dat"

# Load curves from Zeng et al. 2016
zeng_mass,zeng_purefe,zeng_rock,zeng_50wat,zeng_100wat,zeng_earth = np.loadtxt(path_zeng,delimiter="\t",skiprows=1,unpack=True)

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
line_zeng, = plt.plot(x, rad_zeng(x,init_cmf_zeng),lw=2,color='brown',alpha=0.5,zorder=-10)

# Zeng fixed
zeng_width  = 1.0
# plt.plot(list_zeng_masses,list_zeng_fe_radii,linewidth=zeng_width,color='black',zorder=5)
plt.plot(list_zeng_masses,list_zeng_ea_radii,linewidth=zeng_width,color='brown',alpha=0.75,zorder=6)
plt.plot(list_zeng_masses,list_zeng_mg_radii,linewidth=zeng_width,color='grey',alpha=0.75,zorder=7)
# plt.plot(zeng_mass,zeng_50wat,linewidth=zeng_width,color='darkturquoise',alpha=0.75,zorder=-8)
# plt.plot(zeng_mass,zeng_100wat,linewidth=zeng_width,color='darkturquoise',alpha=0.75,zorder=-9)

# Zeng fixed labels
plt.text(0.55, 0.60, '100% Core',fontsize=5,
        color='black',rotation=10)
        
plt.text(2.1, 1.26, 'Earth-like',fontsize=5,
        color='brown',rotation=26)

plt.text(2.1, 1.36, '100% Mantle',fontsize=5,
        color='dimgrey',rotation=26)

plt.text(2.1, 1.6, '50% Liquid H2O',fontsize=5,
        color='teal',rotation=28)

plt.text(2.1, 1.77, '100% Liquid H2O',fontsize=5,
        color='teal',rotation=32)


# ------------------------------------------------ #
# ADD SOLAR SYSTEM

# Data for masses and radii of planets (in Earth masses and Earth radii)
ss_planet_data = {
    "Mercury": (0.055, 0.383),
    "Venus": (0.815, 0.949),
    "Earth": (1.0, 1.0),
    "Mars": (0.107, 0.532),
    "Jupiter": (317.8, 11.21),
    "Saturn": (95.2, 9.45),
    "Uranus": (14.6, 4.01),
    "Neptune": (17.2, 3.88),
}

# Alchemy symbols for planets
ss_alchemy_symbols = {
    "Mercury": "☿",
    "Venus": "♀",
    "Earth": "⊕",
    "Mars": "♂",
    "Jupiter": "♃",
    "Saturn": "♄",
    "Uranus": "♅",
    "Neptune": "♆",
}

# Extracting data
ss_planets = list(ss_planet_data.keys())
ss_masses = [ss_planet_data[planet][0] for planet in ss_planets]
ss_radii = [ss_planet_data[planet][1] for planet in ss_planets]
ss_symbols = [ss_alchemy_symbols[planet] for planet in ss_planets]

for i in range(len(ss_planets)):
    text = plt.text(ss_masses[i], ss_radii[i], ss_alchemy_symbols[ss_planets[i]], \
            fontsize=15, ha='center', va='center', color='black')
    #text.set_path_effects([patheffects.withStroke(linewidth=1, foreground='black')])


# ------------------------------------------------ #
# KDE PLOT

i_color = sns.cm.crest
f_color = sns.cm.flare


# unicorn_killer = period_gp > 1

# # Plot KDE for Dataset 1
sns.kdeplot(x=mp_start/M_e, y=rp_start/R_e, cmap=i_color, log_scale=True, \
            label="Initial", levels=5)
sns.kdeplot(x=mp_start/M_e, y=rp_start/R_e, log_scale=True, fill=True, \
            cmap=i_color, levels=5, alpha=0.25, legend=True)

# Plot KDE for Dataset 2 #log_scale=True, 
sns.kdeplot(x=mp_end/M_e, y=rp_end/R_e, cmap=f_color, log_scale=True, \
            label="Final", levels=5)
sns.kdeplot(x=mp_end/M_e, y=rp_end/R_e, log_scale=True, fill=True, \
            cmap=f_color, levels=5, alpha=0.25,legend=True)


# Planets
list_exo_mpe = [abs(list_exo_mpe2), list_exo_mpe1]
list_exo_rpe = [abs(list_exo_rpe2), list_exo_rpe1]

# Exoplanets from the catalog file
catalog_points = plt.errorbar(list_exo_mp,list_exo_rp,
            yerr=list_exo_rpe,
            xerr=list_exo_mpe,
            fmt='o',zorder=-10,elinewidth=1,
            c="slategray",label='Exoplanet',alpha=0.15)  #slategrey
# plt.show()
# plt.axhline(1.8, color='black', linestyle='--', linewidth=.75) #1.8 R_e radius valley



# Add labels and title
plt.xlabel(r'Mass [$M_\oplus$]')
plt.ylabel(r'Radius [$R_\oplus$]')
plt.xscale('log')
plt.yscale('linear')
plt.xticks(ticks=[1,3,5,10,20],labels=['1','3','5','10','20'])
plt.yticks(ticks=[1,2,3,4],labels=['1','2','3','4'])
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.legend(fontsize=7)

# import matplotlib.patches as mpatches
# handles = [mpatches.Patch(facecolor=sns.cm.i_color(100), label="Initial"),
#            mpatches.Patch(facecolor=sns.cm.f_color(100), label="Final")]  #sns.cm.crest(100)   sns.cm.flare(100)

# plt.legend(handles=handles, fontsize=7)

filename = f'{n_planets}_kde'
plt.savefig(f'./paper/plots/{filename}.png',dpi=600)
print(f'  File saved: {os.getcwd()}/{filename}.png')

# Show the plot
plt.show()
