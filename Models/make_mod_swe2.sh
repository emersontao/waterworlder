# This is a shell script to create the iop interpolator, written in fortran, wrapped for python

# Compile the "mod_fint.f90" fortran module (general "fortran interpolator")
# Will create "mod_fint.o"
gfortran -c -fPIC mod_fint.f90

# Compile the "mod_iop.f95" module (the actual iop interpolator)
# Will create "mod_fint.mod" and "mod_iop.o"
gfortran -c -fPIC mod_swe2.f95

# Make the wrapper using "fortran 2 python"
# Will create "iop_fint.cpython-312-darwin.so" ("312" is the python version)
f2py -c -m swe2025_fint mod_fint.f90 mod_swe2.f95

# ----------------------------------
# Guide to use it in python scripts:
# ----------------------------------
# Add the path of the folder containing "iop_fint.cpython-312-darwin.so" 
# as a location to check when importing modules
# >>> import sys
# >>> sys.path.insert(0, '/Users/aaguichine/Models')

# Import iop_fint
# >>> import iop_fint as iop

# At the beginning of your code, load the data in the RAM
# >>> iop.iop.load_data("/Users/aaguichine/Models/Aguichine2021_mr_all_2024.dat")

# Use it like before. I don't think this one can be parallelized, 
# >>> radius = iop.iop.interpol_iop(cmf,wmf,teq,mass)
