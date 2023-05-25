import pygalaxev
import numpy as np
from scipy.interpolate import splrep, splev
import os
import h5py
import pandas as pd


work_dir = 'TNG_50_library/'
output_number = 17
galaxy_number = 2

tau = 1.
tau_V = 0.15
mu = 0.4
epsilon = 0.

cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)

df = pd.read_csv(work_dir+'library_galaxies.txt',header=None,
                       names=['galaxy number','age end SFR','Mstar end SFR'])
print(df)

#age = 11. # time since beginning of star formation (in Gyr)
age = df['age end SFR'][df['galaxy number'] == galaxy_number].values[0]

# Create the mass normalization models
massname = work_dir+'/'+cspname+'.mass'
d = np.loadtxt(massname)
mass_spline = splrep(d[:, 0], d[:, 10], k=3, s=0) #using the sum of M*_liv+M_rem to renormalize the mass

# extracts SED corresponding to the given age

tmpname = work_dir+'/tmp.in'

oname = work_dir+'/'+cspname+'_age=%06.3f.sed'%age

pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
os.system('$bc03/galaxevpl < %s'%tmpname)

f = open(oname, 'r')
wsed = np.loadtxt(f)
f.close()

wave = wsed[:, 0]
flux = wsed[:, 1]
           
# Renormalize the mass!
logAge = np.log10(age)+9.
mass = splev(logAge, mass_spline)
sed = flux/mass

# store SED to an .hdf5 file

#output_file = h5py.File(oname.replace('.sed', '.hdf5'), 'w')
#wave_dset = output_file.create_dataset('wave', data=wave)
#wave_dset.attrs['units'] = 'Angstrom'

#sed_dset = output_file.create_dataset('Llambda', data=sed)
#sed_dset.attrs['description'] = 'Luminosity density (dL/dlambda) for a 1 Solar Mass (living + remnants) stellar population'
#sed_dset.attrs['units'] = 'L_Sun/Angstrom'

