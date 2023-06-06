import pygalaxev
import numpy as np
from scipy.interpolate import splrep, splev
import os
import h5py
import pandas as pd


work_dir = 'TNG_50_library/'
output_number = 17
galaxy_number = [0,1,2,3,4,5]

tau = 1.
tau_V = 0.15
mu = 0.4
epsilon = 0.

df = pd.read_csv(work_dir+'library_galaxies.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])
for x in df['galaxy number'].to_numpy():
    os.system('\n')
    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)
    print(df)

    #age = 11. # time since beginning of star formation (in Gyr)
    age = df['age end SFR'][df['galaxy number'] == x]
    print(age)
    # extracts SED corresponding to the given age

    tmpname = work_dir+'/tmp.in'

    oname = str(work_dir)+'/'+str(cspname)+'_age=%06.3f.sed'%age

    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)

df = pd.read_csv(work_dir+'library_galaxies_with_BH.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])
for x in df['galaxy number'].to_numpy():
    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_mod'
    tmpname = work_dir+'/tmp.in'
    
    #age = 11. # time since beginning of star formation (in Gyr)
    age = df['age end SFR'][df['galaxy number'] == x]
    print(age)

    oname = str(work_dir)+'/'+str(cspname)+'_age=%06.3f_mod.sed'%age
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)
