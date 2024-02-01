"""
Author: Clara Lilje, adapted from pygalaxev code for personal use.

This file will use the galaxevpl function of the BC03 code to produce the 
SED for a galaxy at a specific age.

"""

import pygalaxev
import numpy as np
from scipy.interpolate import splrep, splev
import os
import h5py
import pandas as pd

#Define parameters of workspace
work_dir = 'TNG_50_library/'

#This is catalogue specific file naming convention
output_number = 11

#Not actually used in this code, but can be used in case one wants to run the code for specific galaxies
#galaxy_number = [0,1,2,3,4,5]

#These are the parameters for the dust absorption model in BC03.
tau = 1.
tau_V = 0.15
mu = 0.4
epsilon = 0.

tau_V = [0.1,0.5,1,2,3,4,5,6]
mu = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 ]
"""
This section covers all galaxies without taking their AGN into consideration.
"""
#Read in the galaxy specific data from the galaxy library that was previously created
df = pd.read_csv(work_dir+'library_galaxies'+str(output_number)+'.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])

#loop through different galaxies, indexing files by galaxy number
for tt in tau_V:
  for mm in mu:
    for x in [384]:
        os.system('\n')
        
        cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'tau='+str(tt)+'mu='+str(mm)
        #print(df)

        #age = 11. # time since beginning of star formation (in Gyr)

        #extracts galaxy age from library and then converts to Gyr
        age = df['age end SFR'][df['galaxy number'] == x]
        #print(age)
        age = age.to_numpy()[0]
        age = age*1e-9


        # extracts SED corresponding to the given age
        #pygalaxev file naming conventions
        tmpname = work_dir+'/tmp.in'

        oname = str(work_dir)+'/'+str(cspname)+'tau='+str(tt)+'mu='+str(mm)+'_age=%06.3f.sed'%age
        #Calls the pygalaxev function to create a file with all the input for galaxevpl
        pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
        os.system('$bc03/galaxevpl < %s'%tmpname)
