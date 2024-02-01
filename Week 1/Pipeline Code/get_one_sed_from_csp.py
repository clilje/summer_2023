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

"""
This section covers all galaxies without taking their AGN into consideration.
"""
#Read in the galaxy specific data from the galaxy library that was previously created
if output_number == 17:
    df = pd.read_csv(work_dir+'library_galaxies.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])

else:
    df = pd.read_csv(work_dir+'library_galaxies'+str(output_number)+'.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])

#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    os.system('\n')
    
    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)
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

    oname = str(work_dir)+'/'+str(cspname)+'_age=%06.3f.sed'%age
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)

"""
This section covers only the galaxies with BH and considers their SEDs with the AGN added.
"""
#Read in the galaxy specific data from the galaxy library that was previously created
if output_number == 17:
    df = pd.read_csv(work_dir+'library_galaxies_with_BH.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])

else:
    df = pd.read_csv(work_dir+'library_galaxies_with_BH'+str(output_number)+'.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])

#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    
    #File naming conventions including _mod, to make sure the AGN contribution is considered
    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_mod'
    tmpname = work_dir+'/tmp.in'
    
    #age = 11. # time since beginning of star formation (in Gyr)
    #extracts galaxy age from library and then converts to Gyr
    age = df['age end SFR'][df['galaxy number'] == x]
    #print(age)
    age = age.to_numpy()[0]
    age = age*1e-9
    
    oname = str(work_dir)+'/'+str(cspname)+'_age=%06.3f_mod.sed'%age
    #print(cspname)
    
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)
