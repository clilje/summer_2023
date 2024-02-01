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
work_dir = 'norm_catalog/'

#This is catalogue specific file naming convention
output_number = 17

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

df = df.nlargest(5, columns=['Mstar end SFR'])

#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    os.system('\n')
    
    #cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)
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

    # Normal dust
    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_elines'
    oname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_elines.sed'
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)

    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_everything'
    oname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_everything.sed'
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)
    

    # medium dust
    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_elines_moredust'
    oname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_elines_moredust.sed'
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)

    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_everything_moredust'
    oname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_everything_moredust.sed'
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)

    # Big dust
    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_elines_moredust_case2'
    oname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_elines_moredust_case2.sed'
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)

    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_everything_moredust_case2'
    oname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_everything_moredust_case2.sed'
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)



    # No dust
    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_elines_nodust'
    oname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_elines_nodust.sed'
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)

    cspname = 'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_everything_nodust'
    oname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_everything_nodust.sed'
    #Calls the pygalaxev function to create a file with all the input for galaxevpl
    pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
    os.system('$bc03/galaxevpl < %s'%tmpname)