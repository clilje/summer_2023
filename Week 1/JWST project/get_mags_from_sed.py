"""
Author: Clara Lilje, adapted from pygalaxev code for personal use.


This file reads in a csp sed file and given a redshift calculates the magnitude of this 
sed in a series of filters.

"""

import numpy as np
import h5py
import pygalaxev_cosmology
from pygalaxev_cosmology import L_Sun, Mpc, c as csol
from scipy.interpolate import splrep, splev, splint
import os
import pandas as pd



#Parameters for file naming conventions
work_dir = 'TNG_50_library/'
output_number = 17
#galaxy_number = [0,1,2,3,4,5]

#Dust absorption parameters as a reminder
tau = 1.
tau_V = 0.15
mu = 0.4
epsilon = 0.

#cosmological parameters from TNG 50
h0 = 67.74
omega_m = 0.3089
omega_lambda = 0.6911

#Numbers of JWST filters as taken from BC03 overview file
filter_list = [278,279,280,281,282,283,286,291]

def get_mags(age_today, filters, filename,outname, input_tmpname='tmp.in', output_tmpname='mySSP'):
    """
    This function uses the mm_evolution function from the BC03 code to produce magnitudes in given filters.
    """

    #create temporary file which contains the input for the function
    tmpname = work_dir+'/%s'%input_tmpname

    inputlines = []

    # input cosmology
    inputlines.append(str(h0)+','+str(omega_m)+','+str(omega_lambda)+'\n')
    inputlines.append(str(age_today)+'\n')
    inputlines.append('%s\n'%filename)

    #input chosen filters
    filterstring =','.join([str(x) for x in filters])
    inputlines.append(filterstring+'\n')

    f = open(tmpname, 'w')
    f.writelines(inputlines)
    f.close()

    # Run bc03
    os.system('$bc03/mm_evolution < %s'%tmpname)
    #os.system('csp < %s'%tmpname)
    #os.system('/src/csp_galaxev < %s'%tmpname)


    # Keep the output *.multi_mag_AB file as a *.magnitudes file
    os.system('scp %s.multi_mag_AB %s.magnitudes'%(outname, outname))
    os.system('rm -f %s.multi_mag_vega'%( outname))
    os.system('rm -f %s.multi_mag_AB'%( outname))
    # Clean up
    os.system('rm -f %s/%s*'%(work_dir, output_tmpname))


#Specific to dictionary, should be extended if other output numbers get used
if output_number == 17:
    redshift = 5

#Read library file to get info on all galaxies (all considered without a BH)
df = pd.read_csv(work_dir+'library_galaxies.txt',header=None,
                    names=['galaxy number','age end SFR','Mstar end SFR','BH mass','BH Lum'])


#Loop through all galaxies by index
for x in df['galaxy number'].to_numpy():

    #Name of SED file to be considered
    cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(x)
    
    #get age from library
    age = df['age end SFR'][df['galaxy number'] == x].values[0]
    #mstar = df['Mstar end SFR'][df['galaxy number'] == x].values[0]
    isedname = str(cspname)+'.ised'
    
    """
    In this section we calculate that the galaxy has exactly the age it has when it stops star forming at redshift 5.
    This is SPECIFIC TO OUTPUT NUMBER 17. 
    NEEDS TO BE ADJUSTED FOR OTHER REDSHIFTS/OUTPUT NUMBERS.
    """
    #Only works for redshift 5 for the moment
    if redshift == 5:
        #cosmological calculation
        #age_today = (13.798e+09-1.173e+09) + age

        #lookback time to redshift 5 as taken from BC03 with above cosmological parameters
        age_today = 12.636e+09 + age
    
    #convert age to Gyr
    age_today = age_today*1e-09

    #define output name for magnitude file    
    outname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(x)

    #call function to apply mm_evolution to the given .ised file with given filters at given age
    get_mags(age_today, filter_list, isedname,outname)
