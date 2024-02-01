"""
Author: Clara Lilje, adapted from pygalaxev code for personal use.


This file reads in a csp sed file and given a redshift calculates the magnitude of this 
sed in a series of filters.

"""

from doctest import OutputChecker
import numpy as np
import h5py
import pygalaxev
import pygalaxev_cosmology
from pygalaxev_cosmology import L_Sun, Mpc, c as csol
from scipy.interpolate import splrep, splev, splint
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mtp
import math
import time
import numpy as np
from numpy import *
from collections import Counter
from scipy import interpolate
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
import pylab
from matplotlib import rc
from matplotlib.colors import LogNorm
import scipy
from scipy.optimize import curve_fit

#Parameters defining file structure
work_dir = 'norm_catalog/'

#This code is specific to output 17
output_number = 17

#optional in case one wants specific galaxies
galaxy_number = [2]

#define list of filters as per BC03 filter dictionary
filter_list = [276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292]

#Dust parameters as a reminder
tau = 1.
tau_V = 3
mu = 0.4
epsilon = 0.

#cosmological parameters from Illustris TNG-50
h0 = 67.74
omega_m = 0.3089
omega_lambda = 0.6911


### Parameters for the AGN S
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



def normalize_mag(i, age_gal_i,mstar_TNG):
    """
    This function adds an appropriate spectrum to the galaxy seds in the .ised_ASCII file.
    Therefore the galaxy sed file needs to be converted to an ascii file first. (you can use the make_ascii 
    function before)
    i = galaxy ID
    age_gal_i = age at which galaxy is done with SFRH
    bh_mass = mass of bh in galaxy (needs to be non-zero) in Msun
    bh_lum = Luminosity of BH in galaxy (needs to be non-zero) in log_10(Luminosity) in erg/s
    """

    """
    This section finds the ages at which the sed is recorded in the ascii file.
    """
    #select only row 1, in this row the ages are recorded
    age_row = [0]
    #read in row 1
    ages = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in age_row)
    #select first value in row, this shows how many different ages are present in this file
    ages = ages[0].str.split(expand=True, n=1)
    #record number of ages in this specific file (changes between 238 and 263)
    num_ages = int(ages[0][0])

    #read in the ages row again, for a clean start
    ages = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in age_row)
    #seperate into expected number of different ages
    ages = ages[0].str.split(expand=True, n=num_ages)
    #drop first value, this is just the number of values in this row
    ages = ages.drop(0,axis=1)
    ages = np.float64(ages.to_numpy()[0])
    

    """
    In this section we read in how many different wavelengths are covered in this file
    """
    #Wavelengths are recorded in row 7
    specific_rows = [6]
    #read in row 7
    wavelength = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust.ised_ASCII', header=None, skiprows = lambda x: x not in specific_rows)
    #same procedure as above for age
    #split of first value, this is how many wavelength values are present in this file
    wavelength = wavelength[0].str.split(expand=True, n=1)
    #record the number of columns this should correspond to 
    num_col = int(wavelength[0][0])


    #re-read in the wavelength row
    wavelength = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust.ised_ASCII', header=None, skiprows = lambda x: x not in specific_rows)
    #split into expected number of different values
    wavelength = wavelength[0].str.split(expand=True, n=num_col)
    #drop first column, this only contains the number of different wavelength values
    wavelength = wavelength.drop(0, axis=1)
    #transpose dataframe to convert to numpy
    wavelength = wavelength.T
    #get the wavelengths as x-array, this is to then give to the AGN function
    xvals = wavelength[0].to_numpy()
    xvals = np.float64(xvals)


    """
    In this section we seperate out all the rows which contain SEDs from the .ASCII file
    """
    #in case one wants to select specifically the row with the age
    #index_row = bisection(ages, age_gal_i)
    #rows_age = [index_row+5, index_row+6, index_row+7]

    #These are all the rows which do NOT contain SEDs, for reference, for a normal 282 line file
    headerplustail_rows = [0,1,2,3,4,5,6,281,280,279,278,277,276,275,274,273,272,271,270]
    
    #trials
    #df = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows=headerplustail_rows)
    #df = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII', header=None, skiprows=lambda x: x not in rows_age)
    
    #skip bottom 12 and top 7 rows to only select the SEDs, same in every file, regardless how many timesteps
    df = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust.ised_ASCII',delimiter='/s', header=None, skiprows=7, skipfooter=12)
    #split every row into the expected number of columns, so we can seperate out only the values that are luminosities
    df = df[0].str.split(expand=True, n=num_col+1)
    
    #Drop first and last column, these do not contain luminosity values
    #But keep them stored to add them back on to re-construct the file
    first_col = df[0]
    last_col = df[num_col+1]

    df = df.drop(0, axis=1)
    df = df.drop((num_col+1), axis=1)

    #transform each value to float and transpose for better editing (pandas likes columns more than rows)
    df = df.astype(float)
    df = df.T


    """
    In this section the BC03 mass and the tng mass are compared.
    """
    df_mass = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(int(i))+'_lowmass_evenmoredust.mass',delimiter='/s', header=None, skiprows=30)
    df_mass = df_mass[0].str.split(expand=True)
    age = df_mass[0].astype('float64').to_numpy()
    mstar_live = df_mass[5].astype('float64').to_numpy()    #live
    sfr_live = df_mass[9].astype('float64').to_numpy()
    index = np.where(sfr_live==0)[0][0]
    #make normalization constant
    mstar_BC03 = mstar_live[index-1]
    print(mstar_BC03)
    print(mstar_TNG)

    # add the AGN spectra to the data frame as new column
    #df['normalisation'] = np.ones(len(df[0].to_numpy()))*(mstar_TNG/mstar_BC03)        #in units of L_sun /A
    norm = (mstar_TNG/mstar_BC03)
    print(norm)
    df = df.astype(float)*norm
    #loop over all columns (this corresponds to number of ages recorded in the file)
    #for y in range(int(len(ages))-1):
    #   
    #   df[y] = df.loc[:,[y,'normalisation']].multiply(axis=1)
    
    #drop AGN SED from Data frame again
    #df = df.drop('normalisation', axis=1)

    """
    In this section we reconstruct the data frame so that we can reconvert the ASCII file in the same format
    """
    #re-transpose so we have original columns again
    df = df.T

    #re-insert the two columns we previously removed
    df.insert(0, 0, first_col, True)
    df.insert(num_col+1, num_col+1, last_col, True)

    #join all columns to be one column of string per row, with appropriate delimiter to match ASCII file
    middle = df.astype(str)
    middle = middle.apply(lambda x: '       '.join(x.dropna()), axis=1)

    ###this section gets the "header" section of the file in its undisturbed form, to join together
    #trial
    #skipfooter = 282-rows_age[0] -(262-int(first_col[0]))
    #front = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skipfooter= 282-rows_age[0])
    #get first 7 lines
    front = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in  [0,1,2,3,4,5,6])
    
    ###this section gets the "footer" section of the file in its undisturbed form, to join together
    back = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust.ised_ASCII',delimiter='/s', header=None, skiprows = int(len(ages))+7)
    
    #trial
    #back = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in [281,280,279,278,277,276,275,274,273,272,271,270])
    #back = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows=rows_age[-1]+1)


    ### now we join together header, main body and footer to one data frame and write to modified ascii file
    frames = [front,middle,back]
    result = pd.concat(frames)
    result.to_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust.ised_ASCII',header=None,index=None)



def bisection(array,value):
    """
    This is a bisection function taken from stackexchange
    """
    '''Given an ``array`` , and given a ``value`` , returns an index j such that ``value`` is between array[j]
    and array[j+1]. ``array`` must be monotonic increasing. j=-1 or j=len(array) is returned
    to indicate that ``value`` is out of range below and above respectively.'''
    n = len(array)
    if (value < array[0]):
        return -1
    elif (value > array[n-1]):
        return n
    jl = 0# Initialize lower
    ju = n-1# and upper limits.
    while (ju-jl > 1):# If we are not yet done,
        jm=(ju+jl) >> 1# compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm# and replace either the lower limit
        else:
            ju=jm# or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
    if (value == array[0]):# edge cases at bottom
        return 0
    elif (value == array[n-1]):# and top
        return n-1
    else:
        return jl



def make_ascii(output_number, galaxy_number):
    """
    This function calls the bc03 ascii_ised function to transform an .ised file for a 
    given galaxy id and output number into an ASCII file.
    """
    #cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'.ised'
    os.system('$bc03/ascii_ised norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_lowmass_evenmoredust.ised')




def make_ised(output_number, galaxy_number):
    """
    This function calls the bc03 bin_ised function to transform an .ised_ASCII file for a 
    given galaxy id and output number into an .ised file.
    """    
    #cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_mod_nodust.ised_ASCII'
    os.system('$bc03/bin_ised norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_lowmass_evenmoredust.ised_ASCII')
    os.system('rm -f norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_lowmass_evenmoredust.ised_ASCII')
    os.system('rm -f norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_lowmass_evenmoredust.ised_ASCII')


#Specific to dictionary, should be extended if other output numbers get used
if output_number == 17:
    redshift = 5
    #Read library file to get info on all galaxies (all considered without a BH)
    df = pd.read_csv(work_dir+'library_galaxies'+str(output_number)+'_lowmass.txt',header=None,
                        names=['galaxy number','age end SFR','Mstar end SFR','BH mass','BH Lum'])
if output_number == 11:
    redshift = 7
    #Read library file to get info on all galaxies (all considered without a BH)
df = pd.read_csv(work_dir+'library_galaxies'+str(output_number)+'_lowmass.txt',header=None,
                 names=['galaxy number','age end SFR','Mstar end SFR','BH mass','BH Lum'])


#Loop through all galaxies by index
for x in df['galaxy number'].to_numpy():
    #for x in [39361]:
    #Name of SED file to be considered
    cspname = 'norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_lowmass_evenmoredust'
    age = df['age end SFR'].astype(int)[df['galaxy number'] == x].values[0]
    print(age)
    mstar = df['Mstar end SFR'][df['galaxy number'] == x].values[0]
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
    if redshift == 7:
        age_today = 13.047e+09 + age

    #convert age to Gyr
    age_today = age_today*1e-09
    
    
    outname = 'norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_lowmass_evenmoredust'

    #transform .ised file from this galaxy to ascii
    make_ascii(output_number, x)

    #create modified .ascii file which adds spectrum of AGN
    normalize_mag(x,age,mstar)

    #convert modified .ascii file to .ised for further use
    make_ised(output_number, x)

    #call function to apply mm_evolution to the given .ised file with given filters at given age
    get_mags(age_today, filter_list, isedname,outname)
