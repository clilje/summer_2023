"""
Author: Clara Lilje, adapted from pygalaxev code for personal use.


This file reads in a csp sed file and adds the appropriate AGN SED
to the galaxy SED using the method from Volonteri 2017.
And given a redshift calculates the magnitude of this 
sed in a series of filters.


"""

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
work_dir = 'TNG_50_library/'

#This code is specific to output 17
output_number = 17

#optional in case one wants specific galaxies
galaxy_number = [0,1,2,3]

#define list of filters as per BC03 filter dictionary
filter_list = [278,279,280,281,282,283,286,291]


#Dust parameters as a reminder
tau = 1.
tau_V = 0.15
mu = 0.4
epsilon = 0.

#cosmological parameters from Illustris TNG-50
h0 = 67.74
omega_m = 0.3089
omega_lambda = 0.6911


### Parameters for the AGN SED model by Volonteri 2017
nu_IR=0.01*3.2899*1e15 #0.01 Ryd to Hz
alpha_UV=-0.5
alpha_x=-1.0

# More parameters for the AGN SED model
# taken from Table1 of Thomas+16 paper
Ax2=-0.18
Ax1=0.0
Ax0=0.59
B=0.250
b3x1=0.034
b3x0=-0.019

b1x2=0.391
b1x1=-0.83
b1x0=0.61
ax2=-0.3
ax1=0.4
ax0=-0.82




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


def S(L_BH, M_BH):
  """
  Function to calculate the eddington ratio in MSun**(-1).
  Used in f1 and f2 to determine the Epeak as in the Thomas / Volonteri model.
  """
  # M_BH in Msun
  # L_BH/L_Edd the Eddington ratio 
  # Here our L is in real units. (erg/s)
  #L_Edd = M_BH*4.*np.pi*cte_G*cte_m_p*cte_c/cte_sigma_t
  L_Edd = 1.26*1e38 * M_BH # in erg/s
  S     = np.log10((L_BH/L_Edd)/M_BH)
  return S

def R(r_corona, M_BH):
  """
  returns log of r_corona in units of r_g
  used in f1 and f2 to determine Epeak as in Thomas / Volonteri model
  """
  # r_corona the coronal radius (inner edge of the visible portion of the disk). Should be in units of r_g
  # r_g the gravitational radius, r_g = cte_G*M_BH/cte_c**2
  #r_g   = cte_G*M_BH/cte_c**2
  #R     = np.log10(r_corona/r_g)
  R = np.log10(r_corona)
  return R

def f1(L_BH, M_BH, r_corona):
  """
  Used to determine Epeak (Thomas, 2016)
  """
  f1    = Ax2*R(r_corona, M_BH)*R(r_corona, M_BH) + Ax1*R(r_corona, M_BH) + Ax0 + B*S(L_BH,M_BH)
  return f1

def f2(L_BH, M_BH, r_corona):
  """
  Used to determine Epeak (Thomas, 2016)
  """
  f2    = (b3x1*R(r_corona, M_BH)+b3x0)*np.power(S(L_BH,M_BH)+6,3) + (b1x2*R(r_corona, M_BH)*R(r_corona, M_BH) + b1x1*R(r_corona, M_BH) + b1x0)*(S(L_BH, M_BH)+6) + (ax2*R(r_corona, M_BH)*R(r_corona, M_BH) + ax1*R(r_corona, M_BH) + ax0)
  return f2

def log10_E_peak(L_BH,M_BH,r_corona):
  """
  Function to determine log(Epeak) of the Big Blue bump as from Thomas 2016
  E_peak in keV
  L_BH in erg/s
  M_BH in Msun
  r_corona in r_g
  """
  # E_peak used to parametrize both the location in energy space and the shape of the Big Blue Bump of the disk emission.
  # E_peak in keV according to Thomas
  log10_E_peak = max(f1(L_BH, M_BH, r_corona),f2(L_BH, M_BH, r_corona))
  return log10_E_peak


def alpha_ox(L_BH,M_BH):
  """
  Slope of connection between two different power laws (?)
  as per thomas 2016
  """
  L_Edd     = 1.26*10**38*M_BH # in erg/s
  alpha_ox  =  -0.13*np.log(M_BH)+0.15*np.log10(L_BH/L_Edd)-0.45
  return alpha_ox

def F_nu(nu,L_BH,M_BH,r_corona):
  """
  Calculate flux of AGN as with the model of Volonteri 2016
  nu in Hz
  L_bH in erg/s
  M_bh in Msun
  r_corona in r_g
  """
  #This line is meant to convert from E_peak in kev to nu_peak in Hz
  numax_peak = 2.417e17*10**log10_E_peak(L_BH,M_BH,r_corona)    # in Hz 

  nu1=4.8359e17 # in Hz
  nu2=1.1992e15 # in Hz

  param = (nu1**(alpha_UV) * np.exp(-nu1/numax_peak) * np.exp(-nu_IR/nu1)-
           10**(alpha_ox(L_BH,M_BH)*2.6056) * nu2**alpha_UV * np.exp(-nu2/numax_peak)*
           exp(-nu_IR/nu2))/(10**(alpha_ox(L_BH,M_BH)*2.6056) * nu2**
                             alpha_x - nu1**alpha_x)
  if param < 0:
    param = (10**(alpha_ox(L_BH,M_BH)*2.6056) * nu2**alpha_UV * exp(-nu2/numax_peak)*
             exp(-nu_IR /nu2))/(nu1**alpha_x)


  #PARAM = 0. ## set to 0 for now.
  PARAM = param
  
  F_nu  = nu**alpha_UV * np.exp(-nu / numax_peak) * np.exp(-nu_IR / nu)  + PARAM*nu**alpha_x
  ##F_nu[np.where(nu<(0.1*3.2899e15))[0]] = nu[np.where(nu<(0.1*3.2899e15))[0]]**alpha_UV * np.exp(-nu[np.where(nu<(0.1*3.2899e15))[0]] / numax_peak) * np.exp(-nu_IR / nu[np.where(nu<(0.1*3.2899e15))[0]])
  for ii in range(len(nu)):
    if nu[ii]<3.2899e14: ### last term of Eq set to zero below 1.36 eV = 912 nm
      F_nu[ii]  = nu[ii]**alpha_UV * np.exp(-nu[ii] / numax_peak) * np.exp(-nu_IR / nu[ii])
  return F_nu

def F_lambda(F_nu,nu):
  # from F_nu to F_lambda
  #from (erg/s/Hz) to (erg/s/A)
  F_lambda = F_nu*(3e18)/((3e18/nu)**2)
  return F_lambda


def make_spectrum_AGN(wavelength,bh_mass,bh_lum):
    """
    This function creates the spectrum of the AGN using the Volonteri 2016 method
    The input is a given wavelength array in Armstrong
    BH_mass in Msun
    L_BH in log_10(BH_lum) in erg/s
    """
    L_BH  = 10**bh_lum # in erg/s
    M_BH  = bh_mass # in Msun
    L_Edd = 1.26*10**38*M_BH # in erg/s
    f_Edd_1 = np.log10(L_BH/L_Edd)
    
    nu = np.divide(3e18,wavelength) # from AA to Hz

    #Set radius for the moment
    r_corona = 10.

    #calculate flux values over given wavelength
    to_print_F_nu=F_nu(nu,L_BH,M_BH,r_corona) #in erg/s/Hz

    #calculate location of peak in Hz, from location of peak in kev
    nu_peak = 10**(log10_E_peak(L_BH,M_BH,r_corona))*2.42e17

    ## normalization using window of 6 orders of magnitude around the peak and the given luminosity
    norm_1=(L_BH)/scipy.integrate.quadrature(F_nu,nu_peak*1e-3,nu_peak*1e3,args=(L_BH,M_BH,r_corona))[0]

    #Return luminosity normalized flux in erg/s/A
    return(F_lambda(norm_1*to_print_F_nu, nu))




def add_spectrum(i, age_gal_i,bh_mass,bh_lum):
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
    ages = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in age_row)
    #select first value in row, this shows how many different ages are present in this file
    ages = ages[0].str.split(expand=True, n=1)
    #record number of ages in this specific file (changes between 238 and 263)
    num_ages = int(ages[0][0])

    #read in the ages row again, for a clean start
    ages = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in age_row)
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
    wavelength = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII', header=None, skiprows = lambda x: x not in specific_rows)
    #same procedure as above for age
    #split of first value, this is how many wavelength values are present in this file
    wavelength = wavelength[0].str.split(expand=True, n=1)
    #record the number of columns this should correspond to 
    num_col = int(wavelength[0][0])


    #re-read in the wavelength row
    wavelength = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII', header=None, skiprows = lambda x: x not in specific_rows)
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
    df = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows=7, skipfooter=12)
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
    In this section the AGN SED is created and added to every single SED in the file
    """
    #make AGN SED
    AGN = np.array(make_spectrum_AGN(xvals,bh_mass,bh_lum))
    # add the AGN spectra to the data frame as new column
    df['AGN'] = np.divide(AGN,3.826e+33)        #in units of L_sun /A

    #loop over all columns (this corresponds to number of ages recorded in the file)
    for y in range(int(len(ages))-1):
       
       df[y] = df.loc[:,[y,'AGN']].sum(axis=1)
    
    #drop AGN SED from Data frame again
    df = df.drop('AGN', axis=1)

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
    front = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in  [0,1,2,3,4,5,6])
    
    ###this section gets the "footer" section of the file in its undisturbed form, to join together
    back = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows = int(len(ages))+7)
    
    #trial
    #back = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in [281,280,279,278,277,276,275,274,273,272,271,270])
    #back = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows=rows_age[-1]+1)


    ### now we join together header, main body and footer to one data frame and write to modified ascii file
    frames = [front,middle,back]
    result = pd.concat(frames)
    result.to_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'_mod.ised_ASCII',header=None,index=None)



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
    os.system('$bc03/ascii_ised TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'.ised')




def make_ised(output_number, galaxy_number):
    """
    This function calls the bc03 bin_ised function to transform an .ised_ASCII file for a 
    given galaxy id and output number into an .ised file.
    """    
    #cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_mod.ised_ASCII'
    os.system('$bc03/bin_ised TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_mod.ised_ASCII')


"""
Main section of the code calling the functions used above
"""

#specific to output number 17, should be extended if extended use is required
if output_number == 17:
    redshift = 5

#read in library file which includes all galaxies with a BH
df = pd.read_csv(work_dir+'library_galaxies_with_BH.txt',header=None,
                    names=['galaxy number','age end SFR','Mstar end SFR','BH mass','BH Lum'])


galaxy_number = df['galaxy number'].to_numpy()
#galaxy_number = [2]

#loop over all galaxies from the library by ID
for x in galaxy_number:
    
    #get key galaxy values from the library file
    age = df['age end SFR'][df['galaxy number'] == x].values[0]
    mstar = df['Mstar end SFR'][df['galaxy number'] == x].values[0]
    bh_mass = df['BH mass'][df['galaxy number'] == x].values[0]
    bh_lum = df['BH Lum'][df['galaxy number'] == x].values[0]

    #transform .ised file from this galaxy to ascii
    make_ascii(output_number, x)

    #create modified .ascii file which adds spectrum of AGN
    add_spectrum(x,age,bh_mass,bh_lum)

    #convert modified .ascii file to .ised for further use
    make_ised(output_number, x)

    #name of modified file
    cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_mod'

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

    #convert age to Gyr as requested in input
    age_today = age_today*1e-09

    #define output of modified magnitude file
    outname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_mod'
    get_mags(age_today, filter_list, isedname,outname)
