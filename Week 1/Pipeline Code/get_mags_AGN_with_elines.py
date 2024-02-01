"""
Author: Clara Lilje, adapted from pygalaxev code for personal use.


This file reads in a csp sed file and adds the appropriate AGN SED
to the galaxy SED using the method from Volonteri+17.
And given a redshift calculates the magnitude of this 
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
import bisect

#Parameters defining file structure
work_dir = 'norm_catalog/'

#This code is specific to output 17
output_number = 17

#optional in case one wants specific galaxies
galaxy_number = [39361]

#define list of filters as per BC03 filter dictionary
filter_list = [276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292]

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


#convert line width
cte_c_kms=2.99792458e5 ## km s-1





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


### Useful functions
def sigma_to_fwhm(sigma):
    """
    Convert from Gaussian sigma to FWHM
    """
    # can be simplified with FWHM = sigma * 2.355
    # sigma is the standard deviation.
    return sigma * 2.0 * np.sqrt(2.0 * np.log(2.0))


def fwhm_to_sigma(fwhm):
    """
    Convert FWHM of 1D Gaussian to sigma
    """
    return fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def height_to_amplitude(height, sigma):
    """
    Convert height of a 1D Gaussian to the amplitude
    """
    return height * sigma * np.sqrt(2 * np.pi)


def amplitude_to_height(amplitude, sigma):
    """
    Convert amplitude of a 1D Gaussian to the height
    """
    return amplitude / (sigma * np.sqrt(2 * np.pi))

def restframe_wl(x, z):
    """
    Transform a given spectrum x into the restframe, given the redshift
    Input:
        x: observed wavelengths of a spectrum, in Angstrom or nm for example
        z: redshift
    Return:
        restframe spectrum
    """
    return x / (1.0 + z)


def restframe_flux(x, z):
    """
    Transform a given spectrum flux density x into the restframe, given the
    redshift
    Input:
        x: observed flux density or standard deviation of a spectrum, in
        erg/s/cm^2/A for example
        z: redshift
    Return:
        restframe spectrum flux density
    """
    return x * (1.0 + z)

def gaussian(x,peak, sigma, height_flux):
    return height_flux*np.exp(-4*np.log(2)*np.power(x - peak, 2.) / (np.power(sigma, 2.))/(sigma*np.sqrt(np.pi/np.log(2))/2.))

def get_elines(wavelength,Z,bh_lum,x_id,n_h,log_Us,a):
  #print(Z)
    
  Zcode_dic = {0.0001:'Z0001',0.001:'Z001', 0.0002:'Z0002', 0.002:'Z002', 0.004:'Z004', 0.0005:'Z0005', 0.006:'Z006', 0.008:'Z008', 0.014:'Z014', 0.017:'Z017', 0.02:'Z020', 0.03:'Z030', 0.04:'Z040', 0.05:'Z050', 0.06:'Z060', 0.07:'Z070'}

  #choose the name of the metallicity that most closely corresponds
  res = bisect.bisect_left(list(Zcode_dic.keys()), Z)
  if res == 7:
      Zcode  = list(Zcode_dic.values())[res-1]
  elif res == 1:
      Zcode = list(Zcode_dic.values())[0]
  else:
      if abs((list(Zcode_dic.keys())[res]-Z)) < abs(list(Zcode_dic.keys())[res-1]-Z):
          Zcode = list(Zcode_dic.values())[res]
      else:
          Zcode = list(Zcode_dic.values())[res-1]


  # -----------Initial parameters
  # Z=0.008 (metallicity of the gas taken half solar)
  logUs, xid, nh, alpha,  AGN_line_OII_3727, AGN_line_Hbeta, AGN_line_OIII_4959, AGN_line_OIII_5007, AGN_line_OI_6300, AGN_line_NII_6548, AGN_line_Halpha,  AGN_line_NII_6584, AGN_line_SII_6717, AGN_line_SII_6731, AGN_line_NV_1240, AGN_line_CIV_1548, AGN_line_CIV_1551, AGN_line_HeII_1640, AGN_line_OIII_1661, AGN_line_OIII_1666, AGN_line_SiIII_1883, AGN_line_SiIII_1888, AGN_line_CIII_1907, AGN_line_CIII_1910 =np.loadtxt('Emission_lines/nlr_nebular_'+str(Zcode)+'.txt',skiprows=1,unpack=True)
  line_name=[AGN_line_OII_3727, AGN_line_Hbeta, AGN_line_OIII_4959, AGN_line_OIII_5007, AGN_line_OI_6300, AGN_line_NII_6548, AGN_line_Halpha, AGN_line_NII_6584, AGN_line_SII_6717, AGN_line_SII_6731, AGN_line_NV_1240, AGN_line_CIV_1548, AGN_line_CIV_1551, AGN_line_HeII_1640, AGN_line_OIII_1661, AGN_line_OIII_1666, AGN_line_SiIII_1883, AGN_line_SiIII_1888, AGN_line_CIII_1907, AGN_line_CIII_1910]
  line_wavelength=[3727,4861,4959,5007,6300,6548,6563, 6584,6717,6731,1240,1548,1551,1640,1661,1666,1883,1888,1907,1910]
  modifiedarray_flux=[0.]*len(wavelength)
  ### FIX WITH PROPER FLUX AND SELECT LINE!!!!
  #line_flux=[1e44]*len(line_wavelength) #[line_name[pp][0]*1e47/1e45 for pp in range(len(line_wavelength))] # the amplitude of the line, to be taken from Gutkin+16
  #for xx in range(len(logUs)):
    #if logUs[xx] == round(log_Us*2)/2 and alpha[xx]== a and xid[xx] == x_id and nh[xx] == n_h:
        #line_flux=[line_name[pp][xx]/1e45*(10**bh_lum)/line_wavelength[pp] for pp in range(len(line_wavelength))] 
        #print(line_flux)
        #line_flux=[line_name[pp][xx]*3.826e+33*sfr/line_wavelength[pp] for pp in range(len(line_wavelength))] 
        #print(line_flux)
  line_sigma=300 #in km/s ## corresponds to sigma (standard deviation of the gaussian)
  # Full Width at Half Maximum (FWHM) is the width measured at half level between the continuum and the peak of the line.



  # -----------Initial arrays from BC03
  array_wavelength=wavelength ## should remain untouched
  #array_flux=flux            ## should remain untouched
  #print('ini',len(array_flux))
  #print(array_wavelength)
  # Modified arrays accounting for the emission lines
  modifiedarray_flux=[0.]*len(wavelength) #np.copy(array_flux)

  # -----------Loop over the different emission lines
  #print(line_name)
  #print(line_name[:,[0]])
  #print(len(line_name[:,[0]]))
  for ppp in np.arange(len(line_wavelength)):
    #print(ppp)
    # Convert line width from km/s to AA
    #line_width = line_wavelength[ppp] * line_sigma * 1e3 / cte_c_kms
    
    for xx in range(len(logUs)):
        if logUs[xx] == round(log_Us*2)/2 and alpha[xx]== a and xid[xx] == x_id and nh[xx] == n_h:
            line_width = line_wavelength[ppp] * line_sigma / (2.355*cte_c_kms)
            line_height=line_name[ppp][xx]/1e45*(10**bh_lum)
                    #print(line_flux)
            #print('line_width in AA',line_width)
            #print(line_name[ppp])
            # Identify the closest wavelength (in array_wavelength) of the emission line.
            #difference_array = np.absolute(array_wavelength*0.0001*(1+z)-line_wavelength[ppp]*0.0001*(1+z))
            #index_wavelength=array_wavelength[difference_array.argmin()]
            #index_flux=flux[difference_array.argmin()]

            # Select a small wavelength range (in array_wavelength) around the emission line.
            filter = (array_wavelength > line_wavelength[ppp]-4*line_width) & (array_wavelength < line_wavelength[ppp]+4*line_width)
            #array_inter_wavelength=array_wavelength[filter]
            #array_inter_flux=array_flux[filter]

            #print('line = ',line_wavelength[ppp])
            #print('line_width=',line_width)
            #print('end',len(array_inter_flux))
            #print(filter)
            print(line_width,line_height)
            for kk in range(len(filter)):
                if filter[kk]==True:
                    #print('hi')
                    modifiedarray_flux[kk] += gaussian(array_wavelength[kk],line_wavelength[ppp], line_width, line_height)
          #if modifiedarray_flux[kk]!=array_flux[kk]:
          #  print(modifiedarray_flux[kk],array_flux[kk])



    ###CHECK IF YOU NEED TO CONVERT THIS TO RESTFRAME
  return(modifiedarray_flux)


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




def add_spectrum(i, age_gal_i,bh_mass,bh_lum,met_gal,av_sfr,x_id,n_h,log_Us,a):
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
    ages = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust_BH_cont.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in age_row)
    #select first value in row, this shows how many different ages are present in this file
    ages = ages[0].str.split(expand=True, n=1)
    #record number of ages in this specific file (changes between 238 and 263)
    num_ages = int(ages[0][0])

    #read in the ages row again, for a clean start
    ages = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust_BH_cont.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in age_row)
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
    wavelength = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust_BH_cont.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in specific_rows)
    #same procedure as above for age
    #split of first value, this is how many wavelength values are present in this file
    wavelength = wavelength[0].str.split(expand=True, n=1)
    #record the number of columns this should correspond to 
    num_col = int(wavelength[0][0])


    #re-read in the wavelength row
    wavelength = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust_BH_cont.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in specific_rows)
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
    df = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust_BH_cont.ised_ASCII',delimiter='/s', header=None, skiprows=7, skipfooter=12)
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
    #AGN = np.array(make_spectrum_AGN(xvals,bh_mass,bh_lum))
    # add the AGN spectra to the data frame as new column
    #df['AGN'] = np.divide(AGN,3.826e+33)        #in units of L_sun /A

    elines = get_elines(xvals,met_gal,bh_lum,x_id,n_h,log_Us,a)
    df['emission lines'] = np.divide(elines,3.826e+33)        #in units of L_sun /A

    #loop over all columns (this corresponds to number of ages recorded in the file)
    for y in range(int(len(ages))-1):
       
       #df[y] = df.loc[:,[y,'AGN']].sum(axis=1)
       df[y] = df.loc[:,[y,'emission lines']].sum(axis=1)
    
    #drop emission lines SED from Data frame again
    df = df.drop('emission lines', axis=1)
    
    #drop AGN SED from Data frame again
    #df = df.drop('AGN', axis=1)

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
    front = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust_BH_cont.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in  [0,1,2,3,4,5,6])
    
    ###this section gets the "footer" section of the file in its undisturbed form, to join together
    back = pd.read_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust_BH_cont.ised_ASCII',delimiter='/s', header=None, skiprows = int(len(ages))+7)
    
    #trial
    #back = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in [281,280,279,278,277,276,275,274,273,272,271,270])
    #back = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows=rows_age[-1]+1)


    ### now we join together header, main body and footer to one data frame and write to modified ascii file
    frames = [front,middle,back]
    result = pd.concat(frames)
    result.to_csv('norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'_lowmass_evenmoredust_BH_elines.ised_ASCII',header=None,index=None)



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
    os.system('$bc03/ascii_ised norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_lowmass_evenmoredust_BH_cont.ised')




def make_ised(output_number, galaxy_number):
    """
    This function calls the bc03 bin_ised function to transform an .ised_ASCII file for a 
    given galaxy id and output number into an .ised file.
    """    
    #cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_mod.ised_ASCII'
    os.system('$bc03/bin_ised norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_lowmass_evenmoredust_BH_elines.ised_ASCII')
    os.system('rm -f norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_lowmass_evenmoredust_BH_cont.ised_ASCII')
    os.system('rm -f norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_lowmass_evenmoredust_BH_elines.ised_ASCII')


"""
Main section of the code calling the functions used above
"""

if output_number == 17:
    redshift = 5
    #read in library file which includes all galaxies with a BH
if output_number == 11:
   redshift = 7
   #read in library file which includes all galaxies with a BH
df = pd.read_csv('norm_catalog/library_galaxies_with_BH'+str(output_number)+'_lowmass.txt',header=None,
                 names=['galaxy number','age end SFR','Mstar end SFR','BH mass','BH Lum'])

df_mod = pd.read_csv('TNG_50_library/RescalingBHs_case2.txt',header=None,
                    names=['galaxy number','BH mass','BH Lum'],delimiter='  ')
#read in library file which includes all galaxies with a BH
#df = pd.read_csv('TNG_50_library/library_galaxies_with_BH'+str(output_number)+'.txt',header=None,
#                    names=['galaxy number','age end SFR','Mstar end SFR','BH mass','BH Lum'])


galaxy_number = df['galaxy number'].to_numpy()
#galaxy_number = [39361]
xid = 0.1
nh = 100
alpha = -1.7

#loop over all galaxies from the library by ID
for x in galaxy_number:
    
    filename_gal='Files_SFRH_TNG50/SFH_TNG50_out_'+str(output_number)+'_gal_'+str(int(x))+'_lowmass.txt'
    #print filename_gal
    #print(filename_gal)
    if os.path.exists(filename_gal):
        snap,red,galid,bh_mass,bh_acc,L_bol_1,Lbol_2, L_xray_1, L_xray_2, mstar, mgas, sfr, metallicity = np.loadtxt(filename_gal,unpack=True)

    
    met_gal = float(metallicity[0])*0.02
    logUs = 0.8*float(metallicity[0])-3.58
    #print(sfr)
    av_sfr = np.mean(sfr)
    
    #get key galaxy values from the library file
    age = df['age end SFR'][df['galaxy number'] == x].values[0]
    mstar = df['Mstar end SFR'][df['galaxy number'] == x].values[0]
    bh_mass = df['BH mass'][df['galaxy number'] == x].values[0]
    bh_lum = df['BH Lum'][df['galaxy number'] == x].values[0]
    #bh_mass = df_mod['BH mass'][df['galaxy number'] == x].values[0]
    #bh_lum = df_mod['BH Lum'][df['galaxy number'] == x].values[0]

    #bh_mass = 10*bh_mass
    #transform .ised file from this galaxy to ascii
    make_ascii(output_number, x)

    #create modified .ascii file which adds spectrum of AGN
    add_spectrum(x,age,bh_mass,bh_lum,met_gal,av_sfr,xid,nh,logUs,alpha)

    #convert modified .ascii file to .ised for further use
    make_ised(output_number, x)

    #name of modified file
    cspname = 'norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_lowmass_evenmoredust_BH_elines'

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

    #convert age to Gyr as requested in input
    age_today = age_today*1e-09

    #define output of modified magnitude file
    outname = 'norm_catalog/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_lowmass_evenmoredust_BH_elines'
    get_mags(age_today, filter_list, isedname,outname)
