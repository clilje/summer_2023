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

# this code reads in a CSP SED file and, given a galaxy redshift, calculates the magnitudes in a series of filters

work_dir = 'TNG_50_library/'
#pygalaxevdir = os.environ.get('PYGALAXEVDIR')
#filtdir = pygalaxevdir+'/filters/'

#redshift = 5
#Dlum = pygalaxev_cosmology.Dlum(redshift) # luminosity distance in Mpc

#log_mstar = 11.
output_number = 17
galaxy_number = [0,1,2,3]
tau = 1.
tau_V = 0.15
mu = 0.4
epsilon = 0.
#cosmological parameters
h0 = 67.74
omega_m = 0.3089
omega_lambda = 0.6911

def get_mags(age_today, filters, filename,outname, input_tmpname='tmp.in', output_tmpname='mySSP'):
    
    tmpname = work_dir+'/%s'%input_tmpname

    inputlines = []
    # input cosmology
    inputlines.append(str(h0)+','+str(omega_m)+','+str(omega_lambda)+'\n')
    inputlines.append(str(age_today)+'\n')
    inputlines.append('%s\n'%filename)

    filterstring =','.join([str(x) for x in filters])
    inputlines.append(filterstring+'\n')

    f = open(tmpname, 'w')
    f.writelines(inputlines)
    f.close()

    # Run bc03
    os.system('$bc03/mm_evolution < %s'%tmpname)
    #os.system('csp < %s'%tmpname)
    #os.system('/src/csp_galaxev < %s'%tmpname)
    # Keep the output *.ised and *.4color files
    os.system('scp %s.multi_mag_AB %s.magnitudes'%(outname, outname))
    os.system('rm -f %s.multi_mag_vega'%( outname))
    os.system('rm -f %s.multi_mag_AB'%( outname))
    # Clean up
    os.system('rm -f %s/%s*'%(work_dir, output_tmpname))


### Parameters for the AGN SED
nu_IR=0.01*3.2899*1e15 #0.01 Ryd to Hz
alpha_UV=-0.5
alpha_x=-1.0

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


def S(L_BH, M_BH):
  # M_BH in Msun
  # L_BH/L_Edd the Eddington ratio 
  # Here our L is in real units.
  #L_Edd = M_BH*4.*np.pi*cte_G*cte_m_p*cte_c/cte_sigma_t
  L_Edd = 1.26*1e38 * M_BH # in erg/s
  S     = np.log10((L_BH/L_Edd)/M_BH)
  #     = np.log10((L_BH)/M_BH)
  return S

def R(r_corona, M_BH):
  # r_corona the coronal radius (inner edge of the visible portion of the disk). Should be in units of r_g
  # r_g the gravitational radius, r_g = cte_G*M_BH/cte_c**2
  #r_g   = cte_G*M_BH/cte_c**2
  #R     = np.log10(r_corona/r_g)
  R = np.log10(r_corona)
  return R

def f1(L_BH, M_BH, r_corona):
  f1    = Ax2*R(r_corona, M_BH)*R(r_corona, M_BH) + Ax1*R(r_corona, M_BH) + Ax0 + B*S(L_BH,M_BH)
  return f1

def f2(L_BH, M_BH, r_corona):
  f2    = (b3x1*R(r_corona, M_BH)+b3x0)*np.power(S(L_BH,M_BH)+6,3) + (b1x2*R(r_corona, M_BH)*R(r_corona, M_BH) + b1x1*R(r_corona, M_BH) + b1x0)*(S(L_BH, M_BH)+6) + (ax2*R(r_corona, M_BH)*R(r_corona, M_BH) + ax1*R(r_corona, M_BH) + ax0)
  return f2

def log10_E_peak(L_BH,M_BH,r_corona):
  # E_peak used to parametrize both the location in energy space and the shape of the Big Blue Bump of the disk emission.
  # E_peak in keV according to Thomas
  #rint('ours',f1(L_BH, M_BH, r_corona),f2(L_BH, M_BH, r_corona))
  log10_E_peak = max(f1(L_BH, M_BH, r_corona),f2(L_BH, M_BH, r_corona))
  return log10_E_peak


def alpha_ox(L_BH,M_BH):
  L_Edd     = 1.26*10**38*M_BH # in erg/s
  alpha_ox  =  -0.13*np.log(M_BH)+0.15*np.log10(L_BH/L_Edd)-0.45
  return alpha_ox

def F_nu(nu,L_BH,M_BH,r_corona):
  #This line is meant to convert from E_peak in kev to nu_peak in Hz
  numax_peak = 2.417e17*10**log10_E_peak(L_BH,M_BH,r_corona)    # in Hz 

  nu1=4.8359e17 # in Hz?
  nu2=1.1992e15 # in Hz?

  param = (nu1**(alpha_UV) * np.exp(-nu1/numax_peak) * np.exp(-nu_IR/nu1)-
           10**(alpha_ox(L_BH,M_BH)*2.6056) * nu2**alpha_UV * np.exp(-nu2/numax_peak)*
           exp(-nu_IR/nu2))/(10**(alpha_ox(L_BH,M_BH)*2.6056) * nu2**
                             alpha_x - nu1**alpha_x)
  if param < 0:
    param = (10**(alpha_ox(L_BH,M_BH)*2.6056) * nu2**alpha_UV * exp(-nu2/numax_peak)*
             exp(-nu_IR /nu2))/(nu1**alpha_x)


  #PARAM = 0. ## set to 0 for now.
  PARAM = param
  #rint('numax_peak=',numax_peak)
  #print('factor Fnu,',np.exp(-nu / numax_peak), np.exp(-nu_IR / nu) )
  
  F_nu  = nu**alpha_UV * np.exp(-nu / numax_peak) * np.exp(-nu_IR / nu)  + PARAM*nu**alpha_x
  ##F_nu[np.where(nu<(0.1*3.2899e15))[0]] = nu[np.where(nu<(0.1*3.2899e15))[0]]**alpha_UV * np.exp(-nu[np.where(nu<(0.1*3.2899e15))[0]] / numax_peak) * np.exp(-nu_IR / nu[np.where(nu<(0.1*3.2899e15))[0]])
  for ii in range(len(nu)):
    if nu[ii]<3.2899e14: ### last term of Eq set to zero below 1.36 eV = 912 nm
      #print('nu= ',nu[ii],3.2899e14)
      F_nu[ii]  = nu[ii]**alpha_UV * np.exp(-nu[ii] / numax_peak) * np.exp(-nu_IR / nu[ii])
  return F_nu

def F_lambda(F_nu,nu):
  # from F_nu to F_lambda
  #from (erg/s/cm^2/Hz) to (erg/s/cm^2/A)
  F_lambda = F_nu*(3e18)/((3e18/nu)**2)
  return F_lambda


def make_spectrum_AGN(wavelength,bh_mass,bh_lum):
   
    L_BH  = 10**bh_lum # in erg/s
    M_BH  = bh_mass # in Msun
    L_Edd = 1.26*10**38*M_BH # in erg/s
    f_Edd_1 = np.log10(L_BH/L_Edd)
    
    nu = np.divide(3e18,wavelength) # from AA to Hz
    r_corona = 10.
    to_print_F_nu=F_nu(nu,L_BH,M_BH,r_corona)
    nu_peak = 10**(log10_E_peak(L_BH,M_BH,r_corona))*2.42e17
    ## normalization
    norm_1=(L_BH)/scipy.integrate.quadrature(F_nu,nu_peak*1e-3,nu_peak*1e3,args=(L_BH,M_BH,r_corona))[0]
    return(F_lambda(norm_1*to_print_F_nu, nu))




def add_spectrum(i, age_gal_i,bh_mass,bh_lum):
    age_row = [0]
    ages = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows = lambda x: x not in age_row)
    ages = ages[0].str.split(expand=True, n=2023)
    ages = ages.drop(0,axis=1)
    ages = np.float64(ages.to_numpy()[0])
    

    index_row = bisection(ages, age_gal_i)
    rows_age = [index_row+5, index_row+6, index_row+7]

    headerplustail_rows = [0,1,2,3,4,5,6,281,280,279,278,277,276,275,274,273,272,271,270]
    #df = pd.read_csv('TNG_50_library\TNG50_out_17_gal_'+str(i)+'.ised_ASCII','/s', header=None, skiprows=headerplustail_rows)
    df = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII', header=None, skiprows=lambda x: x not in rows_age)
    df = df[0].str.split(expand=True, n=2024)
    #print(df.loc[[0]])
    print(rows_age)

    first_col = df[0]
    last_col = df[2024]
    df = df.drop(0, axis=1)
    df = df.drop(2024, axis=1)

    df = df.astype(float)
    df = df.T

    
    specific_rows = [6]
    wavelength = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII', header=None, skiprows = lambda x: x not in specific_rows)
    wavelength = wavelength[0].str.split(expand=True, n=2023)
    wavelength = wavelength.drop(0, axis=1)
    wavelength = wavelength.T
    xvals = wavelength[0].to_numpy()
    xvals = np.float64(xvals)
    
    AGN = np.array(make_spectrum_AGN(xvals,bh_mass,bh_lum))
    print(AGN)
    print(df)
    df[3] = AGN
    df[0] = df.loc[:,[0,3]].sum(axis=1)
    df[1] = df.loc[:,[1,3]].sum(axis=1)
    df[2] = df.loc[:,[2,3]].sum(axis=1)
    df = df.drop(3, axis=1)
    print(df)
    df = df.T
    
    df.insert(0, 0, first_col, True)
    df.insert(2024, 2024, last_col, True)
    middle = df.astype(str)
    middle = middle.apply(lambda x: '       '.join(x.dropna()), axis=1)

    skipfooter = 282-rows_age[0] -(262-int(first_col[0]))
    front = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skipfooter= 282-rows_age[0])

    back = pd.read_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'.ised_ASCII',delimiter='/s', header=None, skiprows=rows_age[-1]+1)

    frames = [front,middle,back]
    result = pd.concat(frames)
    result.to_csv('TNG_50_library/TNG50_out_17_gal_'+str(i)+'_mod.ised_ASCII',header=None,index=None)

def bisection(array,value):
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
    print('hi')
    cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'.ised'
    os.system('$bc03/ascii_ised TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'.ised')


def make_ised(output_number, galaxy_number):
    print('hey')
    cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_mod.ised_ASCII'
    os.system('$bc03/bin_ised TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'_mod.ised_ASCII')

filter_list = [278,279,280,281,282,283,286,291]



if output_number == 17:
    redshift = 5


df = pd.read_csv(work_dir+'library_galaxies_with_BH.txt',header=None,
                    names=['galaxy number','age end SFR','Mstar end SFR','BH mass','BH Lum'])
galaxy_number = df['galaxy number'].to_numpy()
galaxy_number = [0]
for x in galaxy_number:

    #age = 11. # time since beginning of star formation (in Gyr)
    age = df['age end SFR'][df['galaxy number'] == x].values[0]
    mstar = df['Mstar end SFR'][df['galaxy number'] == x].values[0]
    bh_mass = df['BH mass'][df['galaxy number'] == x].values[0]
    bh_lum = df['BH Lum'][df['galaxy number'] == x].values[0]
    make_ascii(output_number, x)
    add_spectrum(x,age,bh_mass,bh_lum)
    make_ised(output_number, x)

    cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_mod'

    isedname = str(cspname)+'.ised'
    print(age)
    if redshift == 5:
        #age_today = (13.798e+09-1.173e+09) + age
        age_today = 12.636e+09 + age
    age_today = age_today*1e-09


    outname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_mod'
    get_mags(age_today, filter_list, isedname,outname)
