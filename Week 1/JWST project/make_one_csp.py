"""
Author: Clara Lilje, adapted from pygalaxev code for personal use.

This file will use the csp function from the BC03 galaxev code to create
a convolution of the metallicities and star formation rate histories from
the TNG-50 simulation. The output is a .ised file including the SEDs at 
multiple different ages of the galaxy and a .mass file which includes 
galaxy and stellar masses at different ages.

"""

import os
import pygalaxev
from pygalaxev_cosmology import uniage
import pandas as pd
import bisect
import numpy as np

# create dictionary of redshifts and universe age for TNG-50 (manually, using the UCLA Ned Wright calculator)
#using cosmological parameters (taken from TNG-50):
h0 = 67.74
omega_m = 0.3089
omega_lambda = 0.6911

redshifts = [5.0,5.23,5.53,5.85,6.01,6.49,7.01,7.24,7.6,8.01,8.45,9.0,9.39,10.0,10.98,11.98,14.99,20.05]
age_uni = [1.173e+9,1.109e+9,1.033e+09,0.962e+09,0.929e+09,0.841e+09,0.761e+09,0.729e+09,0.684e+9,0.638e+09,0.593e+09,0.545e+09,0.515e+9,0.472e+9,0.415e+9,0.368e+9,0.269e+9,0.178e+9]
age = dict(zip(redshifts,age_uni))



# creates one CSP model using galaxev

#These are catalogue specific file naming conventions
output_number = 17

#This decides which galaxies this get applied to, should probably be extended by a loop, checking which files exist
galaxy_number = [0,1,2,3,4,5]

# selects the stellar template library:
# Low-resolution 'Miles' library, Chabrier IMF
ssp_dir = '../Miles_Atlas/Chabrier_IMF/'        #I don't actually use these here
tempname = 'hr_xmiless'

#This just defines the work directory, this does get used
work_dir = './'
ssp_dir = '../../'

#Looping over all galaxies by index as defined above
for i in galaxy_number:

    #file naming convention as in my directories
    filename = 'Files_SFRH_TNG50/SFH_TNG50_out_'+str(output_number)+'_gal_'+str(i)+'.txt'

    #file convention as in the SFH files by Melanie
    df = pd.read_csv(filename,header=None,sep="\s+| |  ",
                        names=['output_number','redshift','galaxy ID',
                                'BH mass','BH accretion rate',
                                'bolLum1','bolLum2', 'lumXray1', 'lumXray2',
                                'Mstar','DMmass','SFR','Z'],engine='python')
    

    """
    Prepare individual SFRH for each galaxy
    """
    #get the important values for making the ASCII file as input for the BC03 code
    redshift = np.flip(df['redshift'].to_numpy())
    sfr = np.flip(df['SFR'].to_numpy())
    data_out = np.zeros((len(redshift),2))
    #defines birth date as age of the universe at the earliest redshift the galaxy appears
    birth = age[redshift[0]]
    x = 0
    while x < len(redshift):
        #prepares array to have age of galaxy and then sfr at that age
        data_out[x][0] = age[redshift[x]]-birth
        data_out[x][1] = sfr[x]
        x += 1

    #convert array to be data frame and then write to a file ready to use for BC03
    df2 = pd.DataFrame(data_out)
    df2.to_csv('TNG_50_library/SFRH_gal_'+str(i)+'.ASCII',header=None,index=False,sep='\t')
    sfh_pars = 'TNG_50_library/SFRH_gal_'+str(i)+'.ASCII'
    

    """
    Select metallicity of initial model
    """
    # Selects the metallicity:
    a = df['Z'].to_numpy()
    #we choose as metallicity the latest metallicity the galaxy has
    Z = a[a!=0][0]
    
    #Z = 0.02 # solar metallicity
    #Dictionary of the metallicities as they correspond within BC03
    Zcode_dic = {0.0001:'m22', 0.0004:'m32', 0.004:'m42', 0.008:'m52', 0.02:'m62', 0.05:'m72', 0.1:'m82'}

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
    #Zcode = Zcode_dic[Z]
    #print(Zcode)
    #print(Z)


    """
    Writes important values to library file
    """
    #This prepares the library for only galaxies that have a BH, to be used later
    if df['BH mass'].to_numpy()[0] != 0:
        with open('TNG_50_library/library_galaxies_with_BH.txt', 'a') as f:
            f.write(str(i)+','+str(data_out[-1][0])+','+str(df['Mstar'].to_numpy()[0])+','+str(df['BH mass'].to_numpy()[0])+','+str(df['bolLum2'].to_numpy()[0])+'\n')

    #write galaxy ID, age at end of SFR, mstar in msun, bhmass, AGN lum for all galaxies
    with open('TNG_50_library/library_galaxies.txt', 'a') as f:
        f.write(str(i)+','+str(data_out[-1][0])+','+str(df['Mstar'].to_numpy()[0])+','+str(df['BH mass'].to_numpy()[0])+','+str(df['bolLum2'].to_numpy()[0])+'\n')




    #These are the dust absorbtion parameters to be used in the csp model

    tau = 1. # exponential star formation history timescale (in Gyr)
    tau_V = 0.15 # dust effective optical depth
    mu = 0.4 # fraction of attenuation due to diffuse interstellar medium
    epsilon = 0. # gas recycling (no recycling if set to zero)



    #notes
    #isedname = ssp_dir+'bc2003_%s_%s_chab_ssp.ised'%(tempname, Zcode)
    #isedname = 'bc2003_lr_BaSeL_m62_chab_ssp.ised'
    #bc2003_hr_xmiless_m22_chab_ssp

    #The name of the original base file to be convolved
    isedname = 'bc2003_hr_xmiless_'+str(Zcode)+'_chab_ssp.ised'
    #outname = 'TNG_50_library/bc03_Z=%6.4f_tau=%5.3f_tV=%5.3f_mu=%3.1f_eps=%5.3f'%(Z, tau, tau_V, mu, epsilon)
    
    #Name of final output files
    outname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(i)

    #Run the pygalaxev method which returns an .ised and a .mass file
    pygalaxev.run_csp_galaxev(isedname, outname, sfh='custom', sfh_pars=sfh_pars, tau_V=tau_V,mu=mu, epsilon=0., work_dir=work_dir)

