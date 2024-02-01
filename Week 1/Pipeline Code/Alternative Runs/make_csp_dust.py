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



#These are the dust absorbtion parameters to be used in the csp model

#tau = 1. # exponential star formation history timescale (in Gyr)
tau_V = 0.15 # dust effective optical depth
mu = 0.4 # fraction of attenuation due to diffuse interstellar medium
epsilon = 0. # gas recycling (no recycling if set to zero)


tau_V = [0.1,0.5,1,2,3,4,5,6]
mu = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 ]
# creates one CSP model using galaxev

#These are catalogue specific file naming conventions
output_number = 11
galaxy_number = [384]
#This decides which galaxies this get applied to, should probably be extended by a loop, checking which files exist
#for output 7
#galaxy_number = [0, 1, 384, 385, 751, 924, 1180, 1433, 1727, 1985, 2314, 2544, 2794, 2904, 3017, 3135, 3325, 3543, 3718, 3924, 4123, 4226, 4613, 4773, 4911, 5025, 5164, 5305, 5444, 5595, 5765, 5891, 6108, 6237, 6335, 6457, 6729, 6992, 7096, 7524, 7611, 7688, 7824, 8268, 8387, 8461, 8720, 9167, 9219, 9273, 9811, 9909, 10105, 10202, 10290, 10385, 10479, 10582, 10690, 10770, 10956, 11015, 11173, 11449, 11569, 12369, 12505, 13844]
"""
# output number 17, reduced
galaxy_number = [0,
 1,
 2,
 3,
 1228,
 2102,
 2103,
 3119,
 3120,
 4232,
 4906,
 4908,
 5663,
 6099,
 6712,
 6713,
 7371,
 7818,
 8136,
 8705,
 9094,
 9600,
 10049,
 10293,
 10294,
 10653,
 10871,
 11124,
 11452,
 11703,
 11873,
 12234,
 12396,
 12561,
 12791,
 13037,
 13613,
 13867,
 14025,
 14236,
 14524,
 14808,
 15025,
 15211,
 15483,
 15716,
 15936,
 16043,
 16281,
 16493,
 16798,
 17006,
 17298,
 17455,
 17856,
 18062,
 18366,
 18812,
 19362,
 19513,
 19668,
 19831,
 20479,
 20684,
 20809,
 20949,
 21127,
 21755,
 21873,
 22093,
 22205,
 22536,
 22679,
 22989,
 23437,
 23634,
 23742,
 23960,
 24120,
 24273,
 24568,
 24708,
 25069,
 25849,
 26142,
 26267,
 26680,
 26806,
 26945,
 27628,
 27888,
 28110,
 28418,
 28802,
 28918,
 29562,
 29672,
 29780,
 30532,
 30933,
 31501,
 32884,
 34320,
 34852,
 36269,
 38051,
 38580,
 38661,
 39737,
 40323,
 41984,
 42817,
 42995,
 43246,
 43748,
 51529]
"""
"""
#output number 17, full stack
galaxy_number = [0, 1, 2, 3, 1228, 1229, 2102, 2103, 3119, 3120, 3121, 3123, 4232, 4233, 4906, 4908, 5663, 6099, 6712, 6713, 6714,
                  7371, 7374, 7818, 8136, 8137, 8141, 8705, 8708, 9094, 9095, 9600, 9601, 10049, 10293, 10294, 10295, 10653, 10871,
                   10872, 11124, 11125, 11452, 11454, 11703, 11873, 11874, 12234, 12396, 12397, 12561, 12791, 13037, 13290,
                     13613, 13867, 14025, 14236, 14524, 14808, 15025, 15211, 15483, 15716, 15936, 16043, 16044, 16281, 16493,
                       16798, 17006, 17298, 17455, 17600, 17856, 18062, 18366, 18522, 18812, 18815, 19070, 19362, 19513, 19668,
                         19831, 20049, 20235, 20236, 20479, 20684, 20809, 20949, 21127, 21548, 21549, 21755, 21873, 22093, 22205,
                           22358, 22536, 22679, 22793, 22989, 23195, 23437, 23438, 23634, 23742, 23960, 24120, 24273, 24401, 24402,
                             24568, 24708, 24884, 25069, 25070, 25225, 25375, 25480, 25691, 25849, 26142, 26267, 26680, 26806, 
                             26945, 27252, 27440, 27628, 27759, 27888, 28035, 28110, 28248, 28418, 28537, 28698, 28802, 28918,
                               29049, 29050, 29562, 29672, 29780, 29877, 30074, 30384, 30532, 30667, 30806, 30933, 31076, 31251, 31394,
                                 31501, 31633, 31784, 31931, 32039, 32192, 32351, 32500, 32590, 32884, 32991, 33172, 33266, 33387, 33555,
                                   33713, 34031, 34146, 34320, 34381, 34516, 34664, 34745, 34852, 34971, 35124, 35235, 35402, 35498, 35793,
                                     35935, 36050, 36121, 36269, 36366, 36487, 36606, 36718, 36824, 36926, 37085, 37201, 37339, 37504, 37621,
                                       37737, 37878, 38051, 38455, 38580, 38661, 38758, 38888, 38963, 39121, 39236, 39361, 39437,
 39536,
 39630,
 39737,
 39826,
 39939,
 40090,
 40194,
 40323,
 40383,
 40617,
 40808,
 40988,
 41098,
 41198,
 41312,
 41407,
 41547,
 41625,
 41693,
 41777,
 41856,
 41984,
 42052,
 42141,
 42250,
 42381,
 42491,
 42604,
 42732,
 42817,
 42896,
 42995,
 43049,
 43116,
 43246,
 43285,
 43419,
 43540,
 43633,
 43748,
 43814,
 43914,
 44040,
 44151,
 44277,
 44417,
 44503,
 44550,
 44613,
 44988,
 45099,
 45204,
 45317,
 45396,
 45485,
 45588,
 45688,
 45921,
 46004,
 46162,
 46238,
 46368,
 46454,
 46556,
 46712,
 46978,
 47050,
 47168,
 47293,
 47462,
 47581,
 47658,
 47763,
 47866,
 47951,
 48058,
 48170,
 48291,
 48389,
 48480,
 48572,
 48658,
 48735,
 48822,
 48916,
 48998,
 49068,
 49163,
 49373,
 49445,
 49527,
 49600,
 49688,
 49899,
 49997,
 50050,
 50290,
 50353,
 50400,
 50479,
 50567,
 50692,
 50859,
 51294,
 51395,
 51529,
 51682,
 51757,
 51955,
 52068,
 52275,
 52370,
 52454,
 52527,
 52603,
 52684,
 52755,
 52833,
 52938,
 53010,
 53099,
 53189,
 53388,
 53485,
 53590,
 53683,
 53759,
 53861,
 53924,
 54005,
 54084,
 54289,
 54360,
 54459,
 54582,
 54631,
 54694,
 54800,
 54887,
 54980,
 55234,
 55317,
 55476,
 55556,
 55650,
 55875,
 56215,
 56388,
 56468,
 56541,
 56616,
 56799,
 56889,
 57109,
 57196,
 57263,
 57353,
 57465,
 57547,
 57597,
 57651,
 57724,
 57806,
 57951,
 58379,
 58447,
 58530,
 58580,
 58639,
 58733,
 58820,
 58860,
 58941,
 58994,
 59095,
 59179,
 59296,
 59363,
 59427,
 59493,
 59580,
 59668,
 59748,
 59822,
 59924,
 59988,
 60233,
 60416,
 60502,
 60568,
 60742,
 60817,
 60891,
 60958,
 61020,
 61115,
 61295,
 61373,
 61440,
 61494,
 61571,
 61641,
 61725,
 61816,
 61865,
 61938,
 61986,
 62062,
 62208,
 62278,
 62446,
 62520,
 62596,
 62673,
 62834,
 62881,
 62948,
 62994,
 63130,
 63292,
 63445,
 63712,
 63986,
 64057,
 64118,
 64179,
 64261,
 64383,
 64440,
 64506,
 64679,
 64706,
 64772,
 65091,
 65173,
 65238,
 65319,
 65385,
 65432,
 65516,
 65597,
 65753,
 66012,
 66283,
 66340,
 66403,
 66748,
 66818,
 66889,
 66998,
 67162,
 67235,
 67379,
 67462,
 67534,
 67604,
 67661,
 67952,
 68082,
 68151,
 68287,
 68413,
 68484,
 68544,
 68733,
 68960,
 69013,
 69203,
 69364,
 69602,
 69871,
 70071,
 70358,
 70822,
 71248,
 71547,
 72546,
 73018,
 73931,
 75070,
 75666,
 77037,
 82571,
 93706]
"""
# selects the stellar template library:
# Low-resolution 'Miles' library, Chabrier IMF
ssp_dir = '../Miles_Atlas/Chabrier_IMF/'        #I don't actually use these here
tempname = 'hr_xmiless'

#This just defines the work directory, this does get used
work_dir = './'
ssp_dir = '../../'
for tt in tau_V:
  for mm in mu:
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
        df2.to_csv('TNG_50_library/'+str(output_number)+'SFRH_gal_'+str(i)+'.ASCII',header=None,index=False,sep='\t')
        """
        sfh_pars = 'TNG_50_library/'+str(output_number)+'SFRH_gal_'+str(i)+'.ASCII'
        

        """
        Select metallicity of initial model
        """
        # Selects the metallicity:
        a = df['Z'].to_numpy()
        #we choose as metallicity the latest metallicity the galaxy has
        Z = float(a[a!=0][0])*0.02
        #print(a,Z)
        
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
        """
        #This prepares the library for only galaxies that have a BH, to be used later
        if df['BH mass'].to_numpy()[0] != 0:
            with open('TNG_50_library/library_galaxies_with_BH'+str(output_number)+'.txt', 'a') as f:
                f.write(str(i)+','+str(data_out[-1][0])+','+str(df['Mstar'].to_numpy()[0])+','+str(df['BH mass'].to_numpy()[0])+','+str(df['bolLum2'].to_numpy()[0])+'\n')

        #write galaxy ID, age at end of SFR, mstar in msun, bhmass, AGN lum for all galaxies
        with open('TNG_50_library/library_galaxies'+str(output_number)+'.txt', 'a') as f:
            f.write(str(i)+','+str(data_out[-1][0])+','+str(df['Mstar'].to_numpy()[0])+','+str(df['BH mass'].to_numpy()[0])+','+str(df['bolLum2'].to_numpy()[0])+'\n')
      """




        #notes
        #isedname = ssp_dir+'bc2003_%s_%s_chab_ssp.ised'%(tempname, Zcode)
        #isedname = 'bc2003_lr_BaSeL_m62_chab_ssp.ised'
        #bc2003_hr_xmiless_m22_chab_ssp

        #The name of the original base file to be convolved
        isedname = 'bc2003_hr_xmiless_'+str(Zcode)+'_chab_ssp.ised'
        #outname = 'TNG_50_library/bc03_Z=%6.4f_tau=%5.3f_tV=%5.3f_mu=%3.1f_eps=%5.3f'%(Z, tau, tau_V, mu, epsilon)
        
        #Name of final output files
        outname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(i)+'tau='+str(tt)+'mu='+str(mm)

        #Run the pygalaxev method which returns an .ised and a .mass file
        pygalaxev.run_csp_galaxev(isedname, outname, sfh='custom', sfh_pars=sfh_pars, tau_V=tt,mu=mm, epsilon=0., work_dir=work_dir)

