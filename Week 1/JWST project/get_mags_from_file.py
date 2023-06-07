"""
Author: Clara Lilje

"""

import numpy as np
import os
import csv
import pandas as pd

#Define parameters of workspace
work_dir = 'JWST project/TNG_50_library/'

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
df = pd.read_csv(work_dir+'library_galaxies.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])

#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    outname = work_dir+'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'.magnitudes'
    df = pd.read_csv(outname,header=None,skiprows=41)
    df = df[0].str.split(expand=True)
    df = df.T.astype(float)

    headers = df.iloc[0].values
    df.columns = headers
    #df.drop(index=0, axis=0, inplace=True)
    to_print = df[5]
    to_print =to_print.to_numpy()
    to_print = np.insert(to_print,0,x)
    with open(work_dir+'magnitudes_galaxies.txt', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(to_print)



"""
This section covers only the galaxies with BH and considers their SEDs with the AGN added.
"""
#Read in the galaxy specific data from the galaxy library that was previously created
df = pd.read_csv(work_dir+'library_galaxies_with_BH.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])
#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    outname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_mod.magnitudes'
    df = pd.read_csv(outname,header=None,skiprows=41)
    df = df[0].str.split(expand=True)
    df = df.T.astype(float)    
    headers = df.iloc[0].values
    df.columns = headers
    #df.drop(index=0, axis=0, inplace=True)
    to_print = df[5]
    to_print =to_print.to_numpy()
    to_print = np.insert(to_print,0,x)
    with open(work_dir+'magnitudes_galaxies_with_BH.txt', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(to_print)