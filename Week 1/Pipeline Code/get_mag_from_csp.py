"""
Author: Clara Lilje

"""

import numpy as np
import os
import csv
import pandas as pd

#Define parameters of workspace
work_dir = 'norm_catalog/'

#This is catalogue specific file naming convention
output_number = 17
if output_number == 17:
    redshift = 5
if output_number == 11:
   redshift = 7
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

df = pd.read_csv(work_dir+'library_galaxies'+str(output_number)+'_lowmass.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])

#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    outname = work_dir+'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_lowmass_evenmoredust.magnitudes'
    df = pd.read_csv(outname,header=None,skiprows=41)
    df = df[0].str.split(expand=True)
    df = df.T.astype(float)

    headers = df.iloc[0].values
    df.columns = headers
    #df.drop(index=0, axis=0, inplace=True)
    to_print = df[redshift]
    to_print =to_print.to_numpy()
    to_print = np.insert(to_print,0,x)
    with open(work_dir+'magnitudes_galaxies'+str(output_number)+'_lowmass_stellarcontinuum_evenmoredust.txt', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(to_print)



"""
This section covers only the galaxies with BH and considers their SEDs with the AGN added.
"""

#Read in the galaxy specific data from the galaxy library that was previously created
df = pd.read_csv(work_dir+'library_galaxies'+str(output_number)+'_lowmass.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])
#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    outname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_lowmass_evenmoredust_elines.magnitudes'
    df = pd.read_csv(outname,header=None,skiprows=41)
    df = df[0].str.split(expand=True)
    df = df.T.astype(float)    
    headers = df.iloc[0].values
    df.columns = headers
    #df.drop(index=0, axis=0, inplace=True)
    to_print = df[redshift]
    to_print =to_print.to_numpy()
    to_print = np.insert(to_print,0,x)
    with open(work_dir+'magnitudes_galaxies'+str(output_number)+'_lowmass_fullgalaxy_evenmoredust.txt', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(to_print)



"""
This section covers all galaxies without taking their AGN into consideration.
"""

#Read in the galaxy specific data from the galaxy library that was previously created

df = pd.read_csv(work_dir+'library_galaxies_with_BH'+str(output_number)+'_lowmass.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])

#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    outname = work_dir+'TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_lowmass_evenmoredust_BH_cont.magnitudes'
    df = pd.read_csv(outname,header=None,skiprows=41)
    df = df[0].str.split(expand=True)
    df = df.T.astype(float)

    headers = df.iloc[0].values
    df.columns = headers
    #df.drop(index=0, axis=0, inplace=True)
    to_print = df[redshift]
    to_print =to_print.to_numpy()
    to_print = np.insert(to_print,0,x)
    with open(work_dir+'magnitudes_galaxies'+str(output_number)+'_lowmass_withTNG50AGN_continuum_evenmoredust.txt', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(to_print)



"""
This section covers only the galaxies with BH and considers their SEDs with the AGN added.
"""

#Read in the galaxy specific data from the galaxy library that was previously created
df = pd.read_csv(work_dir+'library_galaxies_with_BH'+str(output_number)+'_lowmass.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])
#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    outname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_lowmass_evenmoredust_BH_elines.magnitudes'
    df = pd.read_csv(outname,header=None,skiprows=41)
    df = df[0].str.split(expand=True)
    df = df.T.astype(float)    
    headers = df.iloc[0].values
    df.columns = headers
    #df.drop(index=0, axis=0, inplace=True)
    to_print = df[redshift]
    to_print =to_print.to_numpy()
    to_print = np.insert(to_print,0,x)
    with open(work_dir+'magnitudes_galaxies_with_BH'+str(output_number)+'_lowmass_withTNG50AGN_fullAGN_evenmoredust.txt', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(to_print)



"""
This section covers only the galaxies with BH and considers their SEDs with the AGN added.
"""
"""
#Read in the galaxy specific data from the galaxy library that was previously created
df = pd.read_csv(work_dir+'library_galaxies'+str(output_number)+'.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])
#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    outname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_elines_nodust.magnitudes'
    df = pd.read_csv(outname,header=None,skiprows=41)
    df = df[0].str.split(expand=True)
    df = df.T.astype(float)    
    headers = df.iloc[0].values
    df.columns = headers
    #df.drop(index=0, axis=0, inplace=True)
    to_print = df[redshift]
    to_print =to_print.to_numpy()
    to_print = np.insert(to_print,0,x)
    with open(work_dir+'magnitudes_galaxies_with_elines'+str(output_number)+'_nodust.txt', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(to_print)
"""



"""
This section covers only the galaxies with BH and considers their SEDs with the AGN added.
"""
"""
#Read in the galaxy specific data from the galaxy library that was previously created
df = pd.read_csv(work_dir+'library_galaxies'+str(output_number)+'.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])
#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    outname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_everything_moredust.magnitudes'
    df = pd.read_csv(outname,header=None,skiprows=41)
    df = df[0].str.split(expand=True)
    df = df.T.astype(float)    
    headers = df.iloc[0].values
    df.columns = headers
    #df.drop(index=0, axis=0, inplace=True)
    to_print = df[redshift]
    to_print =to_print.to_numpy()
    to_print = np.insert(to_print,0,x)
    with open(work_dir+'magnitudes_galaxies_with_everything'+str(output_number)+'_moredust.txt', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(to_print)

"""
"""
This section covers only the galaxies with BH and considers their SEDs with the AGN added.
"""
"""
#Read in the galaxy specific data from the galaxy library that was previously created
df = pd.read_csv(work_dir+'library_galaxies'+str(output_number)+'.txt',header=None,names=['galaxy number','age end SFR','Mstar end SFR','bh Mass','bh lum'])
#loop through different galaxies, indexing files by galaxy number
for x in df['galaxy number'].to_numpy():
    outname = work_dir+'/TNG50_out_'+str(output_number)+'_gal_'+str(x)+'_slopes_nodust.magnitudes'
    df = pd.read_csv(outname,header=None,skiprows=41)
    df = df[0].str.split(expand=True)
    df = df.T.astype(float)    
    headers = df.iloc[0].values
    df.columns = headers
    #df.drop(index=0, axis=0, inplace=True)
    to_print = df[redshift]
    to_print =to_print.to_numpy()
    to_print = np.insert(to_print,0,x)
    with open(work_dir+'magnitudes_galaxies_with_slopes'+str(output_number)+'_nodust.txt', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(to_print)
"""
