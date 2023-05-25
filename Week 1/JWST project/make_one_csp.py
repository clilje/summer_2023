import os
import pygalaxev
from pygalaxev_cosmology import uniage
import pandas as pd
import bisect
import numpy as np

# create dictionary of redshifts and universe age for TNG-50

#using parameters:
h0 = 67.74
omega_m = 0.3089
omega_lambda = 0.6911

redshifts = [5.0,5.23,5.53,5.85,6.01,6.49,7.01,7.24,7.6,8.01,8.45,9.0,9.39,10.0,10.98,11.98,14.99,20.05]
age_uni = [1.173e+9,1.109e+9,1.033e+09,0.962e+09,0.929e+09,0.841e+09,0.761e+09,0.729e+09,0.684e+9,0.638e+09,0.593e+09,0.545e+09,0.515e+9,0.472e+9,0.415e+9,0.368e+9,0.269e+9,0.178e+9]
age = dict(zip(redshifts,age_uni))

# creates one CSP model using galaxev

# selects the stellar template library:
output_number = 17
galaxy_number = 2
# Low-resolution 'BaSeL' library, Chabrier IMF
ssp_dir = '../Miles_Atlas/Chabrier_IMF/'
tempname = 'lr_BaSel'
work_dir = './'
ssp_dir = '../../'
filename = 'Files_SFRH_TNG50/SFH_TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)+'.txt'
df = pd.read_csv(filename,header=None,sep="\s+| |  ",
                       names=['output_number','redshift','galaxy ID',
                              'BH mass','BH accretion rate',
                              'bolLum1','bolLum2', 'lumXray1', 'lumXray2',
                              'Mstar','DMmass','SFR','Z'],engine='python')
print(df)
redshift = np.flip(df['redshift'].to_numpy())
sfr = np.flip(df['SFR'].to_numpy())
data_out = np.zeros((len(redshift),2))
birth = age[redshift[0]]
x = 0
while x < len(redshift):
    data_out[x][0] = age[redshift[x]]-birth
    data_out[x][1] = sfr[x]
    x += 1


df2 = pd.DataFrame(data_out)
    
df2.to_csv('TNG_50_library/SFRH_gal_'+str(galaxy_number)+'.ASCII',header=None,index=False,sep='\t')
sfh_pars = 'TNG_50_library/SFRH_gal_'+str(galaxy_number)+'.ASCII'
# Selects the metallicity:
a = df['Z'].to_numpy()
Z = a[a!=0][-1]

with open('TNG_50_library/library_galaxies.txt', 'a') as f:
    f.write(str(galaxy_number)+','+str(data_out[-1][0])+','+str(df['Mstar'].to_numpy()[0])+'\n')



#Z = 0.02 # solar metallicity
Zcode_dic = {0.0001:'m22', 0.0004:'m32', 0.004:'m42', 0.008:'m52', 0.02:'m62', 0.05:'m72', 0.1:'m82'}


res = bisect.bisect_left(list(Zcode_dic.keys()), Z)
print(res)
print(len(list(Zcode_dic.keys())))
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
print(Zcode)
print(Z)

tau = 1. # exponential star formation history timescale (in Gyr)
tau_V = 0.15 # dust effective optical depth
mu = 0.4 # fraction of attenuation due to diffuse interstellar medium
epsilon = 0. # gas recycling (no recycling if set to zero)

#isedname = ssp_dir+'bc2003_%s_%s_chab_ssp.ised'%(tempname, Zcode)
#isedname = 'bc2003_lr_BaSeL_m62_chab_ssp.ised'
isedname = 'bc2003_lr_BaSeL_'+str(Zcode)+'_chab_ssp.ised'
#outname = 'TNG_50_library/bc03_Z=%6.4f_tau=%5.3f_tV=%5.3f_mu=%3.1f_eps=%5.3f'%(Z, tau, tau_V, mu, epsilon)
outname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(galaxy_number)
pygalaxev.run_csp_galaxev(isedname, outname, sfh='custom', sfh_pars=sfh_pars, tau_V=tau_V,mu=mu, epsilon=0., work_dir=work_dir)

