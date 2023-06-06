import numpy as np
import h5py
import pygalaxev_cosmology
from pygalaxev_cosmology import L_Sun, Mpc, c as csol
from scipy.interpolate import splrep, splev, splint
import os
import pandas as pd


# this code reads in a CSP SED file and, given a galaxy redshift, calculates the magnitudes in a series of filters

work_dir = 'TNG_50_library/'
#pygalaxevdir = os.environ.get('PYGALAXEVDIR')
#filtdir = pygalaxevdir+'/filters/'

#redshift = 5
#Dlum = pygalaxev_cosmology.Dlum(redshift) # luminosity distance in Mpc

#log_mstar = 11.
output_number = 17
#galaxy_number = [0,1,2,3,4,5]
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





if output_number == 17:
    redshift = 5

df = pd.read_csv(work_dir+'library_galaxies.txt',header=None,
                    names=['galaxy number','age end SFR','Mstar end SFR','BH mass','BH Lum'])


    #age = 11. # time since beginning of star formatio
for x in df['galaxy number'].to_numpy():
    cspname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(x)
    #age = 11. # time since beginning of star formation (in Gyr)
    age = df['age end SFR'][df['galaxy number'] == x].values[0]
    mstar = df['Mstar end SFR'][df['galaxy number'] == x].values[0]
    isedname = str(cspname)+'.ised'
    print(age)
    if redshift == 5:
        #age_today = (13.798e+09-1.173e+09) + age
        age_today = 12.636e+09 + age
    age_today = age_today*1e-09
    filter_list = [278,279,280,281,282,283,286,291]
    outname = 'TNG_50_library/TNG50_out_'+str(output_number)+'_gal_'+str(x)
    get_mags(age_today, filter_list, isedname,outname)
