The code in this folder was used to create galaxy and AGN spectral energy distributions for the paper by Habouzit, Lilje (in prep.) on the basis of BC03.
The files were kept separately to allow for modular use of the pipeline to understand the impact of the different components on magnitude or color values. 
Below is a more detailed list of instructions regarding the files. 
BC03 is required as a standard installation for this code to run.

# Pipeline

1. Get catalogues from TNG-50
    1. Galaxy stellar mass
    2. Galaxy star formation rate
    3. Galaxy metallicity
    4. BH mass
    5. BH bolometric luminosity

2. Use file **make_one_csp.py**
    1. **Inputs:**
        - Output number from catalogue
        - Galaxy ID of all galaxies in THIS snapshot
        - Stellar template library to use
        - Star Formation Rate history file
        - Take catalogue
        - Extract age at which it\'s done star forming
        - count backwards to get age
        - create seperate .ascii file with age and sfr of galaxy
        - Metallicity of galaxy
        - Metallicity at end of SFRH from TNG-50
        - match to the metallicities available in BC03
        - choose closest option
        - Dust parameters (see Paper Habouzit, Lilje (in prep.))
            - tau_V
            - mu

    2. **Operation**
        - extract all info from catalogues
        - make SFRH file
        - choose metallicity
        - write to library files
        - use **csp_galaxev** from bc03
    3. **Outputs:**

         - Library with stellar and Bh mass, age at end of sfrh and BH
        bolometric luminosity for all galaxies
         - seperate library with same columns for galaxies including a BH
         - .ised file for galaxy
         - .mass file, which is equivalent to .4color file from BC03

3. **OPTIONAL: Use get_one_sed_from_csp.py**
    1. **Inputs:**
        - Library files, which contain galaxies from which to create SED
        - .ised file created earlier

    2. **Operation:**
        - Read library file
        - Extract age at end of SFRH
        - call galaxevpl function of BC03 for every galaxy in library
    
    3. **Output:**
    - .sed file for each galaxy in library at age, which galaxy has at end of SFRH
    - Units: \$L\_{\\odot}/\\AA\$ for wavelengths from \$91\\AA\$ to \$160 \\mu m\$

4. Use file **get_mags_from_sed.py**
    1. **Inputs:**
        - Library of galaxies
        - List of filters
        - Redshift at Output number
        - cosmological parameters: (see Paper Habouzit, Lilje (in prep.))
    	- .ised file created earlier

    2. **Operation:**
        - Reads library file
        - Calculates age of galaxy today from age in file and desired redshift
        - Call function mm_evolution from BC03 code with desired cosmology

    3. **Output:**
        - .magnitudes file, which is equivalent to .multi_mag_AB file from BC03
        - 
5. Use file **get_mags_from_sed_with_elines.py**
    1. **Inputs:**
        - Library of galaxies
        - List of filters
        - Redshift at Output number
        - cosmological parameters: (see Paper Habouzit, Lilje (in prep.))
    	- .ised file created earlier

    2. **Operation:**
        - Reads library file
        - Transform .ised file of galaxy to .ascii_ised file using ascii_ised from BC03
        - split .ascii_ised file into sections
        - add narrow line emission lines from grid of models by Feltre, Charlot & Gutkin (2016) to one section
        - merge all sections together and re-write to file
        - transform modified .ascii_ised file to a .ised file using ised_ascii from BC03
        - calculate age of galaxy today
        - Call function mm_evolution from BC03 code with desired cosmology

    3. **Output:**
        - .magnitudes file, which is equivalent to .multi_mag_AB file from BC03

6. Use **get_mags_from_sed_with_BH.py**
    1. **Inputs:**
        - Library of galaxies
        -  List of filters
        - Redshift at Output number
        - cosmological parameters:  (see Paper Habouzit, Lilje (in prep.))
        - Parameters for SED model
        - .ised file created earlier

    2. **Operation:**
        - Read in library file
        - Transform .ised file of galaxy to .ascii_ised file using ascii_ised from BC03
        - add the spectrum of the AGN
        - split .ascii_ised file into sections
        - create spectrum for AGN with Volonteri+17 method and BH parameters over required wavelengths
        - add spectrum to one section
        - merge all sections together and re-write to file
        - transform modified .ascii_ised file to a .ised file using ised_ascii from BC03
        - calculate age of galaxy today
        - call mm_evolution function from BC03 with new modified .ised file and the age of the galaxy today

    3. **Outputs:**
        - .magnitudes file, which is equivalent to .multi_mag_AB file from BC03

7. Use file **get_mags_AGN_with_elines.py**
    1. **Inputs:**
        - Library of galaxies
        - List of filters
        - Redshift at Output number
        - cosmological parameters: (see Paper Habouzit, Lilje (in prep.))
    	- .ised file created earlier

    2. **Operation:**
        - Reads library file
        - Transform .ised file of galaxy continuum and emission lines and AGN continuum  to .ascii_ised file using ascii_ised from BC03
        - split .ascii_ised file into sections
        - add narrow line emission lines from grid of models by Feltre, Charlot & Gutkin (2016) to one section
        - merge all sections together and re-write to file
        - transform modified .ascii_ised file to a .ised file using ised_ascii from BC03
        - calculate age of galaxy today
        - Call function mm_evolution from BC03 code with desired cosmology

    3. **Output:**
        - .magnitudes file, which is equivalent to .multi_mag_AB file from BC03
    
8. **Use file get_mag_from_csp.py**

    1. **Inputs:**
        - Libraries containing all galaxies
        - Magnitude files created earlier

    2. **Operations:**
        - Read galaxy libraries
        - Read magnitude file
        - Find row at correct redshift
        - Write magnitudes to new file with galaxy ID

    3. **Output:**
        - Magnitude catalogues for all galaxies

