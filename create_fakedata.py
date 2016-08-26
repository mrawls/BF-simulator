#from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from PyAstronomy import pyasl
'''
August 2014
by Meredith Rawls
** updated August 2016 **

Applies a binary star RV model to a pair of input (presumably model) spectra
to generate synthetic spectral observations of an SB2 binary.
--> useful to test spectral decomposition or see how well the BF should be working!

There is an option to make a quick plot of the fake SB2 spectra.
Lots of input files and parameters are hardwired below.

Main output file: formatted for use by FDBinary/fd3
    It will say '# (number of spectra+1) X (length of each spectrum)' on the 1st line
    The 1st column is natural-log wavelength from ln(wavemin) to ln(wavemax)
    Subsequent columns are each synthetic SB2 spectrum

Individual output files: one two-column file per spectrum, wavelength and flux.
'''
# A pair of input spectra - should be normalized two-column text files
# These are the theoretical individual spectra of each component of the binary
inspecA = 'lte04800-2.50-0.0.PHOENIX-4400-6000-norm.txt'
inspecB = 'lte05500-2.50-0.0.PHOENIX-4400-6000-norm.txt'
vsiniA = 0
vsiniB = 20
tripleSysFlag = 1 # set to 1 if star A should always have zero RV
                  # this would apply if there is a faint 3rd component causing star 2's RVs

# Light ratios to be applied to stars A and B (should sum to 1)
lrA = 0.5
lrB = 0.5

# This is a file with timestamps of observations in the first column
timefile = 'obstimes.txt'
times = np.loadtxt(timefile, comments='#', usecols=(0,), unpack=True)

# OPTION 1: Use these options to apply an RV curve from REAL DATA
# (to use OPTION 2 instead, don't set rvfile)
#rvfile = '../../KIC_8848288/rvs_revisited3_BF.txt'

# OPTION 2: Use these options to apply an RV curve from a MODEL
# following values are for KIC 8848288 (TYC 3559)
porb = 5.56648
t0 = 2454904.8038
ecc = 0
argper = 0 #degrees!
amp2 = 6. # km/s, star B
amp1 = 144. # km/s, star A (for triple system when star A is constant, actually star C)
            # currently set assuming (mass of star B) / (mass of BD) = 24
            # K1/K2 = M2/M1, so K1 = K2*(M2/M1)

# This is the file that will save the final spectral data (formatted for use with FDBinary/fd3)
FDBinary_file = 'fd3_infile_fakebinary.obs'

# Boundaries of the new wavelength scale that will be created
# It will be evenly spaced in ln-wavelength because fd3 was written by people who like ln
wavemin = 4500 #3900.
wavemax = 5900 #8450.
logwavemin = np.log10(wavemin)
lnwavemin = np.log(np.power(10,logwavemin))
logwavemax = np.log10(wavemax)
lnwavemax = np.log(np.power(10,logwavemax))
deltalogwave = 0.0000035
deltalnwave = np.log(np.power(10,deltalogwave))
lnwave = np.arange(lnwavemin, lnwavemax, deltalnwave)

# Read in the pair of inspecs (text files) for Star A and Star B
waveA, specA = np.loadtxt(inspecA, comments='#', usecols=(0,1), unpack=True)
waveB, specB = np.loadtxt(inspecB, comments='#', usecols=(0,1), unpack=True)

# Broaden the spectra if necessary
# (You need an evenly spaced wavelength grid for rotBroad to work)
linearLD = 1.0 # if vsini != 0, limb darkening affects the line broadening profiles
if vsiniA != 0:
    print('Broadening Star A spectrum...')
    evenwaveA = np.arange(np.min(waveA), np.max(waveA), np.max(waveA[-1]-waveA[-2], waveA[1]-waveA[0]))
    specA = np.interp(evenwaveA, waveA, specA)
    specA = pyasl.rotBroad(evenwaveA, specA, linearLD, vsiniA)
    waveA = evenwaveA
if vsiniB != 0:
    print('Broadening Star B spectrum...')
    evenwaveB = np.arange(np.min(waveB), np.max(waveB), np.max(waveB[-1]-waveB[-2], waveB[1]-waveB[0]))
    specB = np.interp(evenwaveB, waveB, specB)
    specB = pyasl.rotBroad(evenwaveB, specB, linearLD, vsiniB)
    waveB = evenwaveB

# Placeholder lists for RV shifts to be applied
deltav1s = []
deltav2s = []

try: # OPTION 1
    # read RVs straight from a file
    f1 = open(rvfile)

except: # OPTION 2
    # create a model RV curve
    # save the resulting velocity shifts in the deltav1s and deltav2s lists
    print('Creating model RV curve (option 2)')
    ks = pyasl.MarkleyKESolver()
    argper = np.radians(argper)
    for time in times:
        MA = 2*np.pi*(time - t0)/porb #mean anomaly
        EA = ks.getE(MA, ecc) #eccentric anomaly
        cosTA = ((np.cos(EA) - ecc) / (1 - ecc*np.cos(EA)))
        sinTA = (np.sqrt(1-(ecc*ecc))*np.sin(EA)) / (1 - ecc*np.cos(EA))
        TA = np.arctan2(sinTA, cosTA) #true anomaly
        deltav1s.append(-amp1*(np.cos(TA+argper) + ecc*np.cos(argper)))
        deltav2s.append(amp2*(np.cos(TA+argper) + ecc*np.cos(argper)))
#    plt.plot(times, deltav1s, color='0.75', marker='o', mec=red, mfc='r', ls=':', ms=7, mew=0)
#    plt.plot(times, deltav2s, color='0.75', marker='o', mec=yel, mfc='b', ls=':', ms=7, mew=0)
#    plt.ylabel('Model Radial Velocity (km s$^{-1}$)')
#    plt.xlabel('Time')
#    plt.show()

else: # OPTION 1 CONTINUED
    print('Using list of RV values (option 1)')
    f1.close()
    deltav1s, deltav2s = np.loadtxt(rvfile, comments='#', usecols=(3,5), unpack=True)

# Confirm velocity values look correct
if tripleSysFlag == 1: # set star 1's RVs to zero
    deltav1s = [0 for item in range(0, len(deltav1s))]
print('These are the velocities being applied to Star A (km/s): {0}'.format(deltav1s))
print('These are the velocities being applied to Star B (km/s): {0}'.format(deltav2s))

# Apply pairs of delta-vs (shifts) to each set of A and B.
waveAlist = []; specAlist = []
waveBlist = []; specBlist = []
c = 2.99792e5 # in km/s, just like deltav1s & deltav2s

for idx in range(0, len(times)):
    newwaveA = waveA * deltav1s[idx] / c + waveA
    waveAlist.append( newwaveA )    
    newspecA = np.interp(lnwave, np.log(newwaveA), specA)
    specAlist.append( newspecA )
    newwaveB = waveB * deltav2s[idx] / c + waveB
    waveBlist.append( newwaveB )
    newspecB = np.interp(lnwave, np.log(newwaveB), specB)
    specBlist.append( newspecB )

# Create the fake SB2 spectra by combining pairs of A and B.
# Make a quick plot to make sure nothing terrible happened and SB2 spectra exist.
SB2list = []
ploffset = 0
for i in range(0, len(times)):
    SB2list.append( lrA*specAlist[i] + lrB*specBlist[i] )
    plt.plot(np.power(np.exp(1),lnwave), SB2list[i] + ploffset)
    ploffset += 1.0
plt.show()

# Use each SB2 to make a text outfile formatted for use by FDBinary/fd3.
fout = open(FDBinary_file, 'w')
print('# ' + str(len(times)+1) + ' X ' + str(len(lnwave)), file=fout)
for j in range(0, len(lnwave)):
    newstring = str(lnwave[j])
    for i in range(0, len(times)):
        newstring += '\t' + str(SB2list[i][j])
    print(newstring, file=fout)
fout.close()

# Use each SB2 to also make a set of individual text outfiles for use with BF-python.py
for idx, SB2 in enumerate(SB2list):
    waves = np.power(np.exp(1),lnwave)
    onespecfile = open('fakedata'+str(idx)+'.txt', 'w')
    for wave, SB in zip(waves, SB2):
        print(wave, SB, file=onespecfile)
    onespecfile.close()