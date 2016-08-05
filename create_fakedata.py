from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from PyAstronomy import pyasl
'''
August 2014
by Meredith Rawls

Applies a fake binary star RV shift to a star spectrum for funsies.
Or maybe to test FDBinary. Whichever.

Basically, it turns a single star spectrum into a pair of identical stars
that are pretending to orbit each other so I can test spectral decomposition.

There is an option to make a quick plot of the fake SB2 spectra that can be commented out.
Input files and parameters are specified below.
'''

# This is the original non-SB2 (single star) spectrum
star_spectrum = 'arcturus.fits'

# This is a file with timestamps of observations
timefile = 'obstimes.txt'

# Use these options to apply an RV curve from REAL DATA
#gaussfit_file = 'gaussfit_results.txt'
#rvstd = -5.19 # constant RV offset for the star_spectrum, in km/s ... NOT ACTUALLY USED
#bcvs = [4.0466, 6.1887, 10.0119, 13.0188, 12.7790, 12.2969, 11.6241, 9.8934,
#		8.9401, 4.6017, -2.5461, -2.7178, -2.9058, -4.6610, -4.8365, -9.6012,
#		-12.7467, -13.6868, 10.0202, 12.2441, 11.4402, -3.6975, -5.2094]
# ... bcvs are NOT ACTUALLY USED

# Use these options to apply an RV curve from a MODEL
times = []
for line in open(timefile): times.append(float(line))
porb = 171.277
t0 = 321.189576
ecc = 0.36
argper = 18 #degrees! look out!
amp1 = 33.7
amp2 = 33.1
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2

# This is the file that will be created at the end
FDBinary_file = 'infile_fake.obs'

# Boundaries of the new wavelength scale that will be created
logwavemin = 3.591
lnwavemin = np.log(np.power(10,logwavemin))
logwavemax = 3.927
lnwavemax = np.log(np.power(10,logwavemax))
deltalogwave = 0.0000035
deltalnwave = np.log(np.power(10,deltalogwave))
lnwave = np.arange(lnwavemin, lnwavemax, deltalnwave)

# Reads in the star spectrum from a FITS file.
# Copies this so we have two (presently identical) spectra to work with: A and B.
# Removes any constant radial velocity shift from A and B.

hdu = fits.open(star_spectrum)
# spectrum data
spec = hdu[0].data
# header data
head = hdu[0].header
# wavelength scales... 1st is entire spectrum, 2nd is log-wavelength in short region only
wave = np.arange(head['crval1'], head['crval1']+head['cdelt1']*len(spec), head['cdelt1'])


# Reads output file from gaussfit_file to find N pairs of delta-vs.
# (These delta-vs are from some real binary star.)
# OR, if the gaussfit_file hasn't been defined, create a model RV curve instead.

deltavs = []		# all of these are in km/s
deltav1s = []
deltav2s = []
try:
	f1 = open(gaussfit_file)
	f1.close()
	for line in open(gaussfit_file):
		if line[0:6] == 'SHIFT1': deltavs.append(line[7:20])
	for idx, vel in enumerate(deltavs):
		if (idx % 2 == 0): deltav1s.append(float(vel))
		else: deltav2s.append(float(vel))
	if len(deltav1s) != len(deltav2s):
		print('Warning: unequal number of Gaussian peak pairs')
	N = len(deltav1s)

except:
	# create model RV curve instead
	# save the resulting velocity shifts in the deltav1s and deltav2s lists
	N = len(times)
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
#	plt.plot(times, deltav1s, color='0.75', marker='o', mec=red, mfc=red, ls=':', ms=7, mew=0)
#	plt.plot(times, deltav2s, color='0.75', marker='o', mec=yel, mfc=yel, ls=':', ms=7, mew=0)
#	plt.ylabel('Model Radial Velocity (km s$^{-1}$)')
#	plt.xlabel('Time (BJD -- 2454833)')
#	plt.show()


# Applies pairs of delta-vs (shifts) to each set of A and B.

waveAs = []
waveBs = []
specAs = []
specBs = []
c = 2.99792e5 # in km/s, just like deltav1s & deltav2s
for i in range(0, N):
	# unit check: Ang * km/s / km/s + Ang = Ang
	waveAs.append( wave * deltav1s[i] / c + wave )
	waveBs.append( wave * deltav2s[i] / c + wave )
	specAs.append( np.interp(lnwave, np.log(waveAs[i]), spec) )
	specBs.append( np.interp(lnwave, np.log(waveBs[i]), spec) ) 


# Creates the fake SB2 spectra by combining pairs of A and B.
# OPTION: Make a quick plot to make sure nothing terrible happened and SB2 spectra exist.

SB2s = []
ploffset = 0
for i in range(0, N):
	#SB2s.append( (specAs[i] + specBs[i]) / 2 ) # for an equal flux ratio
	SB2s.append( 0.6*specAs[i] + 0.4*specBs[i] ) # for an unequal flux ratio
	plt.plot(np.power(np.exp(1),lnwave), SB2s[i] + ploffset)
	ploffset += 0.5
#plt.show()


# Uses each SB2 to make an infile_fake.obs file, for FDBinary.
#	It will say '# N+1 X len(SB2)' on the 1st line
#	Then, 1st column will be log wavelength from lnwavemin to lnwavemax
#	Subsequent columns will be data from each SB2

fout = open(FDBinary_file, 'w')
print('# ' + str(N+1) + ' X ' + str(len(lnwave)), file=fout)
for j in range(0, len(lnwave)):
	newstring = str(lnwave[j])
	for i in range(0, N):
		newstring += '\t' + str(SB2s[i][j])
	print(newstring, file=fout)
fout.close()