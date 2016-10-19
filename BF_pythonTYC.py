from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.io import fits
from astropy.time import Time
from PyAstronomy import pyasl
from scipy import ndimage
import pandas as pd
import gaussfitter as gf
import BF_functions as bff
'''
Program to extract radial velocities from a double-lined binary star spectrum.
Uses the Broadening Function technique.

Meredith Rawls
2014-2015

This version is especially for TYC 3559!

Based loosely on Rucinski's BFall_IDL.pro, and uses the PyAstronomy tools.
http://www.astro.utoronto.ca/~rucinski/BFdescription.html
http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/svd.html

In practice, you will run this twice: once to do the initial BF, and then again
to properly fit the peaks of each BF with a Gaussian.

INPUT
infiles:    single-column file with one FITS or TXT filename (w/ full path) per line
            1st entry must be for the template star (e.g., arcturus or phoenix model)
            (the same template is used to find RVs for both stars)
            NO comments are allowed in this file
            FUN FACT: unless APOGEE, these should be continuum-normalized to 1 !!!
bjdinfile:     columns 0,1,2 must be filename, BJD, BCV (e.g., from IRAF bcvcorr)
            top row must be for the template star (e.g., arcturus)
            (the 0th column is never used, but typically looks like infiles_BF.txt)
            one line per observation
            comments are allowed in this file using #
gausspars:    your best initial guesses for fitting gaussians to the BF peaks
            the parameters are [amp1, offset1, width1, amp2, offset2, width2]
            the top line is ignored (template), but must have six values
            one line per observation
            comments are allowed in this file using #

OUTPUT
outfile:    a file that will be created with 8 columns: BJD midpoint, orbital phase,
            Kepler BJD, RV1, RV1 error, RV2, RV2 error
bfoutfile:  a file that contains all the BF function data (raw RV, BF, gaussian model)

'''
# (for KIC 8848288, ie TYC 3559)
#infiles =    'infiles.txt'
#bjdinfile =  'bjdfile.txt'
#gausspars =  'gaussfit.txt'
#outfile =    'rvs_revisited3_BF.txt'
#bfoutfile =  'bfoutfile3.txt'

# same as above, but for APOGEE
infiles =    'infiles_apogee.txt'
bjdinfile =  'bjdfile_apogee.txt'
gausspars =  'gaussfit_apogee.txt'
outfile =    'rvs_apogee.txt'
bfoutfile =  'bfoutfile_apogee.txt'

isAPOGEE = True        # toggle to use near-IR stuff, or not
SpecPlot = True         # toggle to plot spectra before BFs, or not
bjdoffset = 2454833.    # difference between real BJDs and 'bjdfunny' (truncated BJDs)
amplimits = [0.8,1, 0,0.2] # limits for gaussian normalized amplitude [min1,max1,min2,max2]
threshold = 10           # margin for gaussian position (raw RV in km/s)
widlimits = [0,15, 0,40]  # limits for gaussian width (km/s) [min1,max1,min2,max2]

period = 5.56648; BJD0 = 2454904.8038 # 8848288 orbital parameters
rvstd = 0; bcvstd = 0 # model template RV is 0

smoothstd = 1.0 #1.5     # stdev of Gaussian to smooth BFs by (~slit width in pixels)
m = 171             # length of the BF (must be longer if RVs are far from 0)
#w00 = 4485; n = 53000; stepV = 1.5
#w00 = 4485; n = 80000; stepV = 1.5 # testing larger, redder wavelength range
w00 = 15145; n = 15000; stepV = 1.5 # APOGEE

rvneg = -74; rvpos = 34; ymin = -0.05; ymax = 1.05 # plot limits
##########

print('Welcome to the Broadening Function party!')
print('')
print('MAKE SURE THIS IS WHAT YOU WANT:')
print('You set Porb = {0} days, BJD0 = {1} days'.format(period, BJD0))

# CREATE NEW SPECTRUM IN LOG SPACE
# This uses w00, n, and stepV, defined above. The new wavelength grid is w1.
# The BF will be evenly spaced in velocity with length m.
# The velocity steps are r (km/s/pix).
w1, m, r = bff.logify_spec(isAPOGEE, w00, n, stepV, m)

# READ IN ALL THE THINGS
specdata = bff.read_specfiles(infiles, bjdinfile, isAPOGEE)
nspec = specdata[0]; filenamelist = specdata[1]
datetimelist = specdata[2]; wavelist = specdata[3]; speclist = specdata[4]

# INTERPOLATE THE TEMPLATE AND OBJECT SPECTRA ONTO THE NEW LOG-WAVELENGTH GRID
# OPTION TO PLOT THIS (commented out for now)
##plt.figure(1)
newspeclist = []
yoffset = 0
if SpecPlot == True:
    plt.axis([w1[0], w1[-1], 0, nspec+3])
    plt.xlabel(r'Wavelength ({\AA})')
for i in range (0, nspec):
    newspec = np.interp(w1, wavelist[i], speclist[i])
    newspeclist.append(newspec)
    if SpecPlot == True:
        plt.plot(w1, newspec+yoffset, label=datetimelist[i].iso[0:10])#, color='b')
    yoffset = yoffset + 1
if SpecPlot == True:
    ##plt.legend()
    plt.show()

# BROADENING FUNCTION TIME
svd = pyasl.SVD()
# Single Value Decomposition
svd.decompose(newspeclist[0], m)
singularvals = svd.getSingularValues()
bflist = []
bfsmoothlist = []
for i in range (0, nspec):
    # Obtain the broadening function
    bf = svd.getBroadeningFunction(newspeclist[i]) # this is a full matrix
    bfarray = svd.getBroadeningFunction(newspeclist[i], asarray=True)
    # Smooth the array-like broadening function
    # 1ST LINE - python 2.7 with old version of pandas; 2ND LINE - python 3.5 with new version of pandas
    #bfsmooth = pd.rolling_window(bfarray, window=5, win_type='gaussian', std=smoothstd, center=True)
    bfsmooth = pd.Series(bfarray).rolling(window=5, win_type='gaussian', center=True).mean(std=smoothstd)
    # The rolling window makes nans at the start because it's a punk.
    for j in range(0,len(bfsmooth)):
        if np.isnan(bfsmooth[j]) == True:
            bfsmooth[j] = 0
        else:
            bfsmooth[j] = bfsmooth[j]
    bflist.append(bf)
    bfsmoothlist.append(bfsmooth)
    
bfnormlist = []
for a in bfsmoothlist:
    bfnormlist.append((a-np.min(a))/(np.max(a)-np.min(a)))

# Obtain the indices in RV space that correspond to the BF
bf_ind = svd.getRVAxis(r, 1) + rvstd - bcvstd

# OPTION TO PLOT THE SINGULAR VALUES TO SEE WHERE THEY AREN'T A MESS
# this probably isn't important, because instead of choosing which values to throw out,
# we use "Route #2" as described by Rucinski and just use the final row of the BF array
# and smooth it with a Gaussian to get rid of noise problems.
# for more info, seriously, read http://www.astro.utoronto.ca/~rucinski/SVDcookbook.html
##plt.figure(2)
#plt.semilogy(singularvals, 'b-')
#plt.xlabel('BF Index')
#plt.ylabel('Singular Values')
#plt.show()

# OPTION TO PLOT THE SMOOTHED BFs
##plt.figure(3)
plt.axis([rvneg, rvpos, -0.2, float(nspec)/2.5])
plt.xlabel('Radial Velocity (km s$^{-1}$)')
plt.ylabel('Broadening Function (arbitrary amplitude)')
yoffset = 0.0
for i in range(1, nspec):
    plt.plot(bf_ind, bfsmoothlist[i]+yoffset, color='b')
    yoffset = yoffset + 0.4
plt.show()

# FIT THE SMOOTHED BF PEAKS WITH TWO GAUSSIANS
# you have to have pretty decent guesses in the gausspars file for this to work.
#bffitlist = bff.gaussparty(gausspars, nspec, filenamelist, bfsmoothlist, bf_ind, threshold)
bffitlist = bff.gaussparty(gausspars, nspec, filenamelist, bfnormlist, bf_ind, amplimits, threshold, widlimits)
rvraw1 = []; rvraw2 = []; rvraw1_err = []; rvraw2_err = []
rvraw1.append(0), rvraw2.append(0), rvraw1_err.append(0), rvraw2_err.append(0)
for i in range(1, len(bffitlist)):
    rvraw1.append(bffitlist[i][0][1]) # [0,1,2] is amp,rv,width for star 1; [4,5,6] is same for star2
    rvraw2.append(bffitlist[i][0][4])
    rvraw1_err.append(bffitlist[i][2][1])
    rvraw2_err.append(bffitlist[i][2][4])

# CALCULATE ORBITAL PHASES AND FINAL RV CURVE
rvdata = bff.rvphasecalc(bjdinfile, bjdoffset, nspec, period, BJD0, rvraw1, rvraw1_err, rvraw2, rvraw2_err, rvstd, bcvstd)
phase = rvdata[0]; bjdfunny = rvdata[1]
rv1 = rvdata[2]; rv2 = rvdata[3]
rv1_err = rvdata[4]; rv2_err = rvdata[5]
g2 = open(outfile, 'w')
print('# RVs calculated with BF_python.py', file=g2)
print('#', file=g2)
print('# Porb = {0} days, BJD0 = {1} days'.format(period, BJD0), file=g2)
print('# Wavelength axis = [{0} - {1}] Angstroms'.format(w1[0], w1[-1]), file=g2)
print('#', file=g2)
print('# Template spectrum (line 0 of infiles):  {0}'.format(filenamelist[0]), file=g2)
print('# RV of template, BCV of template (km/s): {0}, {1}'.format(rvstd, bcvstd), file=g2)
print('#', file=g2)
print('# List of all input spectra (infiles): {0}'.format(infiles), file=g2)
print('# Target BJD and BCV info (bjdinfile): {0}'.format(bjdinfile), file=g2)
print('# Gaussian fit guesses (gausspars):    {0}'.format(gausspars), file=g2)
print('#', file=g2)
print('# BF parameters: w00 = {0}; n = {1}; stepV = {2}'.format(w00, n, stepV), file=g2)
print('# BF parameters: smoothstd = {0}; m = {1}'.format(smoothstd, m), file=g2)
print('# gaussfit: amplimits = {0}; threshold = {1}, widlimits = {2}'.format(amplimits, threshold, widlimits), file=g2)
print('#', file=g2)
print('# time, phase, adjusted_time, RV1 [km/s], error1 [km/s], RV2 [km/s], error2 [km/s]', file=g2)
print('#', file=g2)
for i in range(1, nspec):
    print ('%.9f %.9f %.9f %.5f %.5f %.5f %.5f' % (bjdfunny[i] + bjdoffset, phase[i], bjdfunny[i], 
            rv1[i], rv1_err[i], rv2[i], rv2_err[i]), file=g2)
g2.close()
print('BJD, phase, and RVs written to %s.' % outfile)
print('Use rvplotmaker.py to plot the RV curve.')

try:
    bfout = open(bfoutfile, 'w')
    for idx in range(1, nspec):
        print('###', file=bfout)
        print('# timestamp: {0}'.format(datetimelist[idx]), file=bfout)
        print('# Gaussian 1 [amp, RV +/- err, wid]: [{0:.2f}, {1:.2f} +/- {2:.2f}, {3:.2f}]'.format(bffitlist[i][0][0], rvraw1[i], rvraw1_err[i], bffitlist[i][0][2]), file=bfout)
        print('# Gaussian 2 [amp, RV +/- err, wid]: [{0:.2f}, {1:.2f} +/- {2:.2f}, {3:.2f}]'.format(bffitlist[i][0][3], rvraw2[i], rvraw2_err[i], bffitlist[i][0][5]), file=bfout)
        print('# Uncorrected_RV, BF_amp, Gaussian_fit', file=bfout)
        print('###', file=bfout)
        for vel, amp, modamp in zip(bf_ind, bfsmoothlist[idx], bffitlist[idx][1]):
            print(vel, amp, modamp, file=bfout)
    bfout.close()
except:
    print('No BF outfile specified, not saving BF data to file')
    
# handy little gaussian function maker
def gaussian(x, amp, mu, sig): # i.e., (xarray, amp, rv, width)
    return amp * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# PLOT THE FINAL SMOOTHED BFS + GAUSSIAN FITS IN INDIVIDUAL PANELS
# manually adjust this multi-panel plot based on how many spectra you have
#plt.figure(4)
windowcols = 3 #4                                # how many window columns there should be
#windowrows = 6
windowrows = int([np.rint((nspec-1)/windowcols) if (np.float(nspec-1)/windowcols)%windowcols == 0 else np.rint((nspec-1)/windowcols)+1][0])
xmin = rvneg
xmax = rvpos
#gaussxs = np.arange(-200, 200, 0.1)
fig = plt.figure(1, figsize=(15,10))
fig.text(0.5, 0.04, 'Uncorrected Radial Velocity (km s$^{-1}$)', ha='center', va='center', size=26)
fig.text(0.07, 0.5, 'Broadening Function', ha='center', va='center', size=26, rotation='vertical')
for i in range (1,nspec):
    ax = fig.add_subplot(windowrows, windowcols,i) # out of range if windowcols x windowrows < nspec
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    if i!=1 and i!=5 and i!=9 and i!=13 and i!=17 and i!=21 and i!=25:
        ax.set_yticklabels(())
    #if i!=20 and i!=21 and i!=22 and i!=23 and i!=24 and i!=25:
    if i < nspec-windowrows:
    #if i!=13 and i!=14 and i!=15 and i!=16:
        ax.set_xticklabels(())
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.axis([xmin, xmax, ymin, ymax])
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.text(xmax - 0.25*(np.abs(xmax-xmin)), 0.8*ymax, '%.3f $\phi$' % (phase[i]), size=12)
    plt.text(xmax - 0.35*(np.abs(xmax-xmin)), 0.6*ymax, '%s' % (datetimelist[i].iso[0:10]), size=12)
    #plt.plot(bf_ind, bfsmoothlist[i], color='k', lw=1.5, ls='-', label='Smoothed BF')
    plt.plot(bf_ind, bfnormlist[i], color='k', lw=1.5, ls='-', label='Normalized Smoothed BF')
    plt.plot(bf_ind, bffitlist[i][1], color='b', lw=2, ls='--', label='Two Gaussian fit')
    gauss1 = gaussian(bf_ind, bffitlist[i][0][0], bffitlist[i][0][1], bffitlist[i][0][2])
    gauss2 = gaussian(bf_ind, bffitlist[i][0][3], bffitlist[i][0][4], bffitlist[i][0][5])
    plt.plot(bf_ind, gauss1, color='#e34a33', lw=2, ls='--')#, label='Gaussian fit 1')
    plt.plot(bf_ind, gauss2, color='#fdbb84', lw=2, ls='--')#, label='Gaussian fit 2')
    # OPTION TO PLOT VERTICAL LINE AT ZERO
    plt.axvline(x=0, color='0.75')    
    # print legend
    if i==nspec-1: ax.legend(bbox_to_anchor=(2.6,0.7), loc=1, borderaxespad=0., 
                        frameon=False, handlelength=3, prop={'size':20})
plt.show()