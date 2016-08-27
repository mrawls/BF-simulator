# BF-simulator
### Overview
The goal is to use synthetic spectroscopic binary spectra to see how well we can recover radial velocities, rotation, temperature, etc.

* Given a flux ratio, a pair of PHOENIX spectra, a set of observation times, an a binary orbit (e=0 for now), we will use the broadening function (BF) technique to investigate how well we can recover individual radial velocity measurements, the semiamplitude K, and each star's rotational velocity vsini.
* To accomplish this, we use `create_fakedata.py` to create a simulated time series of double-lined stellar spectra. Comments within this program guide the user to set input spectra, light ratios, vsini, and simulated observation times. The user then chooses either OPTION 1 to read in a list of predetermined radial velocities or OPTION 2 to compute radial velocities based on orbital parameters.
* The output of `create_fakedata.py` is twofold: (1) a series of `fakedata#.txt` two-column text files, one for each simulated SB2 observation (wavelength in Angstroms and normalized flux); and (2) a single `fd3_infile_fakebinary.obs` file formatted for use with the spectral disentangling program FDBinary/fd3 (wavelength in ln(Angstroms) followed by one column per simulated SB2 flux).
* To run, edit the source code as desired, then simply type `python create_fakedata.py`. Dependencies: numpy, matplotlib, astropy, and PyAstronomy.

### Broadening Function example for TYC 3559
For this paper in progress: https://www.overleaf.com/1778823xdcyhk#/4451323/

Real observed spectra of the triple system TYC 3559 (KIC 8848288) are included in this repo. It is a nearly-stationary giant star and a subgiant star with an eclipsing brown dwarf companion. The spectra are a composite of the bright giant and the dimmer, fast-rotating, RV-variable subgiant.

* The program `BF_pythonTYC.py` is a customized copy of https://github.com/mrawls/BF-rvplotter/blob/master/BF_python.py
* The output is `bfoutfile3.txt` (the broadening functions with a two-Gaussian fit) and `rvs_revisited3_BF.txt` (the radial velocities), plus a few plots and some useful info written to stdout. 
* To run, edit `infiles.txt`, `bjdfile.txt`, and `gaussfit.txt` as necessary, then simply type `python BF_pythonTYC.py`. Dependencies: BF_functions (included), numpy, matplotlib, astropy, PyAstronomy, scipy, pandas, gaussfitter (https://github.com/keflavich/gaussfitter).

### Next steps
* Adjust flux ratios (and other parameters?) in `create_fakedata.py`
* Run `BF_pythonTYC.py` for the series of simulated double-lined spectra and inspect the results
* Repeat the above two steps as desired for different sets of simulated observations
* Attempt to run FDBinary/fd3 on one or more sets of simulated spectra and inspect the results (for more on this, see http://sail.zpf.fer.hr/fdbinary/ and https://github.com/mrawls/FDBinary-tools)
