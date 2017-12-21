'''
This extracts estimates of R(q) from sequences of TIFF images produced by PILATUS detectors.  At the moment, you should make sure that you run 

* pilatus-process-scans.py
* pilatus-create-simple-timefile.py

BEFORE running this script.
'''

# basic imports
import os
import sys
import shelve
import numpy as np
import matplotlib as mpl ; mpl.use('Agg') ;
import matplotlib.pyplot as plt

# custom imports
from lmfit import Parameters		# not a standard library -- install using PIP
import pygisaxs.wrappers.madi as madi		# from my (very small) library
import pygisaxs.utilities as util	# from my (very small) library


# fitting function and parameter guesses
fnc2min = util.extract_rofq.ode1sol_plus_constant
localguess = Parameters()
localguess.add('R',    value=-5e-4)
localguess.add('I0',   value= 1e-4, min=0.0)
localguess.add('beta', value= 5e-8, min=0.0)
localguess.add('background', value=0.0, min=0.0)


# extract directory-dependent fitting information from configuration file
# ----------------------------------------------------------------------------
configfile  = "extraction-parameters.csv" 
mydtype=[
	('directory', '|S64'),
	('orientation', '|S8'),
	('angle', '|S8'),
	('start_time', '<i4'),
	('qmin', '<f8'),
	('qmax', '<f8'),
	('beta_guess', '<f8')
]

try:
  config_table = np.loadtxt(
	configfile,
	delimiter=",", 
	usecols=[0,1,2,3,4,5,6], 
        skiprows=1,
	dtype=mydtype
  )
except Exception as EE:
  print EE
  print "There seems to be some kind of problem."
  print "Either you did not create a configuration file,"
  print "or the file is not in the expected format.  Exiting."
  exit()

config_data = dict()
for row in config_table:
  key = row[0]
  entry = dict()

  column = 1
  for dt in mydtype[1:]:
    inner_key = dt[0]
    inner_val = row[column]
    entry[inner_key] = inner_val
    column += 1
    
  config_data[key] = entry
  print key, entry


# in principle this can be run over multiple directories in sequence (just to save typing and waiting)
directories = sys.argv[1:]
for directory in directories:


  # intensity extraction section
  # -----------------------------
  print "-----------------------------\n extracting intensity for directory %s ... \n ----------------------------\n" % (directory); sys.stdout.flush()

  # extract raw data, plot progressions
  tqi = madi.extract_iofqt_simple(directory)
  util.plotting.plot_iofqt_heatmap(tqi, filename="%s/iofqt-heatmap-1-raw.png" % (directory))
  util.plotting.plot_intensity_progression(tqi, filename="%s/iofqt-manyplots-1-raw.png" % (directory))

  # optional shifting step
  #tqi = util.shift.shift_by_center_mass_2D(tqi)

  # scale by peak volume, re-plot progressions
  tqi = util.scale.scale_by_peak_1D(tqi)
  util.plotting.plot_iofqt_heatmap(tqi, filename="%s/iofqt-heatmap-2-scaled.png" % (directory))
  util.plotting.plot_intensity_progression(tqi, filename="%s/iofqt-manyplots-2-scaled.png" % (directory))


  # subtract high-q noise, re-plot progressions
  tqi = util.scale.subtract_vampires_1D(tqi, directory=directory)
  util.plotting.plot_iofqt_heatmap(tqi, filename="%s/iofqt-heatmap-3-denoised.png" % (directory))
  util.plotting.plot_intensity_progression(tqi, filename="%s/iofqt-manyplots-3-denoised.png" % (directory))


  # unpack for later use
  tvals = tqi["tvals"]
  qvals = tqi["qvals"]
  iofqt = tqi["iofqt"]


  # R(q) fitting section
  # -----------------------------

  # get valid times from global configuration data
  tn = len(tvals)
  start_ti = config_data[directory]["start_time"]
  valid_times = range(start_ti, tn)

  # get valid pixels from global configuration data
  qmin = config_data[directory]["qmin"]
  qmax = config_data[directory]["qmax"]
  valid_pixels = np.where( (np.abs(qvals) > qmin) * (np.abs(qvals) < qmax) )[0]

  # get guesses for beta and background 
  betaguess = config_data[directory]["beta_guess"]
  BGguess = np.min(iofqt[:,valid_pixels]) / 4.0

  # global parameters for hierarchical fitting
  globalguess = None
  globalguess = Parameters()				# COMMENT OUT THESE TWO LINES TO FIT ALL PARAMETERS LOCALLY
  globalguess.add('beta', value=betaguess, min=0.0)	# (useful to identify potential global parameters and guesses)
  globalguess.add('background', value=BGguess, min=0.0) # 

  # Now fit the data to our ODE solution
  print "extracting R(q) ... " ; sys.stdout.flush()
  fits  = util.extract_rofq.fit_rofq_hierarchically(tqi, fnc2min, localguess, globalguess, t_indices=valid_times, q_indices=valid_pixels, image_directory=directory)



  # reporting section
  # ----------------------

  # pull data out of the dictionary for ease of reference
  image_directory = "./"
  qfits  = np.array( fits["qvals"].values )
  Rvals  = np.array( fits["R"].values )
  Rstds  = np.array( fits["R"].stderrs )
  Bvals  = np.array( fits["beta"].values )
  Bstds  = np.array( fits["beta"].stderrs )
  Ivals = np.array( fits["I0"].values )
  Istds = np.array( fits["I0"].stderrs )
  chisqr = np.array( fits["chisqr"].values )


  # summary figure
  # -----------------------
  if image_directory != None:
    plt.figure(figsize=(24,13.5))

    plt.subplot(221)
    plt.errorbar(qfits, Rvals, yerr=Rstds, fmt='bs')
    plt.axis((qfits[0], qfits[-1], -0.004, 0.0005))
    plt.title('Fit for R(q) -- %s' % (directory))

    plt.subplot(222)
    logIvals = np.log10(Ivals)
    logIstds = np.log10((Ivals+Istds)/(Ivals-Istds)) / 2.0
    plt.errorbar(qfits, logIvals, yerr=logIstds, fmt='bs')
    Imin = np.percentile(np.log10(Ivals), 05)
    Imax = np.percentile(np.log10(Ivals), 95)
    plt.xlim([qfits[0], qfits[-1]])
    plt.ylim([-8, 0])
    #plt.axis((qfits[0], qfits[-1], Imin, Imax))
    plt.title('log10 of initial intensity I0(q)')

    plt.subplot(223)
    logBvals = np.log10(Bvals)
    logBstds = np.log10((Bvals+Bstds)/(Bvals-Bstds)) / 2.0
    plt.errorbar(qfits, logBvals, yerr=logBstds, fmt='bs')
    Bmin = np.percentile(np.log10(Bvals), 05)
    Bmax = np.percentile(np.log10(Bvals), 95)
    plt.xlim([qfits[0], qfits[-1]])
    plt.ylim([-12, -4])
    #plt.axis((qfits[0], qfits[-1], Bmin, Bmax))
    plt.title('log10 of background noise B(q)')

    plt.subplot(224)
    plt.plot(qfits, np.log10( chisqr ), 'bs')
    plt.xlim([qfits[0], qfits[-1]])
    plt.ylim([-12, -4])
    #plt.axis((qfits[0], qfits[-1], -0.002, 0.0005))
    plt.title('chi squared error(q)')


    plt.savefig("%s/rofq_summary.png" % (directory))
    plt.savefig("%s/rofq_summary.svg" % (directory))
    plt.close()

  print "\n\n"


  # save results
  # --------------------
  results = shelve.open('results-%s.shelf'%(directory))
  results['tvals'] = tvals
  results['qvals'] = qvals
  results['iofqt'] = iofqt

  results['qfits'] = qfits
  results['Rvals'] = Rvals
  results['Rstds'] = Rstds
  results['Bvals'] = Bvals
  results['Bstds'] = Bstds
  results['Ivals'] = Ivals
  results['Istds'] = Istds
  results['chisqr'] = chisqr
  results.close()





