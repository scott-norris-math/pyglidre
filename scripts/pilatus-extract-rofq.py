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
from lmfit import Parameters			# not a standard library -- install using PIP
import pygisaxs.wrappers.pilatus as pil		# from my (very small) library
import pygisaxs.utilities as util


# geometry of interest in Pilatus Files
#pilatus_rectangle = [361,7]
#valid_pixels = range(10,166)+range(195,351)


# fitting function and parameter guesses
fnc2min = util.extract_rofq.ode1sol
localguess = Parameters()
localguess.add('R',    value=-5e-4)
localguess.add('I0',   value= 1e-4)
localguess.add('beta', value= 1e-8)

# construct global guess for hierarchical fitting
globalguess = None
globalguess = Parameters()		# COMMENT OUT THESE TWO LINES TO FIT ALL PARAMETERS LOCALLY
globalguess.add('beta', value=1e-8)	# (useful to identify potential global parameters and guesses)


# in principle this can be run over multiple directories in sequence (just to save typing and waiting)
directories = sys.argv[1:]
for directory in directories:


  # extract the intensity from all of the Pilatus output files
  tqi = pil.extract_iofqt_2D(directory, usecache=False)

  #frame_prefix = "%s/intensity-2D" % (directory)
  #util.plotting.create_2D_intensity_frames(tqi, prefix=frame_prefix)

  # convert the data from 2D to 1D
  tqi = util.collapse.collapse_2D_to_1D(tqi, vspan=7, qmin=-1.5, qmax=1.5)

  # made some basic plots of the raw intensity I(q,t)
  util.plotting.plot_iofqt_heatmap(tqi, filename="%s/iofqt-heatmap-1-raw.png" % (directory))
  util.plotting.plot_intensity_progression(tqi, filename="%s/iofqt-manyplots-1-raw.png" % (directory))

  # optional shifting step
  #tqi = util.shift.shift_by_center_mass_2D(tqi)

  # scale by peak volume, re-plot progressions
  tqi = util.scale.scale_by_all_1D(tqi)
  util.plotting.plot_iofqt_heatmap(tqi, filename="%s/iofqt-heatmap-2-scaled.png" % (directory))
  util.plotting.plot_intensity_progression(tqi, filename="%s/iofqt-manyplots-2-scaled.png" % (directory))

  # subtract high-q noise, re-plot progressions
  tqi = util.scale.subtract_vampires_1D(tqi, directory=directory)
  util.plotting.plot_iofqt_heatmap(tqi, filename="%s/iofqt-heatmap-3-denoised.png" % (directory))
  util.plotting.plot_intensity_progression(tqi, filename="%s/iofqt-manyplots-3-denoised.png" % (directory))

  # set up valid times, pixels, and pixels for global fitting
  valid_times  = range(1,len(tqi["tvals"]))
  valid_pixels = np.where( (np.abs(tqi["qvals"]) > 0.1) * (np.abs(tqi["qvals"]) < 1.25) )[0]
  beta_fit_pixels = np.where( (np.abs(tqi["qvals"]) > 0.15) * (np.abs(tqi["qvals"]) < 0.45) )[0]

  # Now fit the data to our ODE solution
  print "\n\nextracting R(q) ... " ; sys.stdout.flush()
  fits  = util.extract_rofq.fit_rofq_hierarchically(tqi, fnc2min, localguess, globalguess, t_indices=valid_times, q_indices=valid_pixels, image_directory=directory)
  print "R(q) extraction done.";  sys.stdout.flush()


  # pull data out of the dictionary for ease of reference
  image_directory = "./"
  qvals  = np.array( fits["qvals"].values )
  Rvals  = np.array( fits["R"].values )
  Rstds  = np.array( fits["R"].stderrs )
  I0vals = np.array( fits["I0"].values )
  I0stds = np.array( fits["I0"].stderrs )
  chisqr = np.array( fits["chisqr"].values )
  Bvals  = np.array( fits["beta"].values )
  Bstds  = np.array( fits["beta"].stderrs )


  # summary figure
  if image_directory != None:
    plt.figure(figsize=(32,16))

    plt.subplot(221)
    plt.errorbar(qvals, Rvals, yerr=Rstds, fmt='bs')
    plt.axis((qvals[0], qvals[-1], -0.004, 0.0005))
    plt.title('Fit for R(q) -- %s' % (directory))

    plt.subplot(224)
    plt.plot(qvals, np.log10( chisqr ), 'bs')
    #plt.axis((qvals[0], qvals[-1], -0.002, 0.0005))
    plt.title('chi squared error(q)')

    plt.subplot(222)
    logI0vals = np.log10(I0vals)
    logI0stds = np.log10((I0vals+I0stds)/(I0vals-I0stds)) / 2.0
    plt.errorbar(qvals, logI0vals, yerr=logI0stds, fmt='bs')
    Imin = np.percentile(np.log10(I0vals), 05)
    Imax = np.percentile(np.log10(I0vals), 95)
    #plt.axis((qvals[0], qvals[-1], Imin, Imax))
    plt.title('log10 of initial intensity I0(q)')

    plt.subplot(223)
    logBvals = np.log10(Bvals)
    logBstds = np.log10((Bvals+Bstds)/(Bvals-Bstds)) / 2.0
    plt.errorbar(qvals, logBvals, yerr=logBstds, fmt='bs')
    Bmin = np.percentile(np.log10(Bvals), 05)
    Bmax = np.percentile(np.log10(Bvals), 95)
    #plt.axis((qvals[0], qvals[-1], Bmin, Bmax))
    plt.title('log10 of background noise B(q)')

    plt.savefig("%s/rofq_summary.png" % (directory))
    plt.close()

  print "\n\n"

  continue

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



