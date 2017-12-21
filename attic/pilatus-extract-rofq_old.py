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
import matplotlib.pyplot as plt

# custom imports
from lmfit import Parameters		# not a standard library -- install using PIP
import pygisaxs.io.pilatus as pil	# from my (very small) library
import pygisaxs.utilities as util	# from my (very small) library


# geometry of interest in Pilatus Files
pilatus_rectangle = [361,7]
valid_pixels = range(0,166)+range(195,361)

# fitting function and parameter guesses
fnc2min = util.ode1sol
localguess = Parameters()
localguess.add('R',    value=-5e-4)
localguess.add('beta', value= 1e-8, min=0.0)
localguess.add('I0',   value= 1e-4, min=0.0)

# global parameters for hierarchical fitting
globalguess = None
globalguess = Parameters()			# COMMENT OUT THESE TWO LINES TO FIT ALL PARAMETERS LOCALLY
globalguess.add('beta', value=1e-8, min=0.0)	# (useful to identify potential global parameters and guesses)



# in principle this can be run over multiple directories in sequence (just to save typing and waiting)
directories = sys.argv[1:]
for directory in directories:


  # extract the intensity from all of the Pilatus output files
  print "extracting intensity ... " 
  sys.stdout.flush()
  tqi = pil.extract_iofqt(directory, window_shape=pilatus_rectangle)
  print "intensity extraction done."
  sys.stdout.flush()




  # made some basic plots of the raw intensity I(q,t)
  util.plotting.plot_iofqt_heatmap(tqi, filename="%s/iofqt-heatmap-2-corrected.png" % (directory))
  util.plotting.plot_intensity_progression(tqi, filename="%s/iofqt-manyplots-2-corrected.png" % (directory))




  # optional shifting step



  # optional scaling step
  tqi = util.scale.scale_by_peak_2D(tqi, directory_name=directory)


  # made some basic plots of the scaled-shifted intensity I(q,t)
  util.plotting.plot_iofqt_heatmap(tqi, filename="%s/iofqt-heatmap-2-corrected.png" % (directory))
  util.plotting.plot_intensity_progression(tqi, filename="%s/iofqt-manyplots-2-corrected.png" % (directory))
  

  # Now fit the data to our ODE solution
  print "extracting R(q) ... " ; sys.stdout.flush()
  fits  = util.extract_rofq.fit_rofq_hierarchically(tqi, fnc2min, localguess, globalguess, q_indices=valid_pixels, image_directory=None)
  print "R(q) extraction done.";  sys.stdout.flush()



  # pull data out of the dictionary for ease of reference
  image_directory = "./"
  qvals  = np.array( fits["qvals"].values )
  Rvals  = np.array( fits["R"].values )
  Rstds  = np.array( fits["R"].stderrs )
  Bvals  = np.array( fits["beta"].values )
  Bstds  = np.array( fits["beta"].stderrs )
  I0vals = np.array( fits["I0"].values )
  I0stds = np.array( fits["I0"].stderrs )
  chisqr = np.array( fits["chisqr"].values )



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




