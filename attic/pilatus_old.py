import os
import glob
import Image
import numpy as np
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt

import pygisaxs.utilities as util


def extract_iofqt_simple(directory, vspan, generate_pngs=False):


  # get configuration file and extract raw data
  # --------------------------------------------
  configfile  = "%s/pilatus-config.csv" % (directory)
  try:
    config_data = np.loadtxt(
	configfile,
	delimiter=",", 
	usecols=[0,1,2,3,4], 
	dtype=[
		('time', '<f8'),
		('scan', '<i4'),
		('image','<i4'),
		('filename','|S64'),
		('flux','<i8')
	]
    )
  except Exception as EE:
    print EE
    print "There seems to be some kind of problem."
    print "Either you did not create a configuration file,"
    print "or the file is not in the expected format.  Exiting."
    exit()

  # extract raw data
  times = [item[0] for item in config_data]
  files = [item[3] for item in config_data]
  xrays = [item[4] for item in config_data]

  # find unique times
  uniquetimes = np.unique(times)
  utn = len(uniquetimes)

  # get averaged 2d arrays
  # -----------------------
  configpos = 0
  array2dlist = []
  for utk, utime in enumerate(uniquetimes):

    #print "For time", utime
    images  = 0
    fluence = 0
    photons = 0

    array2d = np.zeros((487,195))
    while (times[configpos] == utime):
      pfile    = "%s/%s" % (directory, files[configpos])
      img      = Image.open(pfile)
      array2d += np.transpose(np.asarray(img)) 

      fluence += xrays[configpos]
      photons += np.sum(array2d)
      images  += 1

      # update
      configpos += 1
      if (configpos >= len(times)):
        break

    print "Found unique time %10d. Averaged %3d image files." % (utime, images)
    #array2d /= float(fluence)
    array2d /= float(photons)
    array2dlist.append(array2d)


  # convert to numpy array
  array2d = np.array(array2dlist)

  # if a peak box file was created, then read it
  peakbox = None
  peakwindowfile = "%s/specular-peak-window.txt" % (directory)
  if (os.path.isfile(peakwindowfile)):
    peakbox = np.loadtxt(peakwindowfile)

  # compress in the qz direction
  iofqt = util.convert_2D_to_1D(array2d, average_vspan=vspan, search_window=peakbox)

  # Identify the qx location of the specular peak (rough approx)
  totalI = np.sum(iofqt, axis=0)
  peak   = np.argmax(totalI)

  # assign q-values 
  pvals  = np.arange(xn) - peak
  qvals  = (2*3.14159 * np.sin(3.14159*pvals/(131*180)) / .124)  # in 1/nm 

  # put the results into a dictionary and return it
  results = dict()
  results["tvals"] = np.array(uniquetimes)
  results["qvals"] = np.array(qvals)
  results["iofqt"] = iofqt
  return results







def extract_iofqt(directory, window_shape, generate_pngs=False):
  #
  # Currently, this function expects to be able to find a file called "pilatus-config.csv" in the directory
  # containing all of the TIFF files.  This is a csv file with one line per image that you care about,
  # whose first five columns must contain the following information associated with that image:
  #
  #   col. 1:  experiment time in seconds (*not* PILATUS internal times)
  #   col. 2:  number of the PILATUS scan
  #   col. 3:  number of the PILATUS sub-scan
  #   col. 4:  TIFF filename
  #   col. 5:  fluence associated with that image
  #
  # This function will process the indicated files, averaging together any images that occur at the same
  # experiment time, and return a dictionary $results containing three elements:
  #
  #   results["tvals"]: times at which data are taken (1D array of unique experiment times -- *not* PILATUS times)
  #   results["qvals"]: q-values at which data are taken (1D array)
  #   results["iofqt"]: the 2D array Intensity(q, t)
  #
  # You need to supply the directory in which the files are located, and also the shape of the window you are
  # going to average over to get the intensity.
  #
  # Finally, there is a diagnostic option generate_pngs, which defaults to False.  
  # This will generate a *lot* of images and should only be used if things aren't going as you think they should.
  #




  # get configuration file and extract raw data
  configfile  = "%s/pilatus-config.csv" % (directory)
  try:
    config_data = np.loadtxt(
			configfile, 
			delimiter=",", 
			usecols=[0,1,2,3,4], 
			dtype=[('time', '<f8'),('scan', '<i4'),('image','<i4'),('filename','|S64'),('flux','<i8')]
			)
  except Exception as EE:
    print EE
    print "There seems to be some kind of problem."
    print "Either you did not create a configuration file,"
    print "or the file is not in the expected format.  Exiting."
    exit()

  # extract raw data
  times = [item[0] for item in config_data]
  files = [item[3] for item in config_data]
  xrays = [item[4] for item in config_data]

  # if a peak box file was created, then read it
  peakbox = None
  peakwindowfile = "%s/specular-peak-window.txt" % (directory)
  if (os.path.isfile(peakwindowfile)):
    peakbox = np.loadtxt(peakwindowfile)

  # Plot of the overall x-ray fluence with time (for diagnostic use)
  plt.figure(1)
  plt.plot(times, xrays)
  plt.title('x-ray fluence over run duration')
  plt.savefig("%s/xray-fluence.png" % (directory))
  plt.close()

  # Allocate storage for 2D intensity array, allowing for repeated measurements
  uniquetimes = np.unique(times)
  utn = len(uniquetimes)
  iofqt = np.empty((utn, window_shape[0]))

  # MAIN LOOP:  loop over contents of time, filename, and xray vectors
  configpos = 0
  for utk, utime in enumerate(uniquetimes):

    #print "For time", utime
    images  = 0
    fluence = 0
    photons = 0

    array2d = np.zeros((487,195))
    while (times[configpos] == utime):
      pfile    = "%s/%s" % (directory, files[configpos])
      img      = Image.open(pfile)
      array2d += np.transpose(np.asarray(img)) 

      fluence += xrays[configpos]
      photons += np.sum(array2d)
      images  += 1

      # update
      configpos += 1
      if (configpos >= len(times)):
        break

    print "Found unique time %10d. Averaged %3d image files." % (utime, images)
    #array2d /= float(fluence)
    array2d /= float(photons)
    iofqt[utk,:] = extract_single_iofq(array2d, window_shape, peakbox, generate_pngs)


  # store and return results
  wwidth = window_shape[0]
  wrange = np.arange(wwidth)
  wmidpt = (wwidth-1)/2
  pixels = wrange - wmidpt

  results = dict()
  results["tvals"] = np.array(uniquetimes)
  results["qvals"] = (2*3.14159*np.sin(3.14159*pixels/(131*180))/.124)  # in 1/nm -- FROM EITAN OR CHARBEL MANY MONTHS AGO.
  results["iofqt"] = iofqt
  return results










# REAL WORK OF IMAGE EXTRACTION IS LOCATED IN THIS FILE

def extract_single_iofq(arr, window_shape, center_window=None, generate_pngs=False):
    #
    # This function reads a single PILATUS tiff file, and extracts an intensity reading by averaging a specified 
    # number of pixels above and below the specular peak.  
    #
    # In general the function first applies hot-pixel removal to deal with slightly faulty detectors, and then
    # finds the brightest pixel in the image and calls it the specular peak.  If you know the specular peak will
    # be the brightest pixel, then you can just let the function find it for you (and it does so on an 
    # image-by-image basis, so it automatically handles slight jitter between readings).
    #
    # However, if the specular peak is **not** the brightest pixel in the image, you can specify a smaller 
    # sub-window in which to look for it, by listing four numbers in a file called "specular-peak-window.txt"
    #


    # obtain a median-filtered version of array (eliminates broken pixels -- borrowed from StackExchange)
    blurred = filters.median_filter(arr, size=2)
    difference = arr - blurred
    threshold = 3*np.std(difference)  					# too large and you don't remove pixels -- too small and you smooth
    hot_pixels = np.nonzero((np.abs(difference[1:-1,1:-1])>threshold) ) # ignore edges for simplicity; we'll discard anyway
    hot_pixels = np.array(hot_pixels) + 1 				# because we ignored the first row and first column
    farr       = np.copy(arr) 						# This is the image with the hot pixels removed
    for x,y in zip(hot_pixels[0], hot_pixels[1]):
      farr[x,y]=blurred[x,y]

    # find the location of the specular peak
    if center_window == None:

      # look over the whole image for the specular peak
      max_coords = np.unravel_index(np.argmax(arr), arr.shape)
      x0 = max_coords[0]
      z0 = max_coords[1]
      #print "Image", pilatus_file, ".  No specular peak window provided -- considering whole image.  Peak found at ", max_coords, " with value", arr[x0,z0]

    else:
      # look only in the provided center_window for the specular peak
      cxmin = center_window[0]
      cxmax = center_window[1]
      cymin = center_window[2]
      cymax = center_window[3]
      cwindow = arr[cxmin:cxmax,cymin:cymax]
      temp_coords = np.unravel_index(np.argmax(cwindow), cwindow.shape)
      max_coords = np.array([cxmin, cymin]) + temp_coords

      x0 = max_coords[0]
      z0 = max_coords[1]
      print "Looking for specular peak in provided window ", center_window, ".  Peak found at ", max_coords, " with value", arr[x0,z0]

      if x0 == cxmin or x0 == cxmax or z0 == cymin or z0==cymax:
        print "WARNING: specular peak found on edge of provided box. \n Cannot confirm that value is a local maximum. \n Identified co-ordinates may be inaccurate."


    # identify the vertical analysis window (specular peak +- half the width/height specified in window_shape)
    zmin = z0-(window_shape[1]-1)/2
    zmax = z0+(window_shape[1]-1)/2
    zwindow = np.arange(zmin, zmax+1)

    # get a z-average of the data
    farr_v_window = farr[:, zmin:zmax+1]
    ivals = np.sum(farr_v_window, axis=1)

    # now identify the center of the average
    x0 = np.argmax(ivals)
    xmin = x0-(window_shape[0]-1)/2  
    xmax = x0+(window_shape[0]-1)/2
    windowed_ivals = ivals[xmin:xmax+1]





    # DIAGNOSTIC IMAGES -- ENABLE IF THINGS AREN'T WORKING
    if (generate_pngs):

      fileparts = str(pilatus_file).split('/')
      just_file = fileparts[-1]
      directory = pilatus_file[0:-len(just_file)]

      pngname = str(just_file).replace('tiff', 'png')

      # RAW IMAGE
      plt.figure(1)
      plt.clf()
      plt.imshow(np.transpose(np.log10(arr)), origin="lower")
      plt.title("log10 of Raw 2D GISAXS output")
      plt.savefig('%s/raw-%s' % (directory, pngname))

      # AFTER HOT-PIXEL REMOVAL
      plt.figure(1)
      plt.clf()
      plt.imshow(np.transpose(np.log10(farr)), origin="lower")
      plt.title("log10 of Smoothed 2D GISAXS output")
      plt.savefig('%s/fixed-%s' % (directory, pngname))
      plt.close()

      # IDENTIFIED WINDOW
      plt.figure(1)
      plt.clf()
      plt.imshow(np.transpose(np.log10(farr_v_window)), origin="lower")
      plt.title("log10 of Smoothed Windowed 2D GISAXS output")
      plt.savefig('%s/fixed-window-%s' % (directory, pngname))
      plt.close()

      # AVERAGE OVER WINDOW
      plt.figure(1)
      plt.clf()
      plt.plot(qvals, np.log10(ivals))
      plt.title("log10 of 1D averaged GISAXS output")
      plt.savefig('%s/intensity-%s' % (directory, pngname))
      plt.close()





    return windowed_ivals



