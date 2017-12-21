import os
import sys
import shelve
import numpy as np
from PIL import Image
import pygisaxs.utilities.smooth as psmooth


def extract_iofqt_2D(directory, generate_pngs=False, imagesize=(487,195), usecache=False):


  print "extracting intensity ... ", 
  sys.stdout.flush()

  if (usecache):
    cachefile = "%s/raw-data.pkl" % (directory)

    if os.path.isfile(cachefile):
      f = shelve.open(cachefile)
      tqi = f["tqi"]
      f.close()

      print "Found cache file %s.  I AM NOT RE-READING THE PILATUS FILES!" % (cachefile)
      return tqi




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

  # pilatus image size
  xn = imagesize[0]
  zn = imagesize[1]

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
  fluencelist = []
  for utk, utime in enumerate(uniquetimes):

    #print "For time", utime
    images  = 0
    fluence = 0

    array2d = np.zeros((487,195))
    while (times[configpos] == utime):
      pfile    = "%s/%s" % (directory, files[configpos])
      img      = Image.open(pfile)
      array2d += np.transpose(np.asarray(img)) 
      fluence += xrays[configpos]
      images  += 1

      # update
      configpos += 1
      if (configpos >= len(times)):
        break

    print "Found unique time %10d. Averaged %3d image files." % (utime, images)
    #cleaned_array2d = psmooth.remove_hot_pixels_2d(array2d)
    array2dlist.append(array2d)
    fluencelist.append(fluence)

  # convert to numpy array
  array2d = np.array(array2dlist)
  fluences = np.array(fluencelist)


  # get total intensity over all times
  totalI = np.sum(array2d, axis=0)


  # if a peak box file was created, then read it
  peakbox = None
  peakwindowfile = "%s/specular-peak-window.txt" % (directory)
  if (os.path.isfile(peakwindowfile)):
    peakbox = np.loadtxt(peakwindowfile)

  # find the location of the specular peak
  if peakbox== None:

    # look over the whole image for the specular peak
    max_coords = np.unravel_index(np.argmax(totalI), totalI.shape)
    x0 = max_coords[0]
    z0 = max_coords[1]
    #print "Image", pilatus_file, ".  No specular peak window provided -- considering whole image.  Peak found at ", max_coords, " with value", arr[x0,z0]

  else:
    # look only in the provided center_window for the specular peak
    cxmin = peakbox[0]
    cxmax = peakbox[1]
    cymin = peakbox[2]
    cymax = peakbox[3]
    cwindow = totalI[cxmin:cxmax,cymin:cymax]
    temp_coords = np.unravel_index(np.argmax(cwindow), cwindow.shape)
    max_coords = np.array([cxmin, cymin]) + temp_coords

    x0 = max_coords[0]
    z0 = max_coords[1]
    print "Looking for specular peak in provided window ", peakbox, ".  Peak found at ", max_coords #, " with value", arr[x0,z0]
    if x0 == cxmin or x0 == cxmax or z0 == cymin or z0==cymax:
      print "WARNING: specular peak found on edge of provided box. \n Cannot confirm that value is a local maximum. \n Identified co-ordinates may be inaccurate."


  # assign q-values 
  xpvals  = np.arange(xn) - x0
  zpvals  = np.arange(zn) - z0
  qxvals   = (2*3.14159 * np.sin(3.14159*xpvals/(131*180)) / .124)  # in 1/nm 
  qzvals   = (2*3.14159 * np.sin(3.14159*zpvals/(131*180)) / .124)  # in 1/nm 
 

  # put the results into a dictionary and return it
  results = dict()
  results["type"]     = "PILATUS"
  results["tvals"]    = np.array(uniquetimes)
  results["qxvals"]   = qxvals
  results["qzvals"]   = qzvals
  results["qvals"]    = qxvals
  results["fluences"] = fluences
  results["iofqt"]    = array2d

  print "intensity extraction done."
  sys.stdout.flush()

  if (usecache):
    cachefile = "%s/raw-data.pkl" % (directory)
    f = shelve.open(cachefile)
    f["tqi"] = results
    f.close()
    print "Stored results in cache file %s." % (cachefile)

  return results

