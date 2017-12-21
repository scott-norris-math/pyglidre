import numpy as np


def extract_iofqt_simple(directory):

  # identify the files of the directory format
  timefile  = "%s/initial_time.txt" % (directory)
  pixelfile = "%s/xcoordinates.txt" % (directory)
  datafile  = "%s/initial_intensity.txt" % (directory)

  # read the files in the directory
  tvals = np.loadtxt(timefile)
  xvals = np.loadtxt(pixelfile)
  iofqt = np.transpose(np.loadtxt(datafile))

  # Identify the location of the specular peak (rough approx)
  totalI = np.sum(iofqt, axis=0)
  peak   = np.argmax(totalI)

  # assign q-values 
  pvals  = np.arange(len(xvals)) - peak
  qvals  = (2*3.14159 * np.sin(3.14159*pvals/(131*180)) / .124)  # in 1/nm

  # put the results into a dictionary and return it
  results = dict()
  results["tvals"] = np.array(tvals)
  results["qvals"] = np.array(qvals)
  results["iofqt"] = iofqt
  return results






def extract_iofqt(directory):

  # identify the files of the directory format
  timefile  = "%s/initial_time.txt" % (directory)
  pixelfile = "%s/xcoordinates.txt" % (directory)
  datafile  = "%s/initial_intensity.txt" % (directory)

  # read the files in the directory
  tvals = np.loadtxt(timefile)
  xvals = np.loadtxt(pixelfile)
  Rvals = np.transpose(np.loadtxt(datafile))

  # identify the cetner and window width for each reading
  min_width = len(xvals)
  max_coords = np.zeros(len(tvals))
  for ii,R in enumerate(Rvals):
    max_coords[ii] = np.argmax(R)
    min_width = min([min_width, max_coords[ii], len(xvals)-max_coords[ii]])
  min_width -= 2
  min_width = int(min_width)

  # build the centered, adjusted R(q,t)
  iofqt = np.zeros((len(tvals), 2*min_width+1))
  for ii,R in enumerate(Rvals):
    xmin = max_coords[ii]-min_width
    xmax = max_coords[ii]+min_width
    iofqt[ii,:] = R[xmin:xmax+1]
    
  # build the qvals
  window = np.arange(-min_width, min_width+1)
  qvals = (2*3.14159*np.sin(3.14159*(window)/(131*180))/.124)  # in 1/nm

  # put the results into a dictionary and return it
  results = dict()
  results["tvals"] = np.array(tvals)
  results["qvals"] = np.array(qvals)
  results["iofqt"] = iofqt
  return results


