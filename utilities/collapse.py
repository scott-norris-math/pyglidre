
import copy 
import numpy as np



#  ---------------------------------------
#    convert 2D to 1D
#  ---------------------------------------
def collapse_2D_to_1D(mydict, vspan=9, qmin=-np.inf, qmax=np.inf):
  '''
  # This function converts a 2D GISAXS reading to a 1D vertical average centered on the specular peak.
  #
  # In general the function first applies hot-pixel removal to deal with slightly faulty detectors, and then
  # finds the brightest pixel in the image and calls it the specular peak.  If you know the specular peak will
  # be the brightest pixel, then you can just let the function find it for you (and it does so on an 
  # image-by-image basis, so it automatically handles slight jitter between readings).
  #
  # However, if the specular peak is **not** the brightest pixel in the image, you can specify a smaller 
  # sub-window in which to look for it, by listing four numbers in a file called "specular-peak-window.txt"
  '''

  # unpack
  qxvals = mydict["qxvals"]
  qzvals = mydict["qzvals"]
  iofqt  = mydict["iofqt"]

  # find qx range
  qi     = np.where( (qxvals > qmin) * (qxvals < qmax) )[0]
  qimin  = qi[0]
  qimax  = qi[-1]

  # find z location of the specular peak, identify window
  z0     = np.where(qzvals == 0.0)[0][0]
  zmin   = z0-(vspan-1)/2
  zmax   = z0+(vspan-1)/2

  # average over window
  new_qxvals = qxvals[qimin:qimax]
  new_iofqt = np.sum(iofqt[:,qimin:qimax,zmin:zmax+1], axis=2)

  # return copy of data 
  newdict = copy.deepcopy(mydict)
  newdict["qxvals"] = new_qxvals
  newdict["qvals"] = new_qxvals
  newdict["iofqt"] = new_iofqt
  return newdict









