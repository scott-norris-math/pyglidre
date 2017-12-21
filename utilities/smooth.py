import numpy as np
import scipy.ndimage.filters as filters



def remove_hot_pixels_2d(array2d, stencil=2, threshold=3.0):
  '''
  obtain a median-filtered version of array (eliminates broken pixels -- borrowed from StackExchange)
  '''

  blurred    = filters.median_filter(array2d, size=stencil)
  difference = array2d - blurred
  threshold  = 3*np.std(difference)  						# too large and you don't remove pixels -- too small and you smooth
  hot_pixels = np.nonzero((np.abs(difference[1:-1,1:-1]) > threshold) ) 	# ignore edges for simplicity; we'll discard anyway
  print hot_pixels
  hot_pixels = np.array(hot_pixels) + 1 					# because we ignored the first row and first column
  array2dcp  = np.copy(array2d) 						# This is the image with the hot pixels removed

  for x,y in zip(hot_pixels[0], hot_pixels[1]):
    array2dcp[x,y] = blurred[x,y]

  return array2dcp

