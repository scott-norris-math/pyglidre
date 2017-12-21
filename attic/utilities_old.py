'''
This file converts the greyscale TIF "images" produced by PILATUS 
into heatmap-colored .png images for further 
'''
import os
import sys
import copy
import shelve
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filters
import scipy.interpolate as interp

# This library needs to be installed using pip
from lmfit import minimize, Parameters, Parameter, report_fit
from NestedLmfit import *

# Another file from PyGisaxs
import pygisaxs.io as io







#  ---------------------------------------
#    convert 2D to 1D
#  ---------------------------------------
def convert_2D_to_1D(array2Dlist, average_vspan=None, search_window=None):
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

  # get dimensions of (t, qx, qz) array
  datadims = np.shape(array2Dlist)
  tn = datadims[0]
  xn = datadims[1]
  zn = datadims[2]

  # obtain a median-filtered version of array (eliminates broken pixels -- borrowed from StackExchange)
  mf_data = np.array(array2Dlist)
  for kk in xrange(tn):
    arr = array2Dlist[kk,:,:]						
    blurred    = filters.median_filter(arr, size=2)
    difference = arr - blurred
    threshold  = 3*np.std(difference)  					# too large and you don't remove pixels -- too small and you smooth
    hot_pixels = np.nonzero((np.abs(difference[1:-1,1:-1])>threshold) ) # ignore edges for simplicity; we'll discard anyway
    hot_pixels = np.array(hot_pixels) + 1 				# because we ignored the first row and first column
    for x,y in zip(hot_pixels[0], hot_pixels[1]):
      mf_data[kk,x,y]=blurred[x,y]

  # build the search window
  if search_window == None:
    print "Looking over the whole image for the specular peak."
    cxmin = 0
    cxmax = xn
    cymin = 0
    cymax = zn
  else:
    print "Looking for specular peak in provided window ", center_window
    cxmin = search_window[0]
    cxmax = search_window[1]
    cymin = search_window[2]
    cymax = search_window[3]

  # find the specular peak
  totalI = sum(mf_data, axis=0)
  cwindow = totalI[cxmin:cxmax,cymin:cymax]
  temp_coords = np.unravel_index(np.argmax(cwindow), cwindow.shape)
  max_coords = np.array([cxmin, cymin]) + temp_coords
  x0 = max_coords[0]
  z0 = max_coords[1]

  # report, with optional warning
  print "Peak found at ", max_coords
  if x0 == cxmin or x0 == cxmax or z0 == cymin or z0==cymax:
    print "WARNING: specular peak found on edge of search window."
    print "Cannot confirm that value is a local maximum."
    print "Identified co-ordinates may be inaccurate."


  # identify the vertical analysis window (specular peak +- half the width/height specified in window_shape)
  zmin  = z0-(average_vspan-1)/2
  zmax  = z0+(average_vspan-1)/2
  iofqt = sum(mf_data[:,:,zmin:zmax+1], axis=2)
  return iofqt

















# -------------------------------------------------
#      data-smoothing functions
# -------------------------------------------------







# fit the specular peak for sub-pixel smoothing 
def simple_quadratic(params, qvals, data):

  A  = params['A'].value
  C  = params['C'].value
  Q0 = params['Q0'].value
  model = A + C*(qvals - Q0)**2
  return (model - data)
  

def simple_gaussian(params, qvals, data):
  Q0 = params['Q0'].value
  A  = params['A'].value
  C  = params['C'].value
  model = A * np.exp(C*(qvals - Q0)**2)
  return (model - data)


def exp2fun(params, t, data):

  A = params['A'].value
  B = params['B'].value
  C = params['C'].value
  I = params['I'].value

  Drt = np.sqrt(B**2-4*A*C)	#; print "Drt = ", Drt
  up  = (-B - Drt) / (2*C)	#; print "up = ", up
  um  = (-B + Drt) / (2*C)	#; print "um = ", um
  F   = (I-up)/(I-um)		#; print "I =", I
  model = (up - um*F*np.exp(-Drt*t)) / (1 - F*np.exp(-Drt*t))
  return (model - data)



def center_iofqt_by_peak(mydict, directory_name=None):

  # unpack storage object
  tvals = mydict["tvals"]
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  # get shape of the array
  ss = np.shape(iofqt)
  tn = ss[0]
  qn = ss[1]

  # area of the specular peak
  qmid = (qn-1)/2 
  qmin = qmid - 5
  qmax = qmid + 5

  # storage for corrected iofqt
  new = np.zeros(np.shape(iofqt))

  # walk through each time, correcting
  for kk,oldrow in enumerate(iofqt):

    # construct a guess for the specular peak
    opeak = oldrow[qmid]

    #fnc2min = simple_gaussian
    #myguess = Parameters()
    #myguess.add('A',  value=  opeak, min=0.0)
    #myguess.add('B',  value= -8000, max=0.0)
    #myguess.add('Q0', value=  0.0)

    fnc2min = simple_gaussian
    myguess = Parameters()
    myguess.add('A',  value=  opeak, min=0.0)
    myguess.add('C',  value= -10*opeak, max=0.0)
    myguess.add('Q0', value=  0.0)


    # fit the specular peak
    results1 = minimize(fnc2min, myguess, args=(qvals[qmin:qmax], oldrow[qmin:qmax]), method="nelder")
    results2 = minimize(fnc2min, results1.params, args=(qvals[qmin:qmax], oldrow[qmin:qmax]), **{"xtol":1e-12})
    rparams = results2.params
    A  = rparams["A"]
    C  = rparams["C"]
    Q0 = rparams["Q0"]

    # make a plot of the fit
    plt.figure()
    plt.plot(qvals[qmin:qmax], oldrow[qmin:qmax], 'bs')
    plt.plot(qvals[qmin:qmax], oldrow[qmin:qmax] + fnc2min(rparams, qvals[qmin:qmax], oldrow[qmin:qmax]))
    plt.xlabel('q')
    plt.ylabel('iofq')
    plt.title('specular peak fit for t=%f' %(tvals[kk]))
    plt.savefig("%s/specular-peak-time=%f.png" % (directory_name, tvals[kk]))

    # adjust the value of the data by the fit for the peak
    newrow = interp.griddata(qvals - Q0, oldrow, qvals, method='linear') / A

    # report and store
    print "A = %2.6f; C=%2.6f; C/A=%2.6f; Q0=%2.6f.   oldpeak=%f.   newpeak=%f" % (A, C, C/A, Q0, oldrow[qmid], newrow[qmid])
    new[kk,:] = newrow

  newdict = dict()
  newdict['tvals'] = mydict['tvals']
  newdict['qvals'] = mydict['qvals']
  newdict['iofqt'] = new
  return newdict






def smooth_iofqt_by_peak(mydict, directory_name=None, peakradius=3):

  # unpack storage object
  tvals = mydict["tvals"]
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  # get shape of the array
  ss = np.shape(iofqt)
  tn = ss[0]
  qn = ss[1]

  # at each time, count the photons in the specular peak (assume it is 9 pixels wide)
  peakcenter = np.where(qvals==0.0)[0][0]
  qlo = peakcenter-peakradius
  qhi = peakcenter+peakradius
  peak_photons = np.sum(iofqt[:,qlo:qhi+1], axis=1)
  
  # divide each intensity function by the number of photons in the peak
  new = np.zeros(np.shape(iofqt))
  for kk,oldrow in enumerate(iofqt):
    new[kk,:] = oldrow / peak_photons[kk]

  newdict = dict()
  newdict['tvals'] = mydict['tvals']
  newdict['qvals'] = mydict['qvals']
  newdict['iofqt'] = new
  return newdict






def smooth_iofqt_by_fitted_wings(mydict, directory_name=None):

  tvals = mydict["tvals"]
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  ss = np.shape(iofqt)
  tn = ss[0]
  qn = ss[1]

  q40 = (qn-1)/2-12#np.percentile(range(qn), 40)
  q60 = (qn-1)/2+12#np.percentile(range(qn), 60)
  all_photons = np.sum(iofqt[:,0:q40], axis=1) + np.sum(iofqt[:,q60:], axis=1)

  # build a Parameters object
  params = Parameters()
  params.add('A', value= 1e-6, min=0.0)
  params.add('B', value= 0.00)
  params.add('C', value=-0.02, max=0.0)
  params.add('I', value= 0.02, min=0.0)
  result   = minimize(exp2fun, params, args=(tvals, all_photons), method="nelder")
  result   = minimize(exp2fun, params, args=(tvals, all_photons))
  smoothed = all_photons + exp2fun(params, tvals, all_photons)
  #print smoothed 

  p = plt.figure()
  ax = p.add_subplot(111)
  plt.plot(tvals, all_photons, 'bs')
  plt.plot(tvals, smoothed, 'r')
  plt.text(0.5, 0.8, 
    'A=%.3e \n B=%.3e \n C=%.3e \n I=%.3e' % (params['A'].value, params['B'].value, params['C'].value, params['I'].value)
     , transform=ax.transAxes)
  #plt.show()
  plt.savefig("%s/smoothing-fit.png" % (directory_name))
  plt.close()

  #blurred = filters.median_filter(all_photons, size=3)
  #difference = all_photons - blurred
  #threshold = 2*np.std(difference)
  #hot_pixels = np.nonzero((np.abs(difference[1:-1])>threshold) )  # ignore edges
  #hot_pixels = np.array(hot_pixels) + 1 # because we ignored the first row and first column
  #farr       = np.copy(all_photons)     # This will be the image with the hot pixels removed
  #for x in hot_pixels:
  #  farr[x]=blurred[x]

  new = np.zeros(np.shape(iofqt))
  for kk,oldrow in enumerate(iofqt):
    new[kk,:] = oldrow / all_photons[kk] * smoothed[kk]

  newdict = dict()
  newdict['tvals'] = mydict['tvals']
  newdict['qvals'] = mydict['qvals']
  newdict['iofqt'] = new
  return newdict









# -----------------------------------------
#    Local fit routine
# -----------------------------------------



# define objective function: returns the array to be minimized
def expfun(params, t, data):

  R  = params['R'].value
  B  = params['beta'].value
  I0 = params['I0'].value
  model = I0*np.exp(2.0*R*t) + B/(2.0*R)*(np.exp(2.0*R*t) - 1.0)
  return (model - data)


def quadfun(params, t, data):

  R  = params['R'].value
  B  = params['beta'].value
  I0 = params['I0'].value
  model = I0 + (B+I0*2.0*R)*t + (2.0*R)/2.0*(B+I0*2.0*R)*t**2 #+ (2.0*R)**2/6.0*(B+I0*2.0*R)*t**3
  return model - data





class PSeries(object):

  def __init__(self):
    self.values  = list()
    self.stderrs = list()

  
def is_good_fit(fit):
  return fit.errorbars



def fit_rofq_locally(mydict, fnc2min, myguess, trange=[None, None], qrange=[None, None], image_directory=None):
  '''
  This function fits each timeseries in I(q,t) to extract R(q).
  It does so using a standard fit, where parameters are optimized
  *locally* (they can vary for each q).  It is useful for getting
  ideas about how to construct a hierarchical fit, which seems to
  produce better results, especially for nearly-constant timeseries.
  '''

  # pull data out of dictionary for ease of reference
  tvals = mydict["tvals"]
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  # get shape of data
  tn   = len(tvals)
  qn   = len(qvals)
  qmid = (qn-1.0)/2.0

  # user-specified time range
  tmin = 0    if trange[0]==None else trange[0]
  tmax = tn   if trange[1]==None else trange[1]

  # user-specified qval range
  qmin = 0    if qrange[0]==None else qrange[0]
  qmax = qn   if qrange[1]==None else qrange[1]


  # perform the initial fit attempt for each time sequence
  # ------------------------------------------------------

  print "initial fits ... ", ; sys.stdout.flush()
  totalvals = 0
  goodvals = 0

  fit_results = list()
  for k,q in enumerate(qvals[qmin:qmax]):

    # identify data to fit
    Ivals = iofqt[tmin:tmax,k]
    trange = tvals[tmin:tmax]

    # call function minimizer
    eparams = myguess
    result1 = minimize(fnc2min, eparams, args=(trange, Ivals), method="nelder")
    result2 = minimize(fnc2min, result1.params, args=(trange, Ivals), **{"xtol":1e-12})
    fit_results.append(result2)

    totalvals += 1
    if is_good_fit(result2): goodvals += 1

  print "%d / %d good fits." % (goodvals, totalvals) ; sys.stdout.flush()


  # now try to fill in holes using fitted values for *left* neighbors
  # -----------------------------------------------------------

  print "attempting fill-ins from left  ... ", ; sys.stdout.flush()
  fixed = 0

  for pp in range(qn)[qmin+1:qmax]:

    kk   = pp - qmin
    fit  = fit_results[kk]
    lfit = fit_results[kk-1]

    if (is_good_fit(fit) == True):
      #print "q=%0.5f is good." % (qvals[kk])
      continue ;  

    else:
      #print "q=%0.5f is poor." % (qvals[kk]),

      if (is_good_fit(lfit) == False):
        #print "Left neighbor is also bad -- giving up."
        continue

      else:
        #print "Left neighbor is good fit -- trying that.", 
        trange = tvals[tmin:tmax]
        Ivals = iofqt[tmin:tmax,kk]
        result1 = minimize(fnc2min, lfit.params, args=(trange, Ivals), method="nelder")
        result2 = minimize(fnc2min, result1.params, args=(trange, Ivals), **{"xtol":1e-12})
        fit_results[kk] = result2

        if (is_good_fit(fit_results[kk]) == True): 
          fixed += 1
          #print "success!"
        #else: 
          #print "still bad -- giving up."

  print "fixed %d qvals" % (fixed) ; sys.stdout.flush()


      

  # now try to fill in holes using fitted values for *right* neighbors
  # -----------------------------------------------------------

  print "attempting fill-ins from right ... ", ; sys.stdout.flush()
  fixed = 0

  for pp in reversed(range(qn)[qmin:qmax-1]):

    kk   = pp - qmin
    fit  = fit_results[kk]
    rfit = fit_results[kk+1]

    if (is_good_fit(fit) == True):
      #print "q=%0.5f is good." % (qvals[kk])
      continue ;  

    else:
      #print "q=%0.5f is poor." % (qvals[kk]), 

      if (is_good_fit(rfit) == False):
        #print "Right neighbor is also bad -- giving up."
        continue

      else:
        #print "Right neighbor is good fit -- trying that.", 
        trange = tvals[tmin:tmax]
        Ivals = iofqt[tmin:tmax,kk]
        result1 = minimize(fnc2min, rfit.params, args=(trange, Ivals), method="nelder")
        result2 = minimize(fnc2min, result1.params, args=(trange, Ivals), **{"xtol":1e-12})
        fit_results[kk] = result2

        if (is_good_fit(fit_results[kk]) == True): 
          fixed += 1
          #print "success!"
        #else: 
          #print "still bad -- giving up."

  print "fixed %d qvals" % (fixed) ; sys.stdout.flush()


  # convert a list of results objects to a dictionary of functions of q
  # ------------------------------------------------------

  # storage allocation
  glob = dict()
  glob["qvals"] = PSeries()
  glob["qvals"].values = qvals[qmin:qmax]
  for key in myguess.keys():
    glob[key] = PSeries()
  glob["chisqr"] = PSeries()

  # extract data in the list of results
  for k,result in enumerate(fit_results):
    for key in myguess.keys():
      pfit = result.params[key]
      glob[key].values.append(pfit.value)
      glob[key].stderrs.append(pfit.stderr)

    glob["chisqr"].values.append(result.chisqr)


  # create a diagnostic plot of each time series and best fit.
  # ------------------------------------------------------
  if image_directory != None:
    print "generating images ... ", ; sys.stdout.flush()

    for k,q,result in zip(range(qn)[qmin:qmax], qvals[qmin:qmax], fit_results):

      # identify data to fit
      Ivals = iofqt[tmin:tmax,k]
      trange = tvals[tmin:tmax]

      # text to write on plot
      plottxt = ''
      for key in myguess.keys():
        pfit = result.params[key]
        plottxt = plottxt + "%s = %.3e +- %.3e \n" % (key, pfit.value, pfit.stderr)
      plottxt = plottxt + "chisqr = %.3e \n" % (result.chisqr)
      plottxt = plottxt + "errorbars = %r" % (result.errorbars)


      # make the plot
      p = plt.figure()
      ax = p.add_subplot(111)
      plt.plot(trange, Ivals, 'ro')
      plt.plot(trange, Ivals + fnc2min(result.params, trange, Ivals), 'b')
      plt.gca().set_ylim(bottom=0)
      plt.title('time series fit for q=%f' % (q))
      plt.text(0.5, 0.2, plottxt, transform=ax.transAxes)
      plt.savefig('%s/fitted-I(t)-pixel=%03d.png' % (image_directory, k))
      plt.close()

    print "done." ; sys.stdout.flush()

  return glob








# -------------------------------------------------
#    hierarchical fit routine
# -------------------------------------------------

# define objective function: returns the array to be minimized
def ode1sol(params, t):

  R  = params['R'].value
  B  = params['beta'].value
  I0 = params['I0'].value
  model = I0*np.exp(2.0*R*t) + B/(2.0*R)*(np.exp(2.0*R*t) - 1.0)
  return model



def ode1sol_fitted_time_scalings(params, t):

  

  pass



def fit_rofq_hierarchically(mydict, fnc2min, localguess, globalguess, t_indices=None, q_indices=None, image_directory=None):

  '''
  This function fits each timeseries in I(q,t) to extract R(q).
  It does so using a hierarchical fit, where some of the parameters
  are optimized *globally* (they must be the same for all q), 
  while others are optimized *locally* (they can vary for each q).
  '''

  # pull data out of dictionary for ease of reference
  tvals = mydict["tvals"]
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  # get shape of data
  tn   = len(tvals)
  qn   = len(qvals)

  # get indices over which to look
  if t_indices==None:  t_indices=range(tn)
  if q_indices==None:  q_indices=range(qn)




  # perform the initial fit attempt for each time sequence
  # ------------------------------------------------------

  print "global fits ... ", ; sys.stdout.flush()
  nodelist = list()
  for k in q_indices:

    q = qvals[k]

    # identify data to fit
    trange = tvals[t_indices]
    Ivals  = np.array(iofqt[t_indices,k])
    guess  = copy.deepcopy(localguess)
    node   = FitNode(guess, fnc2min, trange, Ivals)
    nodelist.append(node)

  if globalguess != None:
    sn = FitSuperNode(nodelist, copy.deepcopy(globalguess))
    globalresult = sn.fit()
    print globalresult.params
  else:
    for node in nodelist:  node.fit()

  goodqi = []
  goodnodes = []
  for k,node in zip(q_indices, nodelist):
    if node.result.errorbars == True:
      goodqi.append(k)
      goodnodes.append(node)
  goodqi = np.array(goodqi)
  

  # convert a list of results objects to a dictionary of functions of q
  # ------------------------------------------------------

  # storage allocation
  glob = dict()
  glob["qvals"] = PSeries()
  glob["qvals"].values = qvals[goodqi]
  for key in localguess.keys():
    glob[key] = PSeries()
  glob["chisqr"] = PSeries()

  # extract data in the list of results
  fit_results = [a.result for a in goodnodes]
  for k,result in enumerate(fit_results):
    for key in localguess.keys():
      pfit = result.params[key]
      glob[key].values.append(pfit.value)
      glob[key].stderrs.append(pfit.stderr)

    glob["chisqr"].values.append(result.chisqr)



  # create a diagnostic plot of each time series and best fit.
  # ------------------------------------------------------

  if image_directory != None:
    print "generating images ... ", ; sys.stdout.flush()

    for k,node in zip(q_indices, nodelist):

      q = qvals[k]

      # identify data to fit
      trange = node.xvals
      Ivals  = node.yvals
      result = node.result

      # text to write on plot
      plottxt = ''
      for key in localguess.keys():
        pfit = result.params[key]
        plottxt = plottxt + "%s = %.3e +- %.3e \n" % (key, pfit.value, pfit.stderr)
      plottxt = plottxt + "chisqr = %.3e \n" % (result.chisqr)
      plottxt = plottxt + "errorbars = %r" % (result.errorbars)


      # make the plot
      p = plt.figure()
      ax = p.add_subplot(111)
      plt.plot(trange, Ivals, 'ro')
      plt.plot(trange, fnc2min(result.params, trange), 'b')
      plt.gca().set_ylim(bottom=0)
      plt.title('time series fit for q=%f' % (q))
      plt.text(0.5, 0.2, plottxt, transform=ax.transAxes)
      plt.savefig('%s/fitted-I(t)-pixel=%03d.png' % (image_directory, k))
      plt.close()

    print "done." ; sys.stdout.flush()

  return glob








# -----------------------------------------------
#    plotting section
# -----------------------------------------------


def plot_iofqt_heatmap(tqi, filename=None):

  # unpack variables from dictionary
  tvals = tqi['tvals']
  qvals = tqi['qvals']
  iofqt = tqi['iofqt']

  # set the filename
  if filename==None:  filename="iofqt-heatmap.png"

  # Create a time/space heatmap of I(q,t).  
  #Works best if there are lots of time data points (continuous bombardment)
  plt.figure()
  plt.imshow(np.log10(tqi['iofqt']), origin="lower")
  plt.xlabel(r"q [nm$^{-1}$]")
  plt.ylabel(r"t [sec]")
  plt.title(r"Heatmap of $\log_{10}$(intensity)")
  plt.savefig(filename)
  plt.close()




def plot_intensity_progression(tqi, filename=None, max_lines=None):

  # unpack variables from dictionary
  tvals = tqi['tvals']
  qvals = tqi['qvals']
  iofqt = tqi['iofqt']

  # set the filename
  if filename==None:  filename="intensity-progression.png"

  # set the maximum number of lines in the plot, and get plot times 
  if max_lines==None:  max_lines = 20
  interval = np.ceil(len(tvals)/float(max_lines))
  tlist = np.arange(0,len(tvals),interval)

  # get list of colors
  cmap = plt.get_cmap('rainbow')
  colors = [cmap(i) for i in np.linspace(0, 1, len(tlist))]

  # Create a single plot with up to 20 I(q) values at different times
  plt.figure()
  for ii,tk in enumerate(tlist):
    plt.plot(qvals, np.log10(iofqt[tk,:]), color=colors[ii], label="t=%d"%( round(tvals[tk]) ))
  plt.xlabel(r"q [nm$^{-1}$]")
  plt.ylabel(r"$\log_{10}$(intensity)")
  plt.title(r"Selected values of $\log_{10}$(intensity)")
  plt.legend(prop={'size':8})
  plt.savefig(filename)
  plt.close()

