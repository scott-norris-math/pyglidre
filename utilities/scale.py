import copy
import numpy as np
import matplotlib as mpl ; mpl.use('Agg') ;
import matplotlib.pyplot as plt

# This library needs to be installed using pip
from lmfit import Parameter, Parameters, minimize


'''
All of these functions, either 1D or 2D, can assume that the data
has been pre-aligned using a shift operator.  The recieved data
structure will indecate co-ordinates of the specular peak (where q=0),
and the peak may be assumed to be centered in the *middle* of that pixel.
'''



# ------------------------------------------------------
#    scale by the instrument-provided x-ray fluence
# ------------------------------------------------------


def scale_by_fluence_2D(mydict):

  # unpack storage object
  iofqt   = mydict["iofqt"]
  fluence = mydict["fluences"]

  # sanity check
  if (fluence==None):
    print "WARNING: You have requested to scale by instrument-recorded fluence."
    print "However, your data structure does not include a fluence reading."
    print "I am returning the data without any scaling applied."
    return mydict

  # divide each intensity function by the instrument-provided fluence
  new = np.empty(np.shape(iofqt))
  for kk,olddata in enumerate(iofqt):
    new[kk,:] = olddata / fluence[kk]

  # return an updated *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt'] = new
  return newdict




def scale_by_fluence_1D(mydict):

  # unpack storage object
  iofqt   = mydict["iofqt"]
  fluence = mydict["fluences"]

  # sanity check
  if (fluence==None):
    print "WARNING: You have requested to scale by instrument-recorded fluence."
    print "However, your data structure does not include a fluence reading."
    print "I am returning the data without any scaling applied."
    return mydict

  # divide each intensity function by the instrument-provided fluence
  new = np.empty(np.shape(iofqt))
  for kk,olddata in enumerate(iofqt):
    new[kk,:] = olddata / fluence[kk]

  # return an updated *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt'] = new
  return newdict





# ----------------------------------------------------------
#    scale the data by the total photon count
# ----------------------------------------------------------


def scale_by_all_2D(mydict):

  # unpack storage object
  iofqt2D = mydict["iofqt2D"]

  # at each time, count the photons in the image
  total_photons = np.sum(iofqt, axis=(1,2))
  
  # divide each intensity function by the total number of photons recorded
  new = np.empty(np.shape(iofqt))
  for kk,olddata in enumerate(iofqt):
    new[kk,:] = olddata / total_photons[kk]

  # return an updated *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt2D'] = new
  return newdict



def scale_by_all_1D(mydict, radius=1.0):

  # unpack storage object
  iofqt = mydict["iofqt"]

  # at each time, count the photons in the image
  total_photons = np.sum(iofqt, axis=1)              ###  FOR SOME REASON THIS IS RETURNING NAN!!!

  # divide each intensity function by the total number of photons recorded
  new = np.empty(np.shape(iofqt))
  for kk,olddata in enumerate(iofqt):
    new[kk,:] = olddata / total_photons[kk]

  # return an updated *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt'] = new
  return newdict




# ---------------------------------------------------------------------
#    scale the data by the photon count in the peak (assume constant)
# ---------------------------------------------------------------------




def scale_by_peak_1D(mydict, peakradius=0.02):

  # unpack storage object
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  # at each time, count the photons in the specular peak
  peak_coords  = np.where(np.abs(qvals) < peakradius)[0]
  peak_photons = np.sum(iofqt[:,peak_coords], axis=1)

  # divide each intensity function by the number of photons in the peak
  new = np.zeros(np.shape(iofqt))
  for kk,oldrow in enumerate(iofqt):
    new[kk,:] = oldrow / peak_photons[kk]

  # return an updated *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt'] = new
  return newdict




# ---------------------------------------------------------------------
#    scale the data by the photon count in the peak (fit to exponential)
# ---------------------------------------------------------------------


def scale_by_fitted_peak_1D(mydict, peakradius=0.0333):

  # unpack storage object
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  # at each time, count the photons in the specular peak
  peak_coords  = np.where(np.abs(qvals) < peakradius)[0]
  peak_photons = np.sum(iofqt[:,peak_coords], axis=1)

  # build a Parameters object for a stretched exponential
  params = Parameters()
  params.add('A', value= 1e-6, min=0.0)
  params.add('B', value= 0.00)
  params.add('C', value=-0.02, max=0.0)
  params.add('I', value= peak_photons[0], min=0.0)

  #fit the function 
  result   = minimize(ode2sol, params, args=(tvals, peak_photons), method="nelder")
  result   = minimize(ode2sol, params, args=(tvals, peak_photons), method="leastsq")
  smoothed = peak_photons  + ode2sol(params, tvals, tail_photons)

  # divide each intensity function by the number of photons in the peak
  new = np.zeros(np.shape(iofqt))
  for kk,oldrow in enumerate(iofqt):
    new[kk,:] = oldrow * (smoothed[kk] / peak_photons[kk]) / np.mean(smoothed)

  # return an updated *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt'] = new
  return newdict










# ---------------------------------------------------------------------
#    scale the data by the photon count in the tails (assume constant)
# ---------------------------------------------------------------------


def ode2sol(params, t, data):
  '''
  solution to the 1st-order quadratic ODE 
  
     dx/dt = A + B x + C x^2 
      x(0) = I
  
  useful when something is almost exponential.
  '''

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





def scale_by_tails_1D(mydict, tailradius=0.1):

  # unpack storage object
  tvals = mydict["tvals"]
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  # at each time, count the photons in the tails
  tail_coords = np.where(np.abs(qvals > tailradius))[0]
  tail_photons = np.sum(iofqt[:,tail_coords], axis=1)

  Iguess = tail_photons[0]

  # build a Parameters object
  params = Parameters()
  params.add('A', value= 1e-6, min=0.0)
  params.add('B', value= 0.00)
  params.add('C', value=-0.02, max=0.0)
  params.add('I', value= Iguess, min=0.0)

  #fit the function 
  result   = minimize(ode2sol, params, args=(tvals, tail_photons), method="nelder")
  result   = minimize(ode2sol, params, args=(tvals, tail_photons), method="leastsq")
  smoothed = tail_photons  + ode2sol(params, tvals, tail_photons)

  # scale the data by the fitting
  new = np.zeros(np.shape(iofqt))
  for kk,oldrow in enumerate(iofqt):
    new[kk,:] = oldrow / tail_photons[kk] * smoothed[kk]

  # return an updated *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt'] = new
  return newdict


  #p = plt.figure()
  #ax = p.add_subplot(111)
  #plt.plot(tvals, all_photons, 'bs')
  #plt.plot(tvals, smoothed, 'r')
  #plt.text(0.5, 0.8, 
  #  'A=%.3e \n B=%.3e \n C=%.3e \n I=%.3e' % (params['A'].value, params['B'].value, params['C'].value, params['I'].value)
  #   , transform=ax.transAxes)
  #plt.savefig("%s/smoothing-fit.png" % (directory_name))
  #plt.close()





# ---------------------------------------------------------------------
#    identify and substract extra constant intensity by looking at the tails
# ---------------------------------------------------------------------


def subtract_vampires_1D(mydict, tailradius=0.75, q_indices=None, directory='.'):

  # unpack storage object
  tvals = mydict["tvals"]
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  if q_indices==None:  q_indices = qvals

  # at each time, count the photons in the tails
  tail_coords = np.where(np.abs(qvals) > tailradius)[0]
  tail_coords = [i for i in tail_coords if i in q_indices]
  tail_photons = np.mean(iofqt[:,tail_coords], axis=1)
  Iguess = tail_photons[0]

  # build a Parameters object
  params = Parameters()
  params.add('A', value= 1e-6, min=0.0, max=np.inf)
  params.add('B', value= -1e-8)
  params.add('C', value= -0.02, min=-np.inf, max=-0.00001)
  params.add('I', value= Iguess, min=0.0, max=np.inf)
  #params.add('I', value= Iguess, vary=False)

  #fit the function 
  #result   = minimize(ode2sol, params, args=(tvals, tail_photons), method="nelder")
  result   = minimize(ode2sol, params, args=(tvals, tail_photons), method="leastsq", **{'xtol':1e-8})
  smoothed = tail_photons  + ode2sol(params, tvals, tail_photons)
  vampires = tail_photons - smoothed

  #vampires -=  np.min(vampires)  # ensure that we are not *adding* to the intensity anywhere

  #p = result.params
  #A = p["A"].value
  #B = p["B"].value
  #C = p["C"].value
  #longtime = (-B + np.sqrt(B**2-4*A*C)) / (2*C)

  # scale the data by the fitting
  new = np.zeros(np.shape(iofqt))
  for kk,oldrow in enumerate(iofqt):
    new[kk,:] = oldrow - vampires[kk]

  # diagnostic plot
  p = plt.figure()
  ax = p.add_subplot(111)
  plt.plot(tvals, tail_photons, 'bs')
  plt.plot(tvals, smoothed, 'r')
  plt.text(0.5, 0.8, 
    'A=%.3e \n B=%.3e \n C=%.3e \n I=%.3e' % (params['A'].value, params['B'].value, params['C'].value, params['I'].value)
     , transform=ax.transAxes)
  plt.savefig("%s/vampire-tail-fit.png"%(directory))
  plt.close()


  # return an updated *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt'] = new
  return newdict






