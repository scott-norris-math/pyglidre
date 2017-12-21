import copy
import numpy as np
import scipy.interpolate as interp



# -----------------------------------------------
#    shift by co-ordinate of max-valued pixel
# -----------------------------------------------


def shift_by_max_value_1D(mydict):
  pass


def shift_by_max_value_2D(mydict):
  pass




# -----------------------------------------------
#    shift by location of center of mass
# -----------------------------------------------


def shift_by_center_mass_2D(mydict, peakradius=0.1):

  # unpack storage object
  tvals  = mydict["tvals"]
  qxvals = mydict["qvals"]
  qzvals = mydict["qzvals"]
  iofqt  = mydict["iofqt"]

  [QX, QZ] = np.meshgrid(qxvals, qzvals, indexing='ij')

  # walk through each time, correcting
  new = np.empty(np.shape(iofqt))
  for kk,oldrow in enumerate(iofqt):

    # find the peak
    peak_coords  = ( QX**2 + QZ**2 < peakradius**2 )

    # find the center of mass
    peak_M0  = np.sum(peak_coords * iofqt[kk,:,:]**2 )
    peak_M1x = np.sum(peak_coords * iofqt[kk,:,:]**2 * QX )
    peak_M1z = np.sum(peak_coords * iofqt[kk,:,:]**2 * QZ )
    x0   = peak_M1x / peak_M0
    z0   = peak_M1z / peak_M0

    print "at time t=%f, x0=%f, z0=%f" % (tvals[kk], x0, z0)

    # adjust the value of the data by the fit for the peak
    newrow = interp.griddata(((QX - x0).ravel(), (QZ - z0).ravel()), oldrow.ravel(), (QX, QZ), method='linear', fill_value=0.0)
    new[kk,:] = newrow

  # return a *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt'] = new
  return newdict



def shift_by_center_mass_1D(mydict, peakradius=0.1):

  # unpack storage object
  tvals = mydict["tvals"]
  qvals = mydict["qvals"]
  iofqt = mydict["iofqt"]

  # walk through each time, correcting
  new = np.empty(np.shape(iofqt))
  for kk,oldrow in enumerate(iofqt):

    # find the center of mass of the peak
    peak_coords  = np.where(np.abs(qvals) < peakradius)[0]
    peak_M0 = np.sum(iofqt[kk,peak_coords]**2 )
    peak_M1 = np.sum(iofqt[kk,peak_coords]**2 * qvals[:,peak_coords] )
    cmass   = peak_M1 / peak_M0

    print "at time t=%f, cmass=%f" % (tvals[kk], cmass)

    # adjust the value of the data by the fit for the peak
    newrow = interp.griddata(qvals - cmass, oldrow, qvals, method='linear', fill_value=0.0)
    new[kk,:] = newrow

  # return a *copy* of the data
  newdict = copy.deepcopy(mydict)
  newdict['iofqt'] = new
  return newdict









# -----------------------------------------------
#    shift by center of fit to peak shape
# -----------------------------------------------


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


def shit_by_peak_fit_1D(mydict, directory_name=None):

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
    newrow = interp.griddata(qvals - Q0, oldrow, qvals, method='linear', fill_value=0.0) / A

    # report and store
    print "A = %2.6f; C=%2.6f; C/A=%2.6f; Q0=%2.6f.   oldpeak=%f.   newpeak=%f" % (A, C, C/A, Q0, oldrow[qmid], newrow[qmid])
    new[kk,:] = newrow

  newdict = dict()
  newdict['tvals'] = mydict['tvals']
  newdict['qvals'] = mydict['qvals']
  newdict['iofqt'] = new
  return newdict

