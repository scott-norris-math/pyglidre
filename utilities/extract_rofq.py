import sys
import copy
import numpy as np
import matplotlib as mpl ; mpl.use('Agg') ;
import matplotlib.pyplot as plt
from pygisaxs.utilities.NestedLmfit import *



# -----------------------------------------
#    Local fit routine
# -----------------------------------------


class PSeries(object):

  def __init__(self):
    self.values  = list()
    self.stderrs = list()




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


# define objective function: returns the array to be minimized
def ode1sol_plus_constant(params, t):

  R  = params['R'].value
  I0 = params['I0'].value
  B  = params['beta'].value
  BG = params['background'].value

  model = BG + I0*np.exp(2.0*R*t) + B/(2.0*R)*(np.exp(2.0*R*t) - 1.0)
  return model




def fit_rofq_hierarchically(mydict, fnc2min, localguess, globalguess, t_indices=None, q_indices=None, globalfit_q_indices=None, image_directory=None):

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

  # check for NaN values
  if np.sum(np.isnan(iofqt)) > 0:
    print "WARNING:  in fit_rofq_hierarchically(), iofqt contains NaN values."
    print "This will probably severely disrupt my ability to fit the data."
    print "Check prior manipulations, or use q_indices to avoid these values."


  # perform the global fits only over q in globalfit_q_indices
  # ------------------------------------------------------
  if globalguess != None and globalfit_q_indices != None:
    print "global fits over globalfit_q_indices ... " ; sys.stdout.flush()

    nodelist = list()
    for k in globalfit_q_indices:
      # identify data to fit
      trange = tvals[t_indices] - tvals[t_indices[0]]
      Ivals  = np.array(iofqt[t_indices,k])
      guess  = copy.deepcopy(localguess)
      node   = FitNode(guess, fnc2min, trange, Ivals)
      nodelist.append(node)

    sn = FitSuperNode(nodelist, copy.deepcopy(globalguess))
    globalresult = sn.fit()
    print globalresult.params

    nodelist = list()
    for k in q_indices:
      # identify data to fit
      trange = tvals[t_indices]
      Ivals  = np.array(iofqt[t_indices,k])
      guess  = copy.deepcopy(localguess)
      node   = FitNode(guess, fnc2min, trange, Ivals)
      nodelist.append(node)

    for node in nodelist:
      for key in globalguess:    
        node.params[key].value = globalresult.params[key].value
        node.params[key].stderr = globalresult.params[key].stderr
        node.params[key].set(vary=False)

      node.fit()


  elif globalguess != None and globalfit_q_indices == None:
    print "global fits over q_indices ... " ; sys.stdout.flush()
    nodelist = list()
    for k in q_indices:
      # identify data to fit
      trange = tvals[t_indices]
      Ivals  = np.array(iofqt[t_indices,k])
      guess  = copy.deepcopy(localguess)
      node   = FitNode(guess, fnc2min, trange, Ivals)
      nodelist.append(node)

    sn = FitSuperNode(nodelist, copy.deepcopy(globalguess))
    globalresult = sn.fit()
    print globalresult.params



  else:
    print "local fits over q_indices ... " ; sys.stdout.flush()
    nodelist = list()
    for k in q_indices:
      # identify data to fit
      trange = tvals[t_indices]
      Ivals  = np.array(iofqt[t_indices,k])
      guess  = copy.deepcopy(localguess)
      node   = FitNode(guess, fnc2min, trange, Ivals)
      nodelist.append(node)
      node.fit()











  goodqi = []
  goodnodes = []
  for k,node in zip(q_indices, nodelist):
    if (True):
    #if node.result.errorbars == True:
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

    k0 = np.where(qvals == 0.0)[0]
    for k,node in zip(q_indices, nodelist):

      if np.mod(k-k0,3) != 0:
        continue

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
      plottxt = plottxt + "\n"
      plottxt = plottxt + "nfev = %r" % (result.nfev)
      plottxt = plottxt + "success = %r" % (result.success)
      plottxt = plottxt + "errorbars = %r" % (result.errorbars)
      print result.message


      # make the plot
      p = plt.figure()
      ax = p.add_subplot(111)
      plt.plot(trange, Ivals, 'ro')
      finet = np.linspace(trange[0], trange[-1], 201)
      plt.plot(finet, fnc2min(result.params, finet), 'b')
      plt.gca().set_ylim(bottom=0)
      plt.title('time series fit for q=%f' % (q))
      plt.text(0.5, 0.2, plottxt, transform=ax.transAxes)
      plt.savefig('%s/fitted-I(t)-pixel=%03d.png' % (image_directory, k))
      plt.close()

    print "done." ; sys.stdout.flush()

  return glob




