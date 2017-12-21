import copy
import numpy as np
from lmfit import Parameter, Parameters, minimize


class FitNode(object):
  '''
  This class stores
  - a single lmfit.Parameters() object
  - a single supplied fitting function 
  - a list of indepdnent variable values
  - a list of dependent variable values
  - a list of uncertanty values (optional)
  and then optimizes the parameters using lmfit.minimize().

  Optionally, the standard errors of the data can be provided to guide the fit. 
  '''

  def __init__(self, params, ffunc, xvals, yvals, uvals=None):

    self.params = params
    self.ffunc  = ffunc
    self.xvals  = xvals
    self.yvals  = yvals
    self.uvals  = uvals
    self.result = None


  def evaluate(self, params=None, xvals=None):
    # evaluate the fitting function for a given set of independent variables (useful for plotting)

    if params == None:  params = self.result.params
    if xvals  == None:  xvals  = self.xvals
    return self.ffunc(params, xvals)


  def residual(self, params=None, xvals=None, yvals=None, uvals=None):
    # simple residual function with optional uncertainty weighting

    if params == None:  params = self.params
    if xvals  == None:  xvals  = self.xvals
    if yvals  == None:  yvals  = self.yvals
    if uvals  == None:  uvals  = self.uvals

    if uvals  == None:
      return (self.ffunc(params, xvals) - yvals)
    else:
      return (self.ffunc(params, xvals) - yvals) / uvals


  def fit(self, simplex_step=True):

    # make a copy of the guess, to avoid changing it
    guess = copy.deepcopy(self.params)

    # optional first simplex estimate
    if simplex_step:
      temp_result = minimize(self.residual, guess, args=(self.xvals, self.yvals, self.uvals), method="nelder")
      guess = copy.deepcopy(temp_result.params)

    # least squares estimate with error bars
    result = minimize(self.residual, guess, args=(self.xvals, self.yvals, self.uvals), method="leastsq")#, **{ 'xtol':1e-8, 'maxfev':500*(len(guess)+1) } )
    self.result = result
    return self.result








class FitSuperNode(object):
  '''
  This class stores
  - a list of FitNodes (or other FitSuperNodes!)
  - a global Parameters() object that applies to all nodes
  and then optimizes the global parameters using lmfit.minimize().
  During this process, the local parameters in each FitNode are also optimized.
  '''

  def __init__(self, nodelist=None, params=None):

    self.nodelist = nodelist
    self.params   = params
    self.result   = None


  def residual(self, params=None):

    if params == None:  params = self.params
    print params

    residual_list = []
    for node in self.nodelist:
      for key in params.keys():				# add global params to each local node
        node.params.add(key, value=params[key].value, vary=False)  
      result = node.fit()				# fit each node
      residual_list.append(np.array(result.residual))	# get local residual

    globalres = np.hstack(residual_list)
    return globalres



  def fit(self, simplex_step=True):

    # quick check for abuse of a SuperFit Node
    if self.params == None:
      for node in self.nodelist: node.fit()
      return None

    # make a copy of initial guess, to avoid changing it      
    guess = copy.deepcopy(self.params)

    # optional first estimate using simplex method
    if simplex_step:
      temp_result = minimize(self.residual, guess, method="nelder")
      guess = copy.deepcopy(temp_result.params)

    # final estimate using least squares for error bars
    result = minimize(self.residual, guess, method="leastsq")#, **{'xtol':1e-8, 'maxfev':500*(len(guess)+1)})
    self.result = result

    # update nodes with fits of global parameters, and return
    for node in self.nodelist:
      for key in self.result.params.keys():
        node.result.params[key].value = self.result.params[key].value
        node.result.params[key].stderr = self.result.params[key].stderr

    return self.result

  
    




