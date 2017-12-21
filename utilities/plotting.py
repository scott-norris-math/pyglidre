import numpy as np
import matplotlib as mpl ; mpl.use('Agg') ;
import matplotlib.pyplot as plt



def create_2D_intensity_frames(tqi, prefix=None):

  # unpack variables from the dictionary
  tvals  = tqi['tvals']
  qxvals = tqi['qxvals']
  qzvals = tqi['qzvals']
  iofqt  = tqi['iofqt']

  # set the filename
  if prefix==None:  prefix="iofq2D-"

  # create the frames
  for kk,(tt,ii) in enumerate(zip(tvals, iofqt)):

    plt.figure()
    plt.imshow(np.log10(np.transpose(ii)), origin="lower")
    plt.xlabel(r"qx [nm$^{-1}$]")
    plt.ylabel(r"qz [nm$^{-1}$]")
    plt.title(r"Heatmap of $\log_{10}$(intensity):  t=%8d sec" % (tt))
    plt.savefig("%s-heatmap-frame-%05d.png"%(prefix, kk))
    plt.close()
  
    plt.figure()
    plt.contour(qxvals, qzvals, np.log10(np.transpose(ii)), np.arange(0,8, 0.2) )
    plt.xlabel(r"qx [nm$^{-1}$]")
    plt.ylabel(r"qz [nm$^{-1}$]")
    plt.title(r"Countour Plot of $\log_{10}$(intensity):  t=%8d sec" % (tt))
    plt.colorbar()
    plt.savefig("%s-contour-frame-%05d.png"%(prefix, kk))
    plt.close()


def plot_iofqt_heatmap(tqi, filename=None):

  # unpack variables from dictionary
  tvals = tqi['tvals']
  qvals = tqi['qvals']
  iofqt = tqi['iofqt']

  # set the filename
  if filename==None:  filename="iofqt-heatmap.png"

  # Create a time/space heatmap of I(q,t).  
  #Works best if there are lots of time data points (continuous bombardment)
  fig = plt.figure()
  ax  = fig.add_subplot(111)
  ax.imshow(np.log10(tqi['iofqt']), origin="lower", extent=[qvals[0], qvals[-1], tvals[0], tvals[-1]], aspect='auto')
  #ax.set_aspect(1)
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
  if max_lines==None:  max_lines = 30
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
  #plt.ylim(-6, 0)
  plt.title(r"Selected values of $\log_{10}$(intensity)")
  plt.legend(prop={'size':8})
  plt.savefig(filename)
  plt.close()

