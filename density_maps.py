import numpy as np
import sys
import os
import pandas as pd
from matplotlib.tri import Triangulation
randn = np.random.randn
from scipy.stats import norm
from scipy.stats import sem
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import plotly.express as px
import matplotlib as mpl
from scipy.interpolate import griddata
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import matplotlib.colors as mcolors

def colouring_stup(clr1,clr2):
  clr1_cmap = plt.cm.get_cmap(clr1)
  clr2_cmap = plt.cm.get_cmap(clr2)
  colors = np.vstack((clr1_cmap(np.linspace(0, 1, 256)), clr2_cmap(np.linspace(0, 1, 256))))
  custom_cmap = mcolors.LinearSegmentedColormap.from_list("mix", colors)
  return custom_cmap

def density_map_prep(lipid_density,tot_density,prot_density):
  load_lipid_data=np.loadtxt(lipid_density)
  load_tot_lipid_data=np.loadtxt(tot_density)
  load_prot_data=np.loadtxt(tot_denstiy)
  init_density=data_list[1:,1:]
  

  

