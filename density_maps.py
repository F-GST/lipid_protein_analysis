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

def lipid_density_map_prep(lipid_density,tot_density):
  load_lipid_data=np.loadtxt(lipid_density)
  load_tot_lipid_data=np.loadtxt(tot_density)
  init_lipid_density=load_lipid_data[1:,1:]
  init_tot_lipid_density=load_tot_lipid_data[1:,1:]
  normalized_lipid_dens= 100*(init_lipid_density  /init_tot_lipid_density)
  return normalized_lipid_dens
  
def prot_density_map_prep(prot_density):
  load_prot_data=np.loadtxt(prot_denstiy)
  prot_density=load_prot_data[1:,1:]
  x=load_prot_data[0,1:] 
  y=load_prot_data[1:,0]
  return x,y,prot_density

def prot_mask(lipid_density_list,prot_density,x,y):
  protein_data=np.zeros_like(prot_density)
  for j in range(len(x)):
    for k in range(len(y)):
      if prot_density[j,k]>0.0000000:
        protein_data[j,k]=1
      else:
        continue 
  for l in range(len(x)):
    for m in range(len(y)):
      for lipid in lipid_density_list:
        if lipid[l,m]> 20:
          protein_data[l,m]=0  
        else:
          continue
  mask = np.zeros_like(protein_data, dtype=bool)
  for n in range(len(x)):
    for o in range(len(y)):
      if protein_data[n,o]>0: 
        mask[n,o] = True 
  return mask
  
  def plt_masked_density(lipid_density,mask,x,y,lipid_name,custom_cmap,norm_max):
    masked_lipid_density=np.ma.masked_array(lipid_density,mask)
    plt.figure(figsize=(20,20))
    a=plt.imshow(masked_lipid_density,origin='lower',interpolation='nearest',cmap=custom_cmap,norm=Normalize(vmin=0, vmax=norm_max))#,extent=[xi.min(),xi.max(),yi.min(),yi.max()]) )
    cbar=plt.colorbar(a)
    cbar.ax.tick_params(labelsize=20)
    plt.xlabel('X [nm]',fontsize = 30,weight='bold')
    plt.ylabel('Y [nm]',fontsize = 30,weight='bold')
    plt.yticks(fontsize=20,weight='bold')
    plt.xticks(fontsize=20,weight='bold')
    plt.title('%s'%(lipid_name),fontsize=30,weight='bold')
    plt.tight_layout()
    plt.savefig('%s.png'%(lipid_name))
    return 
    
  
  

  

