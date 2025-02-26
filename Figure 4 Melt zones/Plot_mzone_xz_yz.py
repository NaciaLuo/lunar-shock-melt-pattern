# import required libraries
import numpy as np
from os import path, mkdir
import matplotlib.pyplot as plt
from matplotlib import ticker
from math import pi, cos, sin
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Set default font size and style
plt.rcParams['font.size'] = '11'
plt.rcParams['font.stretch'] = 'condensed'


# Define a function to plot the impact direction based on impact angle
def arrow(ax, angle):
    r = 0.45
    x = cos(angle*pi/180)*r
    y = sin(angle*pi/180)*r

    arroww = mpatches.FancyArrowPatch((0.0,0.5),(0.0-x, 0.5-y),
                                     arrowstyle='Simple', mutation_scale=22,
                                     fc='k', ec='w', lw=1)
    ax.add_patch(arroww)
    

# Model info
modelname = 'I06_deg90_V15'
I = modelname[1:3] 
deg = int(modelname[7:9]) # impact angle
d = float(I) # projectile diameter (km)
a = d/2 # projectile radius (km)
cppr = 60 # cell per projectile radius

workdir = '/home/xizi/Desktop/data repository/Figure 4 Melt zones/'
savedir = workdir + modelname
if not path.exists(savedir):
    mkdir(savedir)

# projectile tracer grid
proj_xgrid = np.load(savedir + '/Proj_xgrid.npy')
proj_ygrid = np.load(savedir + '/Proj_ygrid.npy')
proj_zgrid = np.load(savedir + '/Proj_zgrid.npy')
proj_trp = np.load(savedir + '/Proj_TrP.npy') # peak shock pressure
proj_melt = np.load(savedir + '/Proj_mfrac.npy') # melt fraction (0-1)

# target tracer grid
tar_xgrid = np.load(savedir + '/Tar_xgrid.npy')
tar_ygrid = np.load(savedir + '/Tar_ygrid.npy')
tar_zgrid = np.load(savedir + '/Tar_zgrid.npy')
tar_trp = np.load(savedir + '/Tar_TrP.npy') # peak shock pressure
tar_melt = np.load(savedir + '/Tar_mfrac.npy') # melt fraction (0-1)


# Create figure
xlim1 = [-2.2, 0.5]
xlim2 = [-1., 1]
width_ratio = (xlim1[1]-xlim1[0])/(xlim2[1]-xlim2[0])
scale = d
h = 2.8
w = h*1.8
fig = plt.figure(figsize=(w,h))
gs = gridspec.GridSpec(nrows=1, ncols=3, width_ratios=[width_ratio,1,0])
cmap = matplotlib.cm.RdYlBu_r


# xz-plane data
data0 = np.ma.array(proj_trp, mask=(proj_trp<0.01))[:,:,0] # Get rid of small values
data1 = np.ma.array(tar_trp, mask=(tar_trp<0.01))[:,:,0] 

# cross-section at y = 0
data3 = proj_melt[:,:,0] 
data4 = tar_melt[:,:,0] 

# create meshgrid for plotting
xx0, yy0  = np.meshgrid(proj_xgrid, proj_zgrid)
xx1, yy1 = np.meshgrid(tar_xgrid, tar_zgrid)

# make the grid cell-centered 
xx3 = (xx0 + 0.5*a/cppr)[:-1, :-1]
yy3 = (yy0 + 0.5*a/cppr)[:-1, :-1]
xx4 = (xx1 + 0.5*a/cppr)[:-1, :-1]
yy4 = (yy1 + 0.5*a/cppr)[:-1, :-1]



ax1 = fig.add_subplot(gs[0, 0], aspect='equal')
ax2 = fig.add_subplot(gs[0, 1], aspect='equal')
ax3 = fig.add_subplot(gs[0, 2]) # for placing the colorbar

ax1.set_xlim(xlim1[0], xlim1[1])
ax2.set_xlim(xlim2[0], xlim2[1])

ax1.set_ylim(-1.5, 1.1)
ax2.set_ylim(-1.5, 1.1)



# Plot the pressure
pcm = ax1.pcolormesh(xx0/scale, yy0/scale , data0, cmap=cmap,
                     vmin=0, vmax=200)
ax1.pcolormesh(xx1/scale, yy1/scale, data1, cmap=cmap,
               vmin=0, vmax=200)


# Plot the melt zone contour
ax1.contour(xx3/scale, yy3/scale, data3, extend='both', 
             levels=[0], colors='gray', linewidths=2,linestyles='--')
ax1.contour(xx4/scale, yy4/scale, data4, extend='both', 
             levels=[0], colors='gray', linewidths=2,linestyles='--')


ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

ax1.set_xlabel(r'x/$d_{proj}$')
ax1.set_ylabel(r'z/$d_{proj}$')

ax1.tick_params(axis='both', which='both', direction='in', right=True, top=True)
ax1.tick_params(axis='both', which='major', length=5)

ax1.text(-2, 0.6, r'{}$\degree$'.format(deg),fontsize=16,fontweight='bold')



# Get where the melt zone width is the widest (yz plane)
y_widths = []
K, I, J = tar_melt.shape
for i in range(I):
    y_width = np.where(tar_melt[:,i,:]!=0)[1]
    if len(y_width) == 0:
        y_widths.append(0)
    else:
        y_width = np.where(tar_melt[:,i,:]!=0)[1].max()
        y_widths.append(y_width)

vline = tar_xgrid[np.where(y_widths==max(y_widths))].mean()/scale
#vline = 0
ax1.axvline(vline, ls='--', color='k', lw=1)

x_i = np.argmin(abs(tar_xgrid/scale-vline))
data = tar_trp[:,x_i,:] # the cross-section at x = vline
xx, yy = np.meshgrid(tar_ygrid, tar_zgrid)

data5 = tar_melt[:,x_i,:] 
xx5 = (xx + 0.5*a/cppr)[:-1, :-1]
yy5 = (yy + 0.5*a/cppr)[:-1, :-1]


# Plot the pressure 
pcm = ax2.pcolormesh(xx/scale, yy/scale , data, cmap=cmap,
                     vmin=0, vmax=200)
ax2.pcolormesh(-xx/scale, yy/scale, data, cmap=cmap,
               vmin=0, vmax=200)


# Plot the melt zone contour
ax2.contour(xx5/scale, yy5/scale, data5, extend='both', 
                  levels=[0], colors='gray', linewidths=2, linestyles='--')
ax2.contour(-xx5/scale, yy5/scale, data5, extend='both', 
             levels=[0], colors='gray', linewidths=2, linestyles='--')


ax2.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax2.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))


ax2.set_xlabel(r'$x/d_{proj}$')
ax2.set_yticklabels('')
ax2.tick_params(axis='both', which='both', direction='in', right=True, top=True)
ax2.tick_params(axis='both', which='major', length=5)

# Plot the surface horizon
ax1.axhline(0, c='k', lw=2, zorder=10)
ax2.axhline(0, c='k', lw=2, zorder=10)
# Plot the impactor
impactor1 = plt.Circle((0,a/scale), a/scale, ls='-',
                    color='k',lw=1.5,fill=False, zorder=10)
impactor2 = plt.Circle((0,a/scale), a/scale, ls='--',
                    color='k',lw=1.5,fill=False, zorder=10)
ax1.add_patch(impactor1)
ax2.add_patch(impactor2) 

# Plot the impact direction
arrow(ax1, deg)


# Plot colorbar
ax3.axes.set_axis_off()
divider = make_axes_locatable(ax3)
cx = divider.append_axes("right", size="15%", pad=0.7)
cb = fig.colorbar(pcm, cax=cx, ticks=[0,50,100,150,200])
cb.set_label('Peak shock pressure (GPa)')


plt.subplots_adjust(wspace=0.05)


#fig.savefig(savedir + '/{}_mzone.png'.format(modelname),dpi=300, bbox_inches='tight')

