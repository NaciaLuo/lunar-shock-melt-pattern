# Import required libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.colors as clr
from math import pi, sin

plt.rcParams['font.size'] = '11'
plt.rcParams['font.stretch'] = 'condensed'
    
# ----- Parameters -----
g=1.62; Ve=2.38e3 #Moon gravity and escape velocity


def grid(boundx, boundy, grid_cell):
    '''
    Grid for calculating the thickness
    '''
    gridx_l = np.arange(0, boundx[0]-grid_cell, -grid_cell)[::-1] 
    gridx_r = np.arange(grid_cell, boundx[1]+grid_cell, grid_cell)

    gridx = np.hstack((gridx_l, gridx_r))
    gridy = np.arange(0, boundy[1]+grid_cell, grid_cell)
    plotx, ploty = np.meshgrid(gridx, gridy)
    
    return gridx, gridy, plotx, ploty


def Melt_pattern(landx, landy, mf, boundx, boundy, grid_cell):
    '''
    Fraction of melt in ejecta
    Using equal-area square bins around the impact point
    to calculate the melt pattern
    
    mf = sum(mf_i)/n 
    '''
    if len(landx) != len(landy):
        print("Error: x and y do not have the same length. Stop.")
        return 

    gridx, gridy, plotx, ploty = grid(boundx, boundy, grid_cell)
    
    pattern = np.zeros([len(gridy)-1,len(gridx)-1])
    N = np.zeros([len(gridy)-1,len(gridx)-1])
    std = np.zeros([len(gridy)-1,len(gridx)-1])
    
    for i, xi in enumerate(gridx[:-1]):
        for j, yj in enumerate(gridy[:-1]):
            bound = ((landx>xi) & (landx<=gridx[i+1]) & \
                             (landy>yj) & (landy<=gridy[j+1]))

            if len(mf[bound].flatten()) < 1: #less than 1 big tracer
                pattern[j, i] = -1.
            else:
                pattern[j, i] = np.mean(mf[bound])
                N[j, i] = len(mf[bound])
                std[j, i] = np.std(mf[bound])
    return pattern, N, std, plotx, ploty


def _interpolate(x, x0, x1, y0, y1):

    if x1 > x0:
        return np.interp(0., [x0, x1], [y0, y1])
    else:
        return np.interp(0., [x1, x0], [y1, y0])
    
    
def Rtc_scaling(L,Vi,g,theta):
    '''
    Transient crater radius from scaling law
    impact onto solid rock; gravity dominant
    assuming same material of target and projectile
    Reference: Earth Impact Effects Program (Collins+2005, Eqn. 21)
    '''
    R_Tc = 0.5*1.161*(L**0.78)*(Vi**0.44)*(g**(-0.22))*sin(theta*pi/180)**(1./3)
    return R_Tc



workdir = '/home/xizi/Desktop/data repository/Figure 9 Melt and ejecta patterns/'

deg = 90
I = 30 # km
U = 15 # km/s
modelname = 'I{:02d}_deg{}_V{}'.format(I, deg, U)

ejmodel = modelname.format(deg)
melt_model = 'CPPR60_melt model'
print "Ejecta model:", ejmodel
print "Melt model:", melt_model

savedir = workdir + ejmodel + '/'

landx = np.load(savedir + 'ej_landx.npy')/1000.
landy = np.load(savedir + 'ej_landy.npy')/1000.
landmf = np.load(savedir + 'ej_mf.npy')

d = float(modelname[1:3])
dx = d/30 # low-res model: CPPR 15
truvol = dx**3 # tracer volume, m^3

# ----- Transient crater -----
crater_data = np.loadtxt(savedir + 'crater_metric.txt', skiprows=1)
Vtc = crater_data[:,1]
Dcr = crater_data[:,-2][Vtc == Vtc.max()][0]
Ddr = crater_data[:,-1][Vtc == Vtc.max()][0]
Rtc = (Dcr+Ddr)*0.25

# Transient crater topography
topox, topoy, topoz = np.load(savedir + 'Integrated_surface.npy')

# Get transient crater outline [x_0, y_0]
x_0 = [] 
y_0 = []
jjmax = 0
for ii in range(topox.shape[0]):
    for jj in range(topox.shape[1]-1):
        if topoz[ii, jj] < 0. and topoz[ii, jj+1] >= 0.:
            y_0.append(_interpolate(
                0., topoz[ii, jj], topoz[ii, jj+1],
                topoy[0, jj], topoy[0, jj+1]
                ))
            x_0.append(topox[ii, 0])
            break 
x_0 = np.array(x_0)
y_0 = np.array(y_0)

x_c = (x_0.min() + x_0.max())/2


Rtc_scaled = Rtc_scaling(I*1000., U*1000.,g, deg)/1000
print "iSALE-derived transient crater radius:", Rtc, 'km'
print "Scaling-law-derived transient crater radius:", round(Rtc_scaled,1), 'km'

landx = landx - x_c # shift by the crater center offset from impact site


# %%
l = 2

grid_cell = d*l
boundx = [min(landx), max(landx)]
boundy = [0, max(landy)]


pattern, N, std, plotx, ploty = Melt_pattern(landx, landy, landmf, boundx, boundy, grid_cell)

# Mask out grids where there are no ejecta
pattern_ma = np.ma.array(pattern, mask=(pattern<0))
N_ma = np.ma.array(N, mask=(pattern<0))
std_ma = np.ma.array(std, mask=(pattern<0))


# Ejecta & melt blanket thicknesses obtained from 'pattern' and 'N'
grid_A = grid_cell**2 # cell area

# Ejecta thickness (containing both melted & unmelted material)
TE = N_ma*truvol/grid_A*1000 # meter

# Melt thickness
ME = TE*pattern_ma
ME_mask = np.ma.array(ME, mask=(ME==0))

# Melt fraction
MF = ME_mask/TE*100


plotxc = plotx[:-1,:-1] + grid_cell/2
plotyc = ploty[:-1,:-1] + grid_cell/2

# Azimuthal angle w.r.t. the origin (in this case the impact point)
gridphi = np.arctan2(plotyc, plotxc)*180/pi
    
# Radial distance w.r.t. the origin
gridr = np.sqrt(plotxc**2 + plotyc**2)



# %%  Plotting: Ejecta thickness

## Values for plotting discrete thickness colormap
levels0 = [0,0.3,1,3,10,30,100]
#levels0 = [0,3,10,30,100,300,1000]



colors = ['#4576b6', '#91bfdc', '#e1f3fa', '#fce091', '#f68d5c', '#d93127']  
norm = clr.BoundaryNorm(levels0, ncolors =260)


xlim = [-40,20]
ylim = [-35,35]

         
fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(111, aspect=1)
ax.set_xlim(xlim)
ax.set_ylim(ylim)

# Ejecta thickness
im=ax.contourf(plotxc/Rtc, (plotyc-grid_cell/2)/Rtc, TE, levels=levels0, 
            colors=colors,extend='max')
ax.contourf(plotxc/Rtc, -(plotyc-grid_cell/2)/Rtc, TE, levels=levels0, 
            colors=colors,extend='max')

# Transient crater
ax.plot((x_0-x_c)/Rtc, y_0/Rtc, '-', c='gray', lw=2)
ax.plot((x_0-x_c)/Rtc, -y_0/Rtc, '-', c='gray', lw=2)

cb=fig.colorbar(im, ax=ax, shrink=0.8, extend='both',ticks=levels0,
               orientation='horizontal')
cb.set_label('Ejecta thickness (m)', labelpad=10)

ax.set_xlabel(r'$x/R_{tc}$, $R_{tc}$'+' = {:.0f} km'.format(Rtc))
ax.set_ylabel(r'$y/R_{tc}$')

ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
ax.tick_params(axis='both', direction='in', length=5)

plt.tight_layout()

fig.savefig(savedir + 'ejecta_thickness.png', dpi=300, bbox_inches='tight')



# %%  Plotting: Melt thickness


## Values for plotting discrete thickness colormap
levels0 = [3e-3, 1e-2, 0.03,1e-1, 0.3, 1, 3] 
#levels0 = [3e-2,1e-1,3e-1,1,3,10,30] 


xlim = [-40,20]
ylim = [-35,35]


fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(111, aspect=1)
ax.set_xlim(xlim)
ax.set_ylim(ylim)



norm = clr.BoundaryNorm(levels0, ncolors =260)

# Melt thickness
im = ax.contourf(plotxc/Rtc, (plotyc-grid_cell/2)/Rtc, ME_mask, levels=levels0,
              cmap='coolwarm', norm=norm, extend='both')

ax.contourf(plotxc/Rtc, -(plotyc-grid_cell/2)/Rtc, ME_mask, levels=levels0,
              cmap='coolwarm', norm=norm, extend='both')


# Transient crater
ax.plot((x_0-x_c)/Rtc, y_0/Rtc, '-', c='gray', lw=2)
ax.plot((x_0-x_c)/Rtc, -y_0/Rtc, '-', c='gray', lw=2)


ax.set_xlabel(r'$x/R_{tc}$, $R_{tc}$'+' = {:.0f} km'.format(Rtc))
ax.set_ylabel(r'$y/R_{tc}$')

ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
ax.tick_params(axis='both', direction='in', length=5)

cb=fig.colorbar(im, ax=ax, shrink=0.8, extend='both',
               orientation='horizontal')
cb.set_label('Melt thickness (m)', labelpad=10)

plt.tight_layout()

fig.savefig(savedir + 'melt_thickness.png', dpi=300, bbox_inches='tight')
         

# %%  Plotting: Melt% in ejecta

## Values for plotting discrete melt percent colormap
mf_levels = [0,1,2,5,10,20,50]

xlim = [-40,20]
ylim = [-35,35]


fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(111, aspect=1)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
    

norm = clr.BoundaryNorm(mf_levels, ncolors = 256)
cmap = cm.coolwarm


im = ax.contourf(plotxc/Rtc, (plotyc-grid_cell/2)/Rtc, MF, levels=mf_levels,
            cmap=cmap, norm=norm, extend='max')
ax.contourf(plotxc/Rtc, -(plotyc-grid_cell/2)/Rtc, MF, 
            levels=mf_levels,
            cmap=cmap, norm=norm, extend='max')


# Transient crater
ax.plot((x_0-x_c)/Rtc, y_0/Rtc, '-', c='gray', lw=2)
ax.plot((x_0-x_c)/Rtc, -y_0/Rtc, '-', c='gray', lw=2)


cb=fig.colorbar(im, ax=ax, shrink=0.8, extend='max',
               orientation='horizontal')
cb.set_label('Melt in ejecta (vol%)', labelpad=10)

ax.set_xlabel(r'$x/R_{tc}$, $R_{tc}$'+' = {:.0f} km'.format(Rtc))
ax.set_ylabel(r'$y/R_{tc}$')


ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
ax.tick_params(axis='both', direction='in', length=5)

plt.tight_layout()

fig.savefig(savedir + 'melt_percent.png', dpi=300, bbox_inches='tight')
