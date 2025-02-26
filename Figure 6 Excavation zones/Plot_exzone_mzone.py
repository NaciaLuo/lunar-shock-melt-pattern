import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as clr


def set_size(w, h, ax=None):
    if not ax: ax = plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)
    
    
plt.rcParams['font.size'] = '11'
plt.rcParams['font.stretch'] = 'condensed'

# Model info
modelname = 'I01_deg45_V15'
I = modelname[1:3]
deg = int(modelname[7:9]) # impact angle
d = float(I)*1000 # projectile diameter (km)
a = d/2 # projectile radius (km)
low_cppr = 15 # cell per projectile radius
cell_size = a/low_cppr

high_cppr = 60

workdir = '/home/xizi/Desktop/data repository/Figure 6 Excavation zones/'

# Ejecta (all)
savedir = workdir + modelname
ej_x00 = np.load(savedir + '/ej_x00.npy')
ej_y00 = np.load(savedir + '/ej_y00.npy')
ej_z00 = np.load(savedir + '/ej_z00.npy')
ej_v = np.load(savedir + '/ej_v.npy')

# Ejecta (escaped)
savedir_es = savedir + '/ejecta_escaped'
ej_v2 = np.load(savedir_es + '/ej_v.npy')
ej_x002 = np.load(savedir_es + '/ej_x00.npy')
ej_y002 = np.load(savedir_es + '/ej_y00.npy')
ej_z002 = np.load(savedir_es + '/ej_z00.npy')

vescape = 2.38e3
eslice = np.where(ej_y00 < 8*cell_size) 
eslice2 = (ej_y002 < 8*cell_size)

x00 = np.concatenate((ej_x00[eslice], ej_x002[eslice2]))
z00 = np.concatenate((ej_z00[eslice], ej_z002[eslice2]))
v00 = np.concatenate((ej_v[eslice], ej_v2[eslice2]))


proj_mfrac = np.load(savedir + '/Proj_mfrac.npy')
proj_xgrid = np.load(savedir + '/Proj_xgrid.npy')*1000
proj_zgrid = np.load(savedir + '/Proj_zgrid.npy')*1000


tar_mfrac = np.load(savedir + '/Tar_mfrac.npy')
tar_xgrid = np.load(savedir + '/Tar_xgrid.npy')*1000
tar_zgrid = np.load(savedir + '/Tar_zgrid.npy')*1000


# melt tracer 
melt_x = np.load(savedir + '/melt_x0.npy')*1000
melt_y = np.load(savedir + '/melt_y0.npy')*1000
melt_z = np.load(savedir + '/melt_z0.npy')*1000
melt_frac = np.load(savedir + '/melt_frac.npy')
plane = (melt_y < cell_size)

# Plot the excavation zone and the melt zone in the xz-plane


# ----- Create a colormap based on velocity -----
colors = ['darkblue', '#036eb7', '#2ca6e0', '#63c4e4', '#d1ecf4'][::-1]
cmap_name = 'custom_blue'
cm = clr.LinearSegmentedColormap.from_list(cmap_name, colors, N=300)

values_v = [50, 100, 200, 400, 1000, 2000, 3000]
norm_v = clr.BoundaryNorm(values_v, ncolors = 256)


pr_xx, pr_zz = np.meshgrid(proj_xgrid, proj_zgrid)
tar_xx, tar_zz = np.meshgrid(tar_xgrid, tar_zgrid)

pr_xx = (pr_xx + 0.5*a/high_cppr)[:-1, :-1]
pr_zz = (pr_zz + 0.5*a/high_cppr)[:-1, :-1]
tar_xx = (tar_xx + 0.5*a/high_cppr)[:-1, :-1]
tar_zz = (tar_zz + 0.5*a/high_cppr)[:-1, :-1]


scale = d
unit = r'$d_{proj}$'

# %%
fig, ax = plt.subplots()
ax.set_xlim(-10, 5)
ax.set_ylim(-2, 1.5)
ax.set_xlabel('x/{}'.format(unit))
ax.set_ylabel('z/{}'.format(unit))
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))

ax.tick_params(axis='both', which='both', direction='out')
ax.tick_params(axis='both', which='major', length=5)


# Melt zone contour
ax.contour(pr_xx/scale, pr_zz/scale, proj_mfrac[:,:,0], [0.], 
           colors='tab:orange', linestyles='--', linewidths=2., zorder=10)
ax.contour(tar_xx/scale, tar_zz/scale, tar_mfrac[:,:,0], [0.], 
           colors='tab:orange', linestyles='--', linewidths=2., zorder=10)

# Melt zone plotted using tracers, color = melt fraction (0-1)
ax.scatter(melt_x[plane]/scale, melt_z[plane]/scale, 
           c=melt_frac[plane], s=1, cmap='OrRd', edgecolor='None',
           vmin=0, vmax=1, 
           zorder=5)

# Excavation zone
pcm = ax.scatter(x00/scale, z00/scale, c=v00, 
           cmap = cm, norm = norm_v, 
           s=1)


# Surface & impactor
ax.axhline(0, c='k', lw=1) # Surface
impactor = plt.Circle((0,a/scale), a/scale, ls='-',
                      edgecolor='k', facecolor='None',
                      lw=0.75, zorder=0)
ax.add_patch(impactor)

# colorbar
cb = fig.colorbar(pcm, ax=ax, extend='both', shrink=0.3, 
                  orientation='horizontal', pad=0.2)
cb.set_label('v (m/s)')



set_size(10, 3)
ax.set_aspect('equal')
#fig.savefig(savedir + '/{}_ezone.png'.format(modelname), dpi=300, bbox_inches='tight')
