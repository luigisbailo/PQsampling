import numpy as np
import matplotlib.pyplot as pl

pl.style.use(['style1'])

f, ((ax1),(ax2)) = pl.subplots(2,1,sharex=True)

f.set_figheight(3.)
f.set_figwidth(3.37)

data = np.loadtxt('plot.out')


# f, ((ax1),(ax2)) = pl.subplots(2,1,sharex=True)
ax1.plot (data[:,0], data[:,1], label = r'$t=0.01$',color = '#3366CC',linewidth=1.)
ax1.plot (data[:,0], data[:,2], '--' ,label = r'$t=0.1$',color = '#3366CC',linewidth=1.)
ax1.plot (data[:,0], data[:,3], label = r'$t=0.99$', color = '#CC6633',linewidth=1.)

ax2.plot (data[:,0], data[:,4],color = '#3366CC',linewidth=1.)
ax2.plot (data[:,0], data[:,5], '--', color = '#3366CC',linewidth=1.)
ax2.plot (data[:,0], data[:,6], color = '#CC6633',linewidth=1.)

ax2.set_xlabel (r'$r$',fontsize=10)

ax1.text(0.4,0.8,r'$g_{\Omega}(r,t|0,0)$',fontsize=10,transform=ax1.transAxes,bbox=dict(facecolor='white', alpha=0.5))
ax2.text (0.375,0.8,r'$g_{\Omega}(r,t|\tau;0,0)$',fontsize=10,transform=ax2.transAxes,bbox=dict(facecolor='white'))

ax1.legend (loc='upper right',fontsize=8)

ax1.tick_params(labelsize=10)
ax2.tick_params(labelsize=10)
ax1.set_ylim (0,5)
ax2.set_ylim (0,5)

pl.tight_layout()

f.subplots_adjust(wspace=0.05)
f.subplots_adjust(hspace=0.1)
pl.subplots_adjust(top=0.89)

pl.savefig ('Fig1.pdf',dpi=600)
