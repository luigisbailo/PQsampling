import numpy as np
import matplotlib.pyplot as pl

pl.style.use(['style1'])

f, ((ax1),(ax2)) = pl.subplots(2,1,sharex=True)

f.set_figheight(4)
f.set_figwidth(3.37)

data = np.loadtxt('plot.out')


# f, ((ax1),(ax2)) = pl.subplots(2,1,sharex=True)
ax1.plot (data[:,0], data[:,1], label = r'$t=0.01$')
ax1.plot (data[:,0], data[:,2], label = r'$t=0.1$')
ax1.plot (data[:,0], data[:,3], label = r'$t=0.99$')

ax2.plot (data[:,0], data[:,4])
ax2.plot (data[:,0], data[:,5])
ax2.plot (data[:,0], data[:,6])

ax2.set_xlabel (r'$r$')

ax1.set_title (r'$p_{\Omega}(r,t|0,0)/S_{\Omega}(t|0,0)$')
ax2.set_title (r'$p_{\Omega}(r,t|\tau;0,0)$')

ax1.legend (loc='upper right')

ax1.set_ylim (0,6)
ax2.set_ylim (0,6)

pl.tight_layout()

# f.subplots_adjust(wspace=0.05)
# f.subplots_adjust(hspace=0.1)
# pl.subplots_adjust(top=0.89)

pl.savefig ('Fig1.eps',dpi=600)
