import numpy as np
import matplotlib.pyplot as pl

pl.style.use(['style1'])

# f, ((ax1),(ax2)) = pl.subplots(2,1,sharex=True)

pl.figure(figsize=(3.37,2))
# f.set_figwidth(3.37)
dataBM = np.loadtxt('plotBM.out')
dataGF = np.loadtxt('plotGF.out')


pl.plot (dataBM[:,0], dataBM[:,1],  '--', label = r'$BM$')
pl.plot (dataGF[:,0], dataGF[:,1], 'v', label = r'$p_{\Omega}(r,t|0,0)/S_{\Omega}(t|0,0)$')
pl.plot (dataGF[:,0], dataGF[:,2], '^',label = r'$p_{\Omega}(r,t|0,0;\tau)$')

pl.xlabel (r'$MSD\thinspace [nm^2]$',fontsize=7)
pl.ylabel (r'$Time\thinspace [ns]$',fontsize=7)

lgd=pl.legend (fontsize=7,ncol=2,bbox_to_anchor=(0.5, 1.4), loc='upper center', columnspacing=0)

pl.tight_layout()
#pl.subplots_adjust(top=0.6)
art = []
art.append(lgd)
pl.savefig ('Fig2.eps',dpi=600, additional_artist=art, bbox_inches="tight")
