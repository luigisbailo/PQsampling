import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['style1'])

# f, ((ax1),(ax2)) = pl.subplots(2,1,sharex=True)

plt.figure(figsize=(3.37,2))
# f.set_figwidth(3.37)
dataBM = np.loadtxt('plotBM.out')
dataGF = np.loadtxt('plotGF_proj.out')


plt.plot (dataBM[:,0], dataBM[:,1],  '--', label = r'$BM$')
plt.plot (dataGF[:,0], dataGF[:,1], '.', label = r'$p_{\Omega}(r,t|0,0)$')

plt.xlabel (r'$MSD\thinspace [nm^2]$',fontsize=7)
plt.ylabel (r'$Time\thinspace [ns]$',fontsize=7)

lgd=plt.legend (fontsize=7,ncol=2,bbox_to_anchor=(0.5, 1.25), loc='upper center', columnspacing=0)

plt.tight_layout()
#pl.subplots_adjust(top=0.6)
art = []
art.append(lgd)
plt.show()
# pl.savefig ('Fig3.eps',dpi=600, additional_artist=art, bbox_inches="tight")
