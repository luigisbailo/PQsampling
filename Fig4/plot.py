import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['style1'])

f, (ax1,ax2) = plt.subplots(2,1,sharex=True)

f.set_figheight(3.37)
f.set_figwidth(3.37)
data_diff = np.loadtxt('plot_diff_proj.out')
data_annih = np.loadtxt('plot_annih_proj.out')

marker_size=6
x = np.arange(0,1000,1)
ax1.plot (x, 6*0.01*x,  '--', color='grey', linewidth=1.)
ax1.plot (data_diff[:,0], data_diff[:,1], 'v', color='#336699', markersize=marker_size, label = r'$g_{\Omega}(r,t|0,0)$')
ax1.plot (data_diff[:,0], data_diff[:,2], '^', color='#CC6633', markersize=marker_size, label = r'$g_{\Omega}(r,t|0,0;\tau)$')
ax1.plot (data_diff[:,0], data_diff[:,3],  'x', color='#446280', markersize=4, label = r'$BD$')

x_arr = np.array([3,7,11,15,19])
ax2.plot (data_annih[x_arr,0], data_annih[x_arr,1], 'v', color='#336699',markersize=marker_size, label = r'$g_{\Omega}(r,t|0,0)$')
ax2.plot (data_annih[x_arr,0], data_annih[x_arr,2], '^', color='#CC6633', markersize=marker_size,label = r'$g_{\Omega}(r,t|0,0;\tau)$')
ax2.plot (data_annih[x_arr,0], data_annih[x_arr,3],  'x', color='#446280', markersize=4, label = r'$BD$')

ax1.set_ylabel (r'$MSD\thinspace [nm^2]$',fontsize=10)
ax2.set_ylabel (r'$N$',fontsize=10)
ax2.set_xlabel (r'$Time\thinspace [ns]$',fontsize=10)
ax2.set_xticks ([0,250,500,750,1000])
ax1.set_yticks ([0,20,40,60])
ax2.set_yticks ([0,2,4,6,8])

lgd=ax1.legend (fontsize=8,ncol=1, loc='upper left', columnspacing=0)

art = []
art.append(lgd)
plt.savefig ('Fig4.pdf',dpi=600, additional_artist=art, bbox_inches="tight")
