import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
sys.path.insert(0, '/Users/ningrsrch/Dropbox/Projects/reduction/')
from position_constraints import *
import scipy as sp
import scipy.special
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from windRoses import *

def weibull(x,mu):
    k = 2.0
    L = mu/(scipy.special.gamma(1.+1./k))
    return k/L*(x/L)**(k-1)*np.exp(-(x/L)**k)

"""colors"""
robin = '#55dfe6'
nectar = '#e3bf54'
tuscan = '#e37854'
royal = '#a054e3'

robin = 'C0'
nectar = 'C1'
tuscan = 'C1'
royal = 'C3'

import matplotlib as mpl
mpl.rc('font', family = 'serif', serif = 'cmr10')
plt.rcParams["mathtext.fontset"] = "dejavuserif"
# robin = '#ACBD78'

fig = plt.figure(figsize=[3.25,2.])
ax1 = plt.subplot2grid((1, 1), (0, 0))

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)


x = np.linspace(0.,25.,1000)

f = weibull(x,5.)
ax1.plot(x,f,color=tuscan,alpha=0.5,linewidth=2,label=r'$U_{\mathrm{mean}}=5$')

f = weibull(x,10.)
ax1.plot(x,f,color=tuscan,alpha=0.5,linestyle='--',linewidth=2,label=r'$U_{\mathrm{mean}}=10$')

ax1.set_xlim(0.,25.)
ax1.set_ylim(0.,0.2)

ax1.set_xlabel('wind speed',fontsize=8)
ax1.set_ylabel('probability',fontsize=8)

ax1.tick_params(axis='both', which='major', labelsize=8)
ax1.tick_params(axis='both', which='minor', labelsize=8)

ax1.legend(prop={'size': 8,'family':'serif'})

plt.subplots_adjust(top = 0.9, bottom = 0.22, right = 0.9, left = 0.22,
            hspace = 0.0, wspace = 0.)

plt.savefig('/Users/ningrsrch/Dropbox/Projects/stanley2019-variable-reduction/paper/paper-figures/weibullR1.pdf',transparent=True)
plt.show()
