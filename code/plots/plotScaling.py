import numpy as np
import matplotlib.pyplot as plt

# ga_1200 = np.loadtxt("./data/ga_multistart_rundata_38turbs_nantucketWindRose_12dirs_BPA_all_tolgens1200.txt")
# ga_200 = np.loadtxt("./data/ga_multistart_rundata_38turbs_nantucketWindRose_12dirs_BPA_all_tolgens200.txt")
# snopt = np.loadtxt("./data/snopt_multistart_rundata_38turbs_nantucketWindRose_12dirs_BPA_all.txt")
#
# # ga_calls = ga_200[:,8]
# ga_calls = ga_1200[:,8]
# snopt_calls = snopt[:,8]
#
# fig = plt.figure(figsize=[5.,2.])
# ax1 = plt.subplot2grid((1, 2), (0, 0))
# ax2 = plt.subplot2grid((1, 2), (0, 1))
#
# ax1.spines['right'].set_visible(False)
# ax1.spines['top'].set_visible(False)
# ax1.yaxis.set_ticks_position('left')
# ax1.xaxis.set_ticks_position('bottom')
# ax2.spines['right'].set_visible(False)
# ax2.spines['top'].set_visible(False)
# ax2.yaxis.set_ticks_position('left')
# ax2.xaxis.set_ticks_position('bottom')
#
# ax1.tick_params(axis='both', which='major', labelsize=8)
# ax1.tick_params(axis='both', which='minor', labelsize=8)
# ax2.tick_params(axis='both', which='major', labelsize=8)
# ax2.tick_params(axis='both', which='minor', labelsize=8)
#
# ax1.set_yscale('log')
# bins = np.logspace(np.log10(1.),np.log10(1.E7),50)
#
# # ax1.hist(snopt_calls,bins=bins,color='C0',alpha=0.5)
# # ax1.hist(ga_calls,bins=bins,color='C3',alpha=0.5)
#
# ax1.boxplot(snopt_calls)
# ax1.boxplot(ga_calls)

import matplotlib as mpl
mpl.rc('font', family = 'serif', serif = 'cmr10')
fig = plt.figure(figsize=[3.,2.5])
ax2 = plt.subplot(111)

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')

ax2.tick_params(axis='both', which='major', labelsize=8)
ax2.tick_params(axis='both', which='minor', labelsize=8)

dimvec2 = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
snopt_grad = np.array([31, 35, 36, 45, 64, 68, 157, 183, 195, 244])
slsqp_grad = np.array([33, 30, 43, 76, 135, 271, 535, 271, 445, 577])
snopt_fd = np.array([102, 248, 375, 1017, 3054, 6315, 60142, 151453, 212337, 776704])
slsqp_fd = np.array([84, 118, 227, 621, 1646, 5324, 16877, 36639, 46939, 152124])
dimvec = [2, 4, 8, 16, 32, 64]
alpso = np.array([1150, 32780, 108040, 488240, 2649760, 12301760])


ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.plot(dimvec2,snopt_grad,'o',color='C0',alpha=0.5)
ax2.plot(dimvec2,snopt_fd,'o',color='C1',alpha=0.5)
ax2.plot(dimvec,alpso,'o',color='C3',alpha=0.5)



# ax1.set_xticks((1E1,1E3,1E5,1E7))
# ax1.set_xticklabels((r'10$^1$',r'10$^3$',r'10$^5$',r'10$^7$'))
# ax1.set_yticks((0,50,100,150))
# ax1.set_yticklabels(('0','50','100','150'))

ax2.set_xticks((1E1,1E2,1E3))
# ax2.set_xticklabels((r'10$^1$',r'10$^2$',r'10$^3$'))
ax2.set_yticks((1E1,1E3,1E5,1E7))
ax2.set_yticklabels((r'10$^1$',r'10$^3$',r'10$^5$',r'10$^7$'))

ax2.text(45, 20, 'analytic gradients',fontsize=8, color='C0')
ax2.text(45, 1*5e2, 'finite difference\ngradients',fontsize=8, color='C1')
ax2.text(30, 7e5, 'gradient-free',fontsize=8, color='C3')

ax2.set_ylabel('number of function\ncalls to optimize',fontsize=8)
ax2.set_xlabel('number of design variables',fontsize=8)

plt.minorticks_off()

plt.subplots_adjust(top = 0.94, bottom = 0.2, right = 0.98, left = 0.25,
            hspace = 0.6, wspace = 0.2)

plt.savefig('/Users/ningrsrch/Dropbox/Projects/reduction/paper/paper-figures/scaling.pdf',transparent=True)
plt.show()
