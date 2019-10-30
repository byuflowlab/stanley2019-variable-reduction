import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
sys.path.insert(0, '/Users/ningrsrch/Dropbox/Projects/reduction/')
from position_constraints import *
import scipy as sp
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

"""colors"""
robin = '#B2DBD5'
nectar = '#EAB364'
tuscan = '#B2473E'
olive = '#ACBD78'

robin = 'C0'
nectar = 'C1'
tuscan = 'C2'
royal = 'C3'

# dirs = np.array([5,10,15,20,25,30,50])
# speeds = np.ones(7,dtype=int)*10
#
#
# # dirs = np.ones(7,dtype=int)*30
# # speeds = np.array([1,2,3,4,5,6,10])
#
#
# nDirs = dirs[0]
# nSpeeds = speeds[0]
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_opt.txt'%(nDirs,nSpeeds))
# opt5 = data[:]/358000.
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_true.txt'%(nDirs,nSpeeds))
# true5 = (data[:]*0.936)/358000.
#
# nDirs = dirs[1]
# nSpeeds = speeds[1]
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_opt.txt'%(nDirs,nSpeeds))
# opt10 = data[:]/358000.
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_true.txt'%(nDirs,nSpeeds))
# true10 = (data[:]*0.936)/358000.
#
# nDirs = dirs[2]
# nSpeeds = speeds[2]
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_opt.txt'%(nDirs,nSpeeds))
# opt15 = data[:]/358000.
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_true.txt'%(nDirs,nSpeeds))
# true15 = (data[:]*0.936)/358000.
#
# nDirs = dirs[3]
# nSpeeds = speeds[3]
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_opt.txt'%(nDirs,nSpeeds))
# opt20 = data[:]/358000.
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_true.txt'%(nDirs,nSpeeds))
# true20 = (data[:]*0.936)/358000.
#
# nDirs = dirs[4]
# nSpeeds = speeds[4]
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_opt.txt'%(nDirs,nSpeeds))
# opt25 = data[:]/358000.
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_true.txt'%(nDirs,nSpeeds))
# true25 = (data[:]*0.936)/358000.
#
# nDirs = dirs[5]
# nSpeeds = speeds[5]
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_opt.txt'%(nDirs,nSpeeds))
# opt30 = data[:]/358000.
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_true.txt'%(nDirs,nSpeeds))
# true30 = (data[:]*0.936)/358000.
#
# nDirs = dirs[6]
# nSpeeds = speeds[6]
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_opt.txt'%(nDirs,nSpeeds))
# opt50 = data[:]/358000.
# data = np.loadtxt('data/bins-convergence/%s_%s/AEP_true.txt'%(nDirs,nSpeeds))
# true50 = (data[:]*0.936)/358000.
#
# arr = np.zeros((7,7))
# for i in range(7):
#     arr[i][j] =
# for i in range(1):
#     # plt.plot(1,true5[i+k],'o',color='C%s'%i)
#     # plt.plot(2,true10[i+k],'o',color='C%s'%i)
#     # plt.plot(3,true15[i+k],'o',color='C%s'%i)
#     # plt.plot(4,true20[i+k],'o',color='C%s'%i)
#     # plt.plot(5,true25[i+k],'o',color='C%s'%i)
#     # plt.plot(6,true30[i+k],'o',color='C%s'%i)
#     # plt.plot(10,true50[i+k],'o',color='C%s'%i)
#     #
#     # plt.plot(1,opt5[i+k],'s',color='C%s'%(i+1))
#     # plt.plot(2,opt10[i+k],'s',color='C%s'%(i+1))
#     # plt.plot(3,opt15[i+k],'s',color='C%s'%(i+1))
#     # plt.plot(4,opt20[i+k],'s',color='C%s'%(i+1))
#     # plt.plot(5,opt25[i+k],'s',color='C%s'%(i+1))
#     # plt.plot(6,opt30[i+k],'s',color='C%s'%(i+1))
#     # plt.plot(10,opt50[i+k],'s',color='C%s'%(i+1))
#
#     plt.plot(1,(true5[i+k]-opt5[i+k])/true5[i+k],'s',color='C%s'%(i+1))
#     plt.plot(2,(true10[i+k]-opt10[i+k])/true10[i+k],'s',color='C%s'%(i+1))
#     plt.plot(3,(true15[i+k]-opt15[i+k])/true15[i+k],'s',color='C%s'%(i+1))
#     plt.plot(4,(true20[i+k]-opt20[i+k])/true20[i+k],'s',color='C%s'%(i+1))
#     plt.plot(5,(true25[i+k]-opt25[i+k])/true25[i+k],'s',color='C%s'%(i+1))
#     plt.plot(6,(true30[i+k]-opt30[i+k])/true30[i+k],'s',color='C%s'%(i+1))
#     plt.plot(10,(true50[i+k]-opt50[i+k])/true50[i+k],'s',color='C%s'%(i+1))
#
# # for i in range(10):
# #     plt.plot(1,true5[i],'o',color='C%s'%i)
# #     plt.plot(2,true10[i],'o',color='C%s'%i)
# #     plt.plot(3,true15[i],'o',color='C%s'%i)
# #     plt.plot(4,true20[i],'o',color='C%s'%i)
# #     plt.plot(5,true25[i],'o',color='C%s'%i)
# #     plt.plot(6,true30[i],'o',color='C%s'%i)
# #     plt.plot(10,true50[i],'o',color='C%s'%i)
#
# # for i in range(10):
# #     plt.plot(5,true5[i],'o',color='C%s'%i)
# #     plt.plot(10,true10[i],'o',color='C%s'%i)
# #     plt.plot(15,true15[i],'o',color='C%s'%i)
# #     plt.plot(20,true20[i],'o',color='C%s'%i)
# #     plt.plot(25,true25[i],'o',color='C%s'%i)
# #     plt.plot(30,true30[i],'o',color='C%s'%i)
# #     plt.plot(50,true50[i],'o',color='C%s'%i)
#
#
# # ax1.text(-2500.,2700.,'c',fontsize=10,color='black')
# # ax2.text(-2500.,2700.,'b',fontsize=10,color='black')
# # ax3.text(-2500.,2700.,'a',fontsize=10,color='black')
# # ax5.text(-2500.,2700.,'d',fontsize=10,color='black')
# # ax5.text(750.,4500.,'e',fontsize=10,color='black')
# # ax5.text(-500.,-4000.,'f',fontsize=10,color='black')
#
#
# # plt.ylim(0.98,1.005)
# # plt.tight_layout()
# # plt.savefig('/Users/ningrsrch/Dropbox/Projects/reduction/paper/paper-figures/layout-variable.pdf')
# plt.show()


dirs = np.array([5,10,15,20,25,30,50])
speeds = np.array([1,2,3,4,5,6,10])

arr_true = np.zeros((10,7,7))
arr_diff = np.zeros((10,7,7))

for i in range(7):
    for j in range(7):
        nDirs = dirs[i]
        nSpeeds = speeds[j]
        data = np.loadtxt('data/bins-convergence/%s_%s/AEP_opt.txt'%(nDirs,nSpeeds))
        opt = data[:]
        data = np.loadtxt('data/bins-convergence/%s_%s/AEP_true.txt'%(nDirs,nSpeeds))
        tr = (data[:]*0.936)

        for k in range(10):
            arr_true[k][i][j] = tr[k]
            arr_diff[k][i][j] = (opt[k]-tr[k])/tr[k]*100.

dirs = np.array([5,10,15,20,25,30,55])
speeds = np.array([1,2,3,4,5,6,11])

DIRS, SPEEDS = np.meshgrid(dirs,speeds)

for i in range(10):
    fig = plt.figure(i)
    # im = plt.gca().imshow(arr_true[i]/358000.,cmap='Reds')
    im = plt.gca().pcolormesh(DIRS-2.5, SPEEDS-0.5, arr_true[i]/358000.,cmap='Reds')
    # im = plt.gca().pcolormesh(DIRS, SPEEDS, arr_diff[i],vmin=-2.,vmax=2.,cmap='Reds',shading='gourand')
    cbar = plt.colorbar(im)
    # cbar.ax.set_yticklabels(['0','1','2','>3'])
    cbar.set_label('normalized "true" AEP', rotation=270,labelpad=15)

    plt.xlabel('# directions')
    plt.ylabel('# speeds')
    plt.title('starting layout %s'%(i+1))

plt.show()
