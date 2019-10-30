import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from math import sin, cos, radians
from plantenergy.utilities import sunflower_points

if __name__=='__main__':

    def sunflower_points(n, alpha=1.0):
        # this function generates n points within a circle in a sunflower seed pattern
        # the code is based on the example found at
        # https://stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle

        def radius(k, n, b):
            if (k + 1) > n - b:
                r = 1. # put on the boundary
            else:
                r = np.sqrt((k + 1.) - 1. / 2.) / np.sqrt(n - (b + 1.) / 2.)  # apply squareroot

            return r

        x = np.zeros(n)
        y = np.zeros(n)

        b = np.round(alpha * np.sqrt(n)) # number of boundary points

        phi = (np.sqrt(5.) + 1.) / 2.  # golden ratio

        for k in np.arange(0, n):

            r = radius(k, n, b)

            theta = 2. * np.pi * (k+1) / phi**2

            x[k] = r * np.cos(theta)
            y[k] = r * np.sin(theta)

        return x, y


    mpl.rc('font', family = 'serif', serif = 'cmr10')
    plt.rcParams["mathtext.fontset"] = "dejavuserif"

    fig = plt.figure(figsize=[4.5,2.5])
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))

    bladeX = np.array([3.,7.,10.,15.,20.,25.,30.,35.,30.,25.,20.,15.,10.,5.,3.,3.])
    bladeY = np.array([0.,0.,0.8,1.5,1.7,1.9,2.1,2.3,2.4,2.4,2.4,2.4,2.4,2.4,2.4,0.])-1.5

    shear = np.array([75,75,75,175,175,175,275,275,275])
    spacing = np.array([0.5,1.0,1.5,0.5,1.0,1.5,0.5,1.0,1.5])

    H = 110.
    D = 130.
    R = 65.
    d = np.array([6.3,5.5,4.])

    circle1 = plt.Circle((0.,H), R, color='C0', fill=False, linestyle = '--', linewidth=1.2*1.5,alpha=0.5)
    ax1.add_artist(circle1)
    circle1 = plt.Circle((0.,H), R, color='C0', fill=False, linestyle = '--', linewidth=1.2*1.5,alpha=0.5)
    ax2.add_artist(circle1)

    c1 = R/35.

    px1 = np.array([0.-d[0]/2,0.-d[1]/2,0.-d[2]/2,0.+d[2]/2,0.+d[1]/2,0.+d[0]/2,0.-d[0]/2])
    py1 = np.array([0,H/2,H-3.*c1,H-3.*c1,H/2,0,0])
    ax1.plot(px1,py1,color='C0', linewidth=1.2*1.5,alpha=0.5)
    ax2.plot(px1,py1,color='C0', linewidth=1.2*1.5,alpha=0.5)

    #add blades
    hub1 = plt.Circle((0.,H), 3*c1, color='C0', fill=False, linewidth=1*1.5,alpha=0.5)
    ax1.add_artist(hub1)
    hub1 = plt.Circle((0.,H), 3*c1, color='C0', fill=False, linewidth=1*1.5,alpha=0.5)
    ax2.add_artist(hub1)

    angle1 = 5.

    blade1X = bladeX*cos(radians(angle1))-bladeY*sin(radians(angle1))
    blade1Y = bladeX*sin(radians(angle1))+bladeY*cos(radians(angle1))

    blade2X = bladeX*cos(radians(angle1+120.))-bladeY*sin(radians(angle1+120.))
    blade2Y = bladeX*sin(radians(angle1+120.))+bladeY*cos(radians(angle1+120.))

    blade3X = bladeX*cos(radians(angle1+240.))-bladeY*sin(radians(angle1+240.))
    blade3Y = bladeX*sin(radians(angle1+240.))+bladeY*cos(radians(angle1+240.))

    ax1.plot(blade1X*c1+0., blade1Y*c1+H, linewidth=1*1.5, color='C0',alpha=0.5)
    ax1.plot(blade2X*c1+0., blade2Y*c1+H, linewidth=1*1.5, color='C0',alpha=0.5)
    ax1.plot(blade3X*c1+0., blade3Y*c1+H, linewidth=1*1.5, color='C0',alpha=0.5)

    angle1 = 65.

    blade1X = bladeX*cos(radians(angle1))-bladeY*sin(radians(angle1))
    blade1Y = bladeX*sin(radians(angle1))+bladeY*cos(radians(angle1))

    blade2X = bladeX*cos(radians(angle1+120.))-bladeY*sin(radians(angle1+120.))
    blade2Y = bladeX*sin(radians(angle1+120.))+bladeY*cos(radians(angle1+120.))

    blade3X = bladeX*cos(radians(angle1+240.))-bladeY*sin(radians(angle1+240.))
    blade3Y = bladeX*sin(radians(angle1+240.))+bladeY*cos(radians(angle1+240.))

    ax2.plot(blade1X*c1+0., blade1Y*c1+H, linewidth=1*1.5, color='C0',alpha=0.5)
    ax2.plot(blade2X*c1+0., blade2Y*c1+H, linewidth=1*1.5, color='C0',alpha=0.5)
    ax2.plot(blade3X*c1+0., blade3Y*c1+H, linewidth=1*1.5, color='C0',alpha=0.5)


    x_sample = np.array([0.69*R,0.,-0.69*R,0.])
    y_sample = np.array([0.,0.69*R,0.,-0.69*R])+H
    ax1.plot(x_sample,y_sample,'o',color='C2',alpha=0.5,markersize=3)

    ax1.text(0.69*R,H+7.,'(0.69,0)',horizontalalignment='center')
    ax1.text(0.,H+0.69*R+7.,'(0,0.69)',horizontalalignment='center')
    ax1.text(-0.69*R,H+7.,'(-0.69,0)',horizontalalignment='center')
    ax1.text(0.,H-0.69*R+7.,'(0,-0.69)',horizontalalignment='center')

    x_sample, y_sample = sunflower_points(100)
    x_sample = x_sample*R
    y_sample = y_sample*R + H
    ax2.plot(x_sample,y_sample,'o',color='C2',alpha=0.5,markersize=3)

    ax1.axis('equal')
    ax1.axis('off')

    ax2.axis('equal')
    ax2.axis('off')

    ax1.set_ylim(-15.,H+R+5)
    ax2.set_ylim(-15.,H+R+5)

    ax1.text(0.,-10.,'optimize',horizontalalignment='center')
    ax2.text(0.,-10.,'evaluate',horizontalalignment='center')

    ax1.text(-R+5.,H+R-10.,'a',fontsize=11,color='black')
    ax2.text(-R+5.,H+R-10.,'b',fontsize=11,color='black')

    plt.subplots_adjust(top = 1., bottom = 0., right = 1, left = 0.,
                hspace = 0, wspace = 0)

    # plt.savefig('rotor_samples.pdf',transparent=True)
    plt.savefig('/Users/ningrsrch/Dropbox/Projects/stanley2019-variable-reduction/paper/paper-figures/rotor_samples.pdf',transparent=True)

    plt.show()
