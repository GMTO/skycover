import matplotlib.pyplot as plt
import sys
import starbase
from numpy import *

if len(sys.argv) < 4:
    print "example:\npython plotgrid.py <probability file> <png output file> <title>"
    sys.exit(0)

tab = starbase.Starbase(sys.argv[1])

if (sys.argv[1].find("dgwf") >= 0):
    ttylim = 16.25
else:
    ttylim = 16.4


fig = plt.figure()
fig.suptitle('Sky Coverage at Galactic Pole', fontsize=14, fontweight='bold')

wfslim = 18.0

gmin = 13
wmin = 13
gmax= 20
wmax= 20
g = arange(gmin, gmax+1)
w = arange(wmin, wmax+1)
X, Y = meshgrid(w, g)

configs = ["dgnf"]
labels  = [sys.argv[3]]

ii={}
axlist=[]
zlist=[]
figlist=[]

levels={}
levels["dgnf"]           = array([0.2,0.4,0.6,0.8,0.9,0.95,0.99,1.0])
levels["folded"]         = array([0.2,0.4,0.6,0.8,0.9,0.95,0.99,1.0])
levels["folded-track60"] = array([0.2,0.4,0.6,0.8,0.9,0.95,0.99,1.0])
levels["phasing"]        = array([0.2,0.4,0.6,0.65,0.7,0.75,0.8])


for i,config in enumerate(configs):
    ii[config]=i
    fig=plt.figure()
    figlist.append(fig)
    ax = fig.add_subplot(1,1,1)
    axlist.append(ax)
    plt.grid(True)
    ax.set_ylabel('TT7 magnitude limit')
    ax.set_xlabel('WFS magnitude limit')
    ax.set_title(labels[i])
    zlist.append(zeros(X.shape))

# Populate data arrays with entries from table
#
for i in range(len(tab[:][1])):
    config = tab[:][1][i]
    wfsmag = int(tab[:][2][i])
    gdrmag = int(tab[:][3][i])
    prob   = float(tab[:][4][i])
    ax     = axlist[ii[config]]
    zlist[ii[config]][gdrmag-gmin][wfsmag-wmin] = prob

# Set any values that are equal to 1 to 1.01
for config in configs:
    a=zlist[ii[config]]
    a[a==1] = 1.01
    a[Y==20] = a[Y==19]
    a[X==20] = a[X==19]

for i,config in enumerate(configs):
    ax=axlist[ii[config]]
    fig=figlist[ii[config]]
    CS=ax.contour(X,Y,zlist[i],levels[config])
    ax.clabel(CS, inline=1, fontsize=10)

    ax.plot(wfslim, ttylim, "b+", markersize=20, markeredgewidth=2)

    fig.savefig(sys.argv[2])
