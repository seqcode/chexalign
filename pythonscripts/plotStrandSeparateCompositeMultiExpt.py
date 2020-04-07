import sys
from numpy import *
import numpy as np
from sklearn import metrics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage

def clamp(val, minimum=0, maximum=255):
    if val < minimum:
        return minimum
    if val > maximum:
        return maximum
    return val

def colorscale(hexstr, scalefactor):
    """
    Scales a hex string by ``scalefactor``. Returns scaled hex string.

    To darken the color, use a float value between 0 and 1.
    To brighten the color, use a float value greater than 1.

    >>> colorscale("#DF3C3C", .5)
    #6F1E1E
    >>> colorscale("#52D24F", 1.6)
    #83FF7E
    >>> colorscale("#4F75D2", 1)
    #4F75D2
    """

    hexstr = hexstr.strip('#')

    if scalefactor < 0 or len(hexstr) != 6:
        return hexstr

    r, g, b = int(hexstr[:2], 16), int(hexstr[2:4], 16), int(hexstr[4:], 16)

    r = clamp(r * scalefactor)
    g = clamp(g * scalefactor)
    b = clamp(b * scalefactor)

    return "#%02x%02x%02x" % (r, g, b)


#############
##main method

#default parameters
left=300
right=300
normalize=False
gsigma=1 # gaussian smoothing parameter

filebase=None # such as "pearson-50-constgap-trnaall_aligned_composite"
try:
    flistname=sys.argv[1]
except IndexError:
    print("please provide path to file list")
    pass

try:
    left=int(sys.argv[2])
except IndexError:
    pass

try:
    right=int(sys.argv[3])
except IndexError:
    pass

try:
    words=sys.argv[4]
    if "normalize" in words:
	normalize=True
except IndexError:
    pass


#output file name
pngoutname="chexalign-out_composite.png"
#list of experiments to take
flist=[]
factorname=[]
with open(flistname) as f:
    for line in f:
	flist.append(line.rstrip())
	factorname.append(line.rstrip().split(".")[1])
nof=len(flist)

#color list
print(random.choice([i for i in '0123456789ABCDEF']))
colortest=[random.choice([i for i in '0123456789ABCDEF']) for j in range(6)]
print(colortest)

color_lst=["#"+''.join([random.choice([i for i in '0123456789ABCDEF']) for j in range(6)]) for i in range(nof)]

#get xaxis from the first file
xaxis,y=loadtxt(flist[0],usecols=(0,1),unpack=True)
pyaxis=[None]*nof
nyaxis=[None]*nof
maxv=[None]*nof

#read files
for i in range(0,nof):
    py,tmpy = loadtxt(flist[i],usecols=(1,2),unpack=True)
    ny = [-float(j) for j in tmpy]
    pyaxis[i]=py
    nyaxis[i]=ny

#smooth
for i in range(0,nof):
    pyaxis[i]=ndimage.gaussian_filter(pyaxis[i],sigma=gsigma)
    nyaxis[i]=ndimage.gaussian_filter(nyaxis[i],sigma=gsigma)
    print(pyaxis[i])
    if normalize:
	htsum=(np.sum(pyaxis[i])-np.sum(nyaxis[i]))/2
	print(htsum)
	print(pyaxis[i])
	pyaxis[i]=pyaxis[i]/htsum
	nyaxis[i]=nyaxis[i]/htsum
	print(pyaxis[i])
    maxv[i]=max(max(pyaxis[i]),max(-nyaxis[i]))*1.1


plt.ioff()
fig, ax = plt.subplots(nof, sharex=True, figsize=(7,15))

for i in range(0,nof):
    ax[i].fill_between(xaxis, 0, pyaxis[i],facecolor=color_lst[i])
    ax[i].fill_between(xaxis, 0, nyaxis[i],facecolor=color_lst[i],alpha=0.6)
    ax[i].plot(xaxis,pyaxis[i],lw=3,color=color_lst[i])
    ax[i].plot(xaxis,nyaxis[i],lw=3,color=color_lst[i])
    plt.draw()
    ax[i].set_xlim([-left,right])
    ax[i].set_ylim([-maxv[i],maxv[i]])
    ax[i].tick_params(axis='both', which='major', labelsize=16)
    ax[i].set_ylabel(factorname[i],fontsize=20)

plt.xlabel("Distance to alignment center",fontsize=20)
#plt.ylabel("Normalized tags",fontsize=20)
plt.xticks(size=16)
plt.tight_layout()
fig.subplots_adjust(hspace=0)
#plt.tick_params(labelsize=25)

plt.savefig(pngoutname)
plt.close(fig)

