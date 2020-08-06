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

filebase=None # such as "pearson-50-constgap-trnaall_aligned_composite"
try:
    filebase=sys.argv[1]
except IndexError:
    print("please provide file base name")
    pass

try:
    outfname=sys.argv[2]
except IndexError:
    print("please provide file base name")
    pass

left=300
right=300
normalize=False


try:
    left=int(sys.argv[3])
except IndexError:
    pass

try:
    right=int(sys.argv[4])
except IndexError:
    pass

try:
    words=sys.argv[5]
    if "normalize" in words:
	normalize=True
except IndexError:
    pass

gsigma=1 # gaussian smoothing parameter


#output file name
pngoutname=outfname+".png"
#list of experiments to take
explist=['Rap1','Sfp1','Ifh1','Fhl1','Hmo1']
#color list
color_lst=["red","grey","green","purple","blue"]
ncolor_lst=[]
for col in color_lst:
    ncolor_lst.append(colorscale(matplotlib.colors.cnames[col],1.5))

flist=[]
for item in explist:
    flist.append(filebase+"."+item+".txt")


print(flist)

nof=len(flist)

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

print(xaxis)
print(pyaxis[0])
print(nyaxis[0])

plt.ioff()
fig, ax = plt.subplots(nof, sharex=True, figsize=(7,15))

for i in range(0,nof):
    ax[i].fill_between(xaxis, 0, pyaxis[i],facecolor=color_lst[i])
    ax[i].fill_between(xaxis, 0, nyaxis[i],facecolor=color_lst[i],alpha=0.6)
    ax[i].plot(xaxis,pyaxis[i],lw=3,color=color_lst[i])
    ax[i].plot(xaxis,nyaxis[i],lw=3,color=color_lst[i])
#    for j,label in enumerate(ax[i].get_yticklabels(which='both')):
#	if not j==1 and not j==len(ax[i].get_yticklabels()) - 2:
#	    label.set_visible(False)
#	else:
#	    label.set_fontsize(25)
    plt.draw()
    ax[i].set_xlim([-left,right])
    ax[i].set_ylim([-maxv[i],maxv[i]])

plt.tight_layout()
fig.subplots_adjust(hspace=0)
plt.tick_params(labelsize=25)
#plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

plt.savefig(pngoutname)
plt.close(fig)

