"""
Produce strand separated heatmaps from aligned profiles
Input .mat format where each row represents different genomic regions and column is a base pair position tag count values
Gaps are represented as a negative value
"""

import glob
import operator
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from colour import Color
import scipy
from scipy import ndimage
import sys
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker

filename = sys.argv[1]

maxv=20 # defualt max value for color saturation
try:
     maxv=float(sys.argv[2])
except IndexError:
    pass

hasRC=""
try:
     hasRC=sys.argv[3]
except IndexError:
    pass

isRC=False
if "rc" in hasRC:
    isRC=True


trimp=0.2 # trim heatmap if less than trimp data available

txt = np.genfromtxt(filename,dtype=None)
wdata=[]
cdata=[]
for rec in txt:
    d = [float(x) if float(x)>=0 else -maxv for x in rec[1].split(",")[:-1]]
    if len(d) > 2:
	if ("watson") in rec[0]:
	    wdata.append(d)
	elif ("crick") in rec[0]:
	    cdata.append(d)	


#make data for display
nrows=len(wdata)
ncols=len(wdata[0])
print("number of rows:"+str(nrows))
print("number of colomns:"+str(ncols))
tfarray=nrows*([True]*ncols,[False]*ncols)

##masked arrays are not used in plots
#wmask=np.array(nrows*([True]*ncols,[False]*ncols), dtype=bool)
#w_data=np.ma.masked_where(wmask,np.repeat(wdata,2,axis=0))

#cmask=np.array(nrows*([False]*ncols,[True]*ncols), dtype=bool)
#c_data=np.ma.masked_where(cmask,np.repeat(cdata,2,axis=0))
#print("my array wdata")
#print(wdata)
#print("cdata")
#print(cdata)

xaxis=[]
for i in range(-ncols/2,ncols/2):
    xaxis.append(i)

#decide what range to plot
#trim data from left
left=0
right=ncols-1
for i in range(0,ncols):
    blank=0
    for j in range(0,nrows):
	if wdata[j][i]<0:
	    blank+=1
    if float(blank)/float(nrows) > (1-trimp):
	left=i
    else:
	break;
#trim data from right
for i in range(ncols-1,-1,-1):
    blank=0
    for j in range(0,nrows):
	if wdata[j][i]<0:
	    blank+=1
    if float(blank)/float(nrows) > (1-trimp):
	right=i
    else:
	break;
print("left is "+str(left))
print("right is "+str(right))

#convert list to numpy arrays
npwdata=np.array(wdata)
npcdata=np.array(cdata)

#color
bcmap = LinearSegmentedColormap.from_list(name="costom",colors =['grey','w',(0,0,1)],N=100) # blue
rcmap = LinearSegmentedColormap.from_list(name="costom",colors =['grey','w',(1,0,0)],N=100) # red

#plot
plt.ioff()

fig, ax = plt.subplots(figsize=(30,30))

if isRC:
    heatmap = sns.heatmap(npwdata[:,right:left:-1],vmin=-maxv,vmax=maxv,cmap=rcmap)
    heatmap = sns.heatmap(npcdata[:,right:left:-1],alpha = 0.4,vmin=-maxv,vmax=maxv,cmap=bcmap)
else:
    heatmap = sns.heatmap(npcdata[:,left:right],vmin=-maxv,vmax=maxv,cmap=rcmap)
    heatmap = sns.heatmap(npwdata[:,left:right],alpha = 0.4,vmin=-maxv,vmax=maxv,cmap=bcmap)
texts = [t.get_text()  for t in ax.get_xticklabels()]
freq=int(texts[1])-int(texts[0])
if isRC:
    ax.set_xticklabels(xaxis[right:left:-1][0::freq])
    plt.savefig(filename.replace(".mat",".rv.png"))
else:
    ax.set_xticklabels(xaxis[left:right][0::freq])
    plt.savefig(filename.replace(".mat",".png"))


plt.close()
