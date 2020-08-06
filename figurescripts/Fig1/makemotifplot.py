"""
Plot motif scores based on aligned tag profile
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

filename = sys.argv[1]

color="green" # defualt max value for color saturation
try:
     color=sys.argv[2]
except IndexError:
    pass

maxv=-1
try:
     maxv=float(sys.argv[3])
except IndexError:
    pass

hasRC=""
try:
     hasRC=sys.argv[4]
except IndexError:
    pass

isRC=False
if "rc" in hasRC:
    isRC=True

maxdefined=True
if maxv==-1:
    maxdefined=False

#1-fraction of data in a column required for plotting
trimp=0.2

negvalue=-18
txt = np.genfromtxt(filename,dtype=None)
data=[]
for rec in txt:
    lst=rec.rstrip().split(",")
    d=[]
    for x in lst:
	if x:
	    if float(x)>=0:
		d.append(float(x))
	    else:
		d.append(negvalue)
    if not maxdefined:
    	currmax=max(d)
    	if currmax > maxv:
	    maxv=currmax
    data.append(d)

#replace maxvalue
#for d in data:
#    for n, i in enumerate(d):
#	if i==negvalue:
#	    d[n]=-maxv

#trim data from left
nrows=len(data)
ncols=len(data[0])

left=0
right=ncols-1
for i in range(0,ncols):
    blank=0
    for j in range(0,nrows):
	if data[j][i]<0:
	    blank+=1
    if float(blank)/float(nrows) > (1-trimp):
	left=i
    else:
	break;
#trim data from right
for i in range(ncols-1,-1,-1):
    blank=0
    for j in range(0,nrows):
	if data[j][i]<0:
	    blank+=1
    if float(blank)/float(nrows) > (1-trimp):
	right=i
    else:
	break;
print("left is "+str(left))
print("right is "+str(right))

#convert to numpy data
npdata=np.array(data)

fcounts=np.sum(npdata, axis=0)
center=np.argmax(fcounts)
print("center: "+str(center))

totalc=0
tenbp=0
for row in npdata:
    for i in range(left,right):
	if row[i] > 0:
	    totalc+=1
    for i in range(center-10,center+10):
	if row[i] > 0:
	    tenbp+=1
print("totalc: "+str(totalc)+".counts at center: "+str(tenbp))

print("max val: "+str(maxv))

#color
bcmap = LinearSegmentedColormap.from_list(name="costom",colors =['grey','w',(0,0,1)],N=100) # blue
rcmap = LinearSegmentedColormap.from_list(name="costom",colors =['grey','w',(1,0,0)],N=100) # red
gcmap = LinearSegmentedColormap.from_list(name="costom",colors =['grey','w',(0,0.7,0)],N=100) # green

cmap=None
if color=="red":
    cmap =rcmap
elif color=="blue":
    cmap = bcmap
else:
    cmap = gcmap    

#plot
plt.ioff()
fig, ax = plt.subplots(figsize=(30,30))
if isRC:
    heatmap = sns.heatmap(npdata[:,right:left:-1],vmin=-maxv,vmax=maxv,cmap=cmap)
    plt.savefig(filename+"_motif.rc.png")
else:
    heatmap = sns.heatmap(npdata[:,left:right],vmin=-maxv,vmax=maxv,cmap=cmap)
    plt.savefig(filename+"_motif.png")


