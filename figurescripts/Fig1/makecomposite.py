import glob
import operator
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns; sns.despine()
from colour import Color
import scipy
from scipy import ndimage
import sys

filename = sys.argv[1]
hasRC=""
try:
     hasRC=sys.argv[2]
except IndexError:
    pass

isRC=False
if "rc" in hasRC:
    isRC=True


trimp=0.2

txt = np.genfromtxt(filename,dtype=None)
w_data=[]
c_data=[]


for rec in txt:
    d = [float(x) if float(x)>=0 else -1 for x in rec[1].split(",")[:-1]]
    if len(d) > 2:
    	if "watson" in rec[0]:
	    w_data.append(d)
	else:
	    c_data.append(d)
	

#decide what range to plot
#trim data from left
nrows=len(c_data)
ncols=len(c_data[0])

left=0
right=ncols-1
for i in range(0,ncols):
    blank=0
    for j in range(0,nrows):
	if w_data[j][i]<0:
	    blank+=1
    if float(blank)/float(nrows) > (1-trimp):
	left=i
    else:
	break;
#trim data from right
for i in range(ncols-1,-1,-1):
    blank=0
    for j in range(0,nrows):
	if w_data[j][i]<0:
	    blank+=1
    if float(blank)/float(nrows) > (1-trimp):
	right=i
    else:
	break;
print("left is "+str(left))
print("right is "+str(right))

#replace -1 to 0
for data in w_data:
    for i, n in enumerate(data):
	if n==-1:
	    data[i]=0
for data in c_data:
    for i, n in enumerate(data):
        if n==-1:
            data[i]=0

w_composite=w_data[0]
for num in range(1,len(w_data)):
    for c, value in enumerate(w_data[num]):
	w_composite[c]+=value
c_composite=c_data[0]
for num in range(1,len(c_data)):
    for c, value in enumerate(c_data[num]):
        c_composite[c]+=value


xaxis=[]
for i in range(-ncols/2,ncols/2):
    xaxis.append(i)

print(len(xaxis))
print(len(w_composite[left:right]))


#smooth
#w_composite=ndimage.gaussian_filter(w_composite,sigma=1)
#c_composite=ndimage.gaussian_filter(c_composite,sigma=1)

plt.ioff()

fig, ax = plt.subplots(figsize=(16,8))

if isRC:
    tmpy=w_composite
    w_composite = [-int(j) for j in tmpy]
    ax.fill_between(xaxis[left:right], 0, c_composite[right:left:-1],facecolor="blue")
    ax.fill_between(xaxis[left:right], 0, w_composite[right:left:-1],facecolor="red")
    ax.plot(xaxis[left:right],c_composite[right:left:-1],lw=2,color="blue")
    ax.plot(xaxis[left:right],w_composite[right:left:-1],lw=2,color="red")
else:
    tmpy=c_composite
    c_composite = [-int(j) for j in tmpy]
    ax.fill_between(xaxis[left:right], 0, w_composite[left:right],facecolor="blue")
    ax.fill_between(xaxis[left:right], 0, c_composite[left:right],facecolor="red")
    ax.plot(xaxis[left:right],w_composite[left:right],lw=2,color="blue")
    ax.plot(xaxis[left:right],c_composite[left:right],lw=2,color="red")


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_linewidth(4)
ax.yaxis.set_tick_params(width=4)

plt.margins(x=0)

if isRC:
    plt.savefig(filename.replace(".mat","_composite.rv.png"))
else:
    plt.savefig(filename.replace(".mat","_composite.png"))

plt.close(fig)
