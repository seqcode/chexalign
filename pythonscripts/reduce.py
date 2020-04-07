import glob
import operator
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import scipy
from scipy import ndimage
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.manifold import TSNE
import sys

#run pca
def pca_red(X):
    pca = PCA(n_components=1)
    pca.fit(X)
    print "variance explained by pca"
    print pca.explained_variance_ratio_
    return pca.transform(X)

#run mds
def mds_red(X):
    mds = MDS(n_components=1)
    return mds.fit_transform(X)

#run tsne
def tsne_red(X):
    tsne = TSNE(n_components=2)
    return tsne.fit_transform(X)

#make plots
def scatter_plot(data,labels,outname):
    plt.ioff()
    fig=plt.figure(figsize=(8,6))
    axes=plt.gca()
    plt.xticks(fontsize=15)
    plt.xlabel("PC1",fontsize=20)
    axes.set_yticklabels([])
    plt.scatter(data,[0]*len(data),s=70)
    for label, x in zip(labels,data):
    	plt.annotate(label,xy=(x,0),xytext=(-5,5),textcoords='offset points')
    plt.savefig(outname)
    plt.close()

#print input for the dimensionality reduction
def printmat(indices,faclst,inmat,outname):
    outs="#\t"
    for ind in indices:
    	outs=outs+str(ind)+"\t"
    outs=outs+"\n"
    for i in range(0,len(inmat)):
    	slist=[str(x) for x in inmat[i]]
    	outs=outs+faclst[i]+"\t"+"\t".join(slist)+"\n"
    with open(outname,'w')as mout:
   	mout.write(outs)

############## 
#main


try:
    flistname=sys.argv[1]
except IndexError:
    print("please provide path to file list")
    pass

flist=[]
with open(flistname) as f:
    for line in f:
        flist.append(line.rstrip())

all_indices=[]
results={}
for f in flist:
    with open(f,'r') as mlf:
	firstline=mlf.readline().rstrip().split("\t")
	for words in firstline:
	    if "XL" in words:
		ind=words.split(',')[1]
		all_indices.append(ind)

tmp_l=list(set(all_indices))
temp_num=sorted([int(x) for x in tmp_l])
indices=[str(x) for x in temp_num]


for f in flist:
    factor=f.split(".")[1]
    d={} #dictionary contains position mapping to column number
    # get column index for positions
    with open(f,'r') as mlf:
	firstline=mlf.readline().rstrip().split("\t")
	for index in indices:
	    for i in range(0,len(firstline)):
		if "XL" in firstline[i]:
		    currcol=firstline[i].split(',')[1]
	            if str(index) == currcol:
		    	d[index]=i
	xl={}
	for index in indices:
	    xl[index]=0
	for line in mlf:
	    back=line.rstrip().split("\t")[2]
	    for index in indices:
		if index in d:
		    xl[index]+=float(line.rstrip().split("\t")[d[index]])
		else:
                    xl[index]+=0
	results[factor]=xl


faclst=[]
mat=[]
for key in results:
    faclst.append(key)
    val=results[key]
    flst=[]
    for i in indices:
	flst.append(val[i])
    mat.append(flst)

#normalize 
nmat=[]
for row in mat:
    nrowmat=[]
    rowsum=0
    for i in row:
	rowsum=rowsum+i
    for i in row:
	nrowmat.append(i/rowsum)
    nmat.append(nrowmat)


#print input matrix for dimentionality reduction
#printmat(indices,faclst,mat,'mat.out')
printmat(indices,faclst,nmat,'nmat.out')

#perform dimentionality reduction
x=np.array(nmat)
pca_mat=pca_red(x)
mds_mat=mds_red(x)
tsne_mat=tsne_red(x)

#make plots
scatter_plot(pca_mat,faclst,"pca_test.png")
scatter_plot(mds_mat,faclst,"mds_test.png")
