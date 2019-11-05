package org.seqcode.projects.chexalign.alignment;

import org.seqcode.genome.location.StrandedPoint;

/**
 * TagProfile is used in to represent individual aligned or unaligned (original) tag profiles.
 * 
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class TagProfile implements Comparable<TagProfile> {
	protected int aLength;
	protected int numCond;
	protected double[][] watsonProfile;			//{Condition; alignment length}
	protected double[][] crickProfile;			//{Condition; alignment length}
	protected double[] controlWatsonProfile;	//{alignment length}
	protected double[] controlCrickProfile;		//{alignment length}
	protected StrandedPoint originalPt=null;	// Original input point
	protected int[] coords;						// Integer arrays indicating positions
	protected int[] gap;
	protected boolean reverse=false;
	protected int pIndex;						// index of point of origin
	protected int id;
	protected int numMembers;
	
	public TagProfile(int al, int nC){
		aLength=al;
		numCond=nC;
		watsonProfile = new double[numCond][aLength];
		crickProfile = new double[numCond][aLength];
		controlWatsonProfile = new double[aLength];
		controlCrickProfile = new double[aLength];
		gap = new int[aLength];
		coords = new int[aLength];
		for (int i=0; i < al ; i++){
			gap[i]=0;
			coords[i]=-1;
		}
	}
	
	public TagProfile(int al, int nC, double[][] watsonProfile, double[][] crickProfile, double[] controlWatsonProfile, double[] controlCrickProfile){
		aLength=al;
		numCond=nC;
		this.watsonProfile = watsonProfile;
		this.crickProfile = crickProfile;
		this.controlWatsonProfile = controlWatsonProfile;
		this.controlCrickProfile = controlCrickProfile;
		gap = new int[aLength];
		coords = new int[aLength];
		for (int i=0; i < al ; i++){
			gap[i]=0;
			coords[i]=-1;
		}
	}
	
	// Accessors
	public double[][] getWatsonProfile(){return watsonProfile;}
	public double[][] getCrickProfile(){return crickProfile;}
	public double[] getControlWatsonProfile(){return controlWatsonProfile;}
	public double[] getControlCrickProfile(){return controlCrickProfile;}
	public int[] getGaps(){return gap;}
	public int[] getCoords(){return coords;}
	public boolean getReverse(){return reverse;}
	public StrandedPoint getStrandedPoint(){return originalPt;}
	public int getIndex(){return pIndex;}
	public int getID(){return id;}
	public int getNumMembers(){return numMembers;}
	public int getAlignLength(){return aLength;}
	
	// Setters
	public void setWatsonProfile(int cond, int index, double val){watsonProfile[cond][index]=val;}
	public void setCrickProfile(int cond, int index, double val){crickProfile[cond][index]=val;}
	public void setControlWatsonProfile(int index, double val){controlWatsonProfile[index]=val;}
	public void setControlCrickProfile(int index, double val){controlCrickProfile[index]=val;}
	public void setCoord(int index, int pos){coords[index]=pos;}
	public void setGapProfile(int[] gaps){this.gap=gaps;}
	public void setCoordinates(int[] positions){coords = positions;}
	public void setStrandedPoint(StrandedPoint pt){originalPt=pt;}
	public void isReverse(){reverse=true;}
	public void setID(int id){this.id=id;}
	public void setIndex(int index){pIndex=index;}
	public void setNumMembers(int numM){numMembers=numM;}
	
	public void setInitialCoords(StrandedPoint pt){
		originalPt=pt;
		if(pt.getStrand()=='+')
			for (int w=0; w < aLength; w++)
				coords[w] = pt.getLocation()-aLength/2+w;
		else
			for (int w=0; w < aLength; w++)
				coords[aLength-w-1]=pt.getLocation()-aLength/2+w;
	}
	
	//Comparable default method
	public int compareTo(TagProfile profile) {
		return (this.id - profile.id);
	}
	
	@Override
	public String toString(){
		return "[id="+this.id+"]";
	}
	
}
