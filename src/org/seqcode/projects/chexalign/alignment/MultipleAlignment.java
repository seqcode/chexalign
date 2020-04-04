package org.seqcode.projects.chexalign.alignment;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Pair;

/**
 * Stores methods for multiple alignment
 */
public class MultipleAlignment {
	
	protected int numAligned;
	protected int alignL;
	protected int numCond;
	protected StrandedPoint[] alignedPoints;
	protected int[] alignedIDs;
	protected TagProfile[] profileAlignment;
	protected boolean hasControl;
	
	public MultipleAlignment(int nA, int aL, int numC, boolean hasControl){
		numAligned=nA;
		alignL=aL;
		numCond=numC;
		profileAlignment = new TagProfile[numAligned];
		for (int i=0; i < numAligned; i++)
			profileAlignment[i] = new TagProfile(alignL, numCond);
		alignedIDs = new int[numAligned];
		this.hasControl=hasControl;
	}
	
	//Accessors
	public int getAlignedLength(){return alignL;}
	public int getNumAligned(){return numAligned;}
	public TagProfile getTagProfile(int index){return profileAlignment[index];}
	
	// Setters
	public void setTagProfile(int index, TagProfile tagProfile){ profileAlignment[index]=tagProfile;}

	public TagProfile Alignment2Profile(){		
		double[][] watsonProfile = new double[numCond][alignL];
		double[][] crickProfile = new double[numCond][alignL];
		double[] controlWatsonProfile = null;
		double[] controlCrickProfile = null;
		int[] gap = new int[alignL];
		for (int z=0; z < alignL; z++){
			for (int c=0; c < numCond; c++){
				watsonProfile[c][z]=0;
				crickProfile[c][z]=0;		
			}	
		}
		
		if (hasControl){
			controlWatsonProfile = new double[alignL];
			controlCrickProfile = new double[alignL];
			for (int z=0; z < alignL; z++){
				controlWatsonProfile[z]=0;
				controlCrickProfile[z]=0;	
			}
		}
		
		for (int z=0; z < alignL; z++){
			for (int x=0; x < numAligned; x++){
				if (profileAlignment[x].getWatsonProfile()[0][z]==-1){
					gap[z]++;
				}else{
					if (hasControl){
						controlWatsonProfile[z]+=profileAlignment[x].getControlWatsonProfile()[z];
						controlCrickProfile[z]+=profileAlignment[x].getControlCrickProfile()[z];
					}
					for (int c=0; c < numCond; c++){
						watsonProfile[c][z]+=profileAlignment[x].getWatsonProfile()[c][z];
						crickProfile[c][z]+=profileAlignment[x].getCrickProfile()[c][z];
					}}}}
		
		TagProfile newProfile= new TagProfile(alignL, numCond, watsonProfile, crickProfile,controlWatsonProfile,controlCrickProfile);
		newProfile.setGapProfile(gap);
		newProfile.setNumMembers(numAligned);
		
		return newProfile;
	}
	
			
	/**
     * Add a single region tag count to a set of aligned tags
     * @param NeedlemanWunschAffine matrix
     * @param Tag profile two
     * @param TagProfile two ID
     * @return New Multiple Alignment
     */
	public MultipleAlignment SingleProfileAddition(NeedlemanWunschAffine nw, TagProfile two){
				
		int aL=nw.getAlignedLength();	// alignment length of new alignment
		
		MultipleAlignment newAlignment = new MultipleAlignment(numAligned+1, aL, numCond, hasControl);
		
		// Set IDs and point location for a new alignment
		for (int a=0; a < numAligned; a++){
			newAlignment.getTagProfile(a).setID(profileAlignment[a].getID());
			newAlignment.getTagProfile(a).setStrandedPoint(profileAlignment[a].getStrandedPoint());
		}
		newAlignment.getTagProfile(numAligned).setID(two.getID());
		newAlignment.getTagProfile(numAligned).setStrandedPoint(two.getStrandedPoint());
		
		// Update exisiting alignment
		ArrayList<Integer> alignSectionOne = nw.getAlignSectionX();	//	aligned section of a set of existing alignment
		ArrayList<Integer> alignSectionTwo = nw.getAlignSectionY();	//	aligned section of a new profile
		
		int last1=-50; int last2=-50;
		int antiZ=0;
		for (int z=aL-1; z>=0; z--){
			if(alignSectionTwo.get(z)==last2 || alignSectionTwo.get(z)==-1){
				// Gap in alignment 2; add in alignment 1's column  only
				for (int a=0; a < numAligned; a++){
					newAlignment.getTagProfile(a).setCoord(antiZ, profileAlignment[a].getCoords()[alignSectionOne.get(z)]);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(a).setWatsonProfile(c, antiZ, profileAlignment[a].getWatsonProfile()[c][alignSectionOne.get(z)]);
						newAlignment.getTagProfile(a).setCrickProfile(c, antiZ, profileAlignment[a].getCrickProfile()[c][alignSectionOne.get(z)]);
					}
					//control experiment
					if (hasControl){
						newAlignment.getTagProfile(a).setControlWatsonProfile(antiZ, profileAlignment[a].getControlWatsonProfile()[alignSectionOne.get(z)]);
						newAlignment.getTagProfile(a).setControlCrickProfile(antiZ, profileAlignment[a].getControlCrickProfile()[alignSectionOne.get(z)]);
					}
				}
				newAlignment.getTagProfile(numAligned).setCoord(antiZ, -1);
				for (int c=0; c < numCond; c++){
					newAlignment.getTagProfile(numAligned).setWatsonProfile(c, antiZ, -1);
					newAlignment.getTagProfile(numAligned).setCrickProfile(c, antiZ, -1);
				}
				if (hasControl){
					newAlignment.getTagProfile(numAligned).setControlWatsonProfile(antiZ, -1);
					newAlignment.getTagProfile(numAligned).setControlCrickProfile(antiZ, -1);
				}
			}else if(alignSectionOne.get(z)==last1 || alignSectionOne.get(z)==-1){
				// Gap in alignment 1; add in alignment 2's column  only
				for (int a=0; a < numAligned; a++){
					newAlignment.getTagProfile(a).setCoord(antiZ, -1);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(a).setWatsonProfile(c, antiZ, -1);
						newAlignment.getTagProfile(a).setCrickProfile(c, antiZ, -1);
					}
					if (hasControl){
						newAlignment.getTagProfile(a).setControlWatsonProfile(antiZ, -1);
						newAlignment.getTagProfile(a).setControlCrickProfile(antiZ, -1);
					}
				}				
				if (nw.isReverse()){	// reverse direction
					newAlignment.getTagProfile(numAligned).setCoord(antiZ,two.getCoords()[two.getAlignLength()-alignSectionTwo.get(z)-1]);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(numAligned).setWatsonProfile(c, antiZ, two.getCrickProfile()[c][two.getAlignLength()-alignSectionTwo.get(z)-1]);
						newAlignment.getTagProfile(numAligned).setCrickProfile(c, antiZ, two.getWatsonProfile()[c][two.getAlignLength()-alignSectionTwo.get(z)-1]);
					}
					if (hasControl){
						newAlignment.getTagProfile(numAligned).setControlWatsonProfile(antiZ, two.getControlCrickProfile()[two.getAlignLength()-alignSectionTwo.get(z)-1]);
						newAlignment.getTagProfile(numAligned).setControlCrickProfile(antiZ, two.getControlWatsonProfile()[two.getAlignLength()-alignSectionTwo.get(z)-1]);
					}
				}else{	// forward direction
					newAlignment.getTagProfile(numAligned).setCoord(antiZ,two.getCoords()[alignSectionTwo.get(z)]);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(numAligned).setWatsonProfile(c, antiZ, two.getWatsonProfile()[c][alignSectionTwo.get(z)]);
						newAlignment.getTagProfile(numAligned).setCrickProfile(c, antiZ, two.getCrickProfile()[c][alignSectionTwo.get(z)]);
					}
					if (hasControl){
						newAlignment.getTagProfile(numAligned).setControlWatsonProfile(antiZ,two.getControlWatsonProfile()[alignSectionTwo.get(z)]);
						newAlignment.getTagProfile(numAligned).setControlCrickProfile(antiZ, two.getControlCrickProfile()[alignSectionTwo.get(z)]);
					}
				}
			}else{
				// No gap; add in both alignments
				for (int a=0; a < numAligned; a++){	
					newAlignment.getTagProfile(a).setCoord(antiZ, profileAlignment[a].getCoords()[alignSectionOne.get(z)]);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(a).setWatsonProfile(c, antiZ, profileAlignment[a].getWatsonProfile()[c][alignSectionOne.get(z)]);
						newAlignment.getTagProfile(a).setCrickProfile(c, antiZ, profileAlignment[a].getCrickProfile()[c][alignSectionOne.get(z)]);
					}
					if (hasControl){
						newAlignment.getTagProfile(a).setControlWatsonProfile(antiZ, profileAlignment[a].getControlWatsonProfile()[alignSectionOne.get(z)]);
						newAlignment.getTagProfile(a).setControlCrickProfile(antiZ, profileAlignment[a].getControlCrickProfile()[alignSectionOne.get(z)]);
					}
				}
				if (nw.isReverse()){	// reverse direction
					newAlignment.getTagProfile(numAligned).setCoord(antiZ,two.getCoords()[two.getAlignLength()-alignSectionTwo.get(z)-1]);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(numAligned).setWatsonProfile(c, antiZ, two.getCrickProfile()[c][two.getAlignLength()-alignSectionTwo.get(z)-1]);
						newAlignment.getTagProfile(numAligned).setCrickProfile(c, antiZ, two.getWatsonProfile()[c][two.getAlignLength()-alignSectionTwo.get(z)-1]);	
					}
					if (hasControl){
						newAlignment.getTagProfile(numAligned).setControlWatsonProfile(antiZ, two.getControlCrickProfile()[two.getAlignLength()-alignSectionTwo.get(z)-1]);
						newAlignment.getTagProfile(numAligned).setControlCrickProfile(antiZ, two.getControlWatsonProfile()[two.getAlignLength()-alignSectionTwo.get(z)-1]);	
					}
				}else{	// forward direction
					newAlignment.getTagProfile(numAligned).setCoord(antiZ, two.getCoords()[alignSectionTwo.get(z)]);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(numAligned).setWatsonProfile(c, antiZ, two.getWatsonProfile()[c][alignSectionTwo.get(z)]);
						newAlignment.getTagProfile(numAligned).setCrickProfile(c, antiZ, two.getCrickProfile()[c][alignSectionTwo.get(z)]);
					}
					if (hasControl){
						newAlignment.getTagProfile(numAligned).setControlWatsonProfile(antiZ, two.getControlWatsonProfile()[alignSectionTwo.get(z)]);
						newAlignment.getTagProfile(numAligned).setControlCrickProfile(antiZ, two.getControlCrickProfile()[alignSectionTwo.get(z)]);
					}
				}
			}
			last1=alignSectionOne.get(z);
			last2=alignSectionTwo.get(z);
			antiZ++;
		}
		
		return newAlignment;
	}
	
	/**
     * Subtract a single profile from a set of alignment
     * @param index of subtracting profile
     * @return New Multiple Alignment
     */
	public MultipleAlignment SingleProfileSubstraction(int removeID){
		// Find the row to remove
		int removeRow=-1;
		for (int i=0; i < numAligned; i ++)
			if (profileAlignment[i].getID()==removeID)
				removeRow=i;	
		if (removeRow==-1){
			System.err.println("Error in Iterative Refinement Multiple Alignment: this profile is not in the current alignment\n");
			System.exit(1);
		}
		
		//Count the columns that have only gaps in non-removeRow positions
		int numBlankCol=0;
		int gapCount=0;
		for (int j=0; j < alignL; j++){
			gapCount=0;
			for (int i=0; i< numAligned; i++){
				if (i!=removeRow && profileAlignment[i].getWatsonProfile()[0][j]==-1)
					gapCount++;
			if (gapCount==numAligned-1)
				numBlankCol++;
			}
		}
		//declare the new alignment
		int newAL=alignL-numBlankCol;			
		MultipleAlignment newAlignment = new MultipleAlignment(numAligned-1, newAL, numCond, hasControl);
		// Copy the relevant information into the new alignment
		int a=0;
		for (int i=0; i < numAligned; i++){
			if (i!=removeRow){
				newAlignment.getTagProfile(a).setID(profileAlignment[i].getID());
				newAlignment.getTagProfile(a).setStrandedPoint(profileAlignment[a].getStrandedPoint());
				a++;
			}
		}
		a=0;		// a counts new rows
		int z=0;	//z counts new columns
		for (int j=0; j < alignL; j++){
			//Is this an empty column?
			gapCount=0;
			a=0;
			for (int i=0; i < numAligned; i++){
				if (i!=removeRow && profileAlignment[i].getWatsonProfile()[0][j]==-1)
					gapCount++;
			}
			if (gapCount < numAligned-1){ //Not an empty column, add the cells in this column
				for (int i=0; i < numAligned; i++){
					if (i!=removeRow){
						newAlignment.getTagProfile(a).setCoord(z, profileAlignment[i].getCoords()[j]);
						for (int c=0; c< numCond; c++){
							newAlignment.getTagProfile(a).setWatsonProfile(c, z, profileAlignment[i].getWatsonProfile()[c][j]);
							newAlignment.getTagProfile(a).setCrickProfile(c, z, profileAlignment[i].getCrickProfile()[c][j]);
						}
						if (hasControl){
							newAlignment.getTagProfile(a).setControlWatsonProfile(z, profileAlignment[i].getControlWatsonProfile()[j]);
							newAlignment.getTagProfile(a).setControlCrickProfile(z, profileAlignment[i].getControlCrickProfile()[j]);
						}
					a++;	
					}}
				z++;
			}
		}
		return newAlignment;		
	}
	
	/**
     * Make aligned tags based on multiple profile addition
     * @param NeedlemanWunschAffine matrix
     * @param a set of alignment A
     * @param a set of alignment B
     * @param anchor index
     * @return mutated array of Alignment[] alignment
     */
	public MultipleAlignment MultipleProfileAddition(NeedlemanWunschAffine nw, MultipleAlignment alignmentTwo){
		int aL=nw.getAlignedLength();	// alignment length of new alignment
		
		MultipleAlignment newAlignment = new MultipleAlignment(numAligned+alignmentTwo.getNumAligned(), aL, numCond, hasControl);
		
		// Copy and set IDs for new alignment
		for (int a=0; a<numAligned; a++){
			newAlignment.getTagProfile(a).setID(profileAlignment[a].getID());
			newAlignment.getTagProfile(a).setStrandedPoint(profileAlignment[a].getStrandedPoint());
		}
		for (int b=0; b < alignmentTwo.getNumAligned(); b++){
			newAlignment.getTagProfile(numAligned+b).setID(alignmentTwo.getTagProfile(b).getID());
			newAlignment.getTagProfile(numAligned+b).setStrandedPoint(alignmentTwo.getTagProfile(b).getStrandedPoint());
		}
		
		// Update existing alignment
		ArrayList<Integer> alignSectionOne = nw.getAlignSectionX();	//	aligned section of a set of existing alignment one
		ArrayList<Integer> alignSectionTwo = nw.getAlignSectionY();	//	aligned section of a set of existing alignment two
		
		int last1=-50; int last2=-50;
		int antiZ=0;
		double currWatsonVal=0; double currCrickVal=0;
		double currControlWatsonVal=0; double currControlCrickVal=0;
		int coord=-1;
		for (int z=aL-1; z>=0; z--){
			if(alignSectionTwo.get(z)==last2 || alignSectionTwo.get(z)==-1){
				// Gap in alignment 2; add in alignment 1's column  only
				for (int a=0; a < numAligned; a++){
					newAlignment.getTagProfile(a).setCoord(antiZ, profileAlignment[a].getCoords()[alignSectionOne.get(z)]);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(a).setWatsonProfile(c, antiZ, profileAlignment[a].getWatsonProfile()[c][alignSectionOne.get(z)]);
						newAlignment.getTagProfile(a).setCrickProfile(c, antiZ, profileAlignment[a].getCrickProfile()[c][alignSectionOne.get(z)]);
					}
					if (hasControl){
						newAlignment.getTagProfile(a).setControlWatsonProfile(antiZ, profileAlignment[a].getControlWatsonProfile()[alignSectionOne.get(z)]);
						newAlignment.getTagProfile(a).setControlCrickProfile(antiZ, profileAlignment[a].getControlCrickProfile()[alignSectionOne.get(z)]);
					}
				}
				for (int b=0; b < alignmentTwo.getNumAligned(); b++){
					newAlignment.getTagProfile(numAligned+b).setCoord(antiZ, -1);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(numAligned+b).setWatsonProfile(c, antiZ, -1);
						newAlignment.getTagProfile(numAligned+b).setCrickProfile(c, antiZ, -1);
					}
					if (hasControl){
						newAlignment.getTagProfile(numAligned+b).setControlWatsonProfile(antiZ,-1);
						newAlignment.getTagProfile(numAligned+b).setControlCrickProfile(antiZ, -1);
					}
				}
			}else if(alignSectionOne.get(z)==last1 || alignSectionOne.get(z)==-1){
				// Gap in alignment 1; add in alignment 2's column  only
				for (int a=0; a < numAligned; a++){
					newAlignment.getTagProfile(a).setCoord(antiZ, -1);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(a).setWatsonProfile(c, antiZ, -1);
						newAlignment.getTagProfile(a).setCrickProfile(c, antiZ, -1);
					}
					if (hasControl){
						newAlignment.getTagProfile(a).setControlWatsonProfile(antiZ, -1);
						newAlignment.getTagProfile(a).setControlCrickProfile(antiZ, -1);
					}
				}
				for (int b=0; b < alignmentTwo.getNumAligned(); b++){
					if (nw.isReverse())
						coord = alignmentTwo.getTagProfile(b).getCoords()[alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];
					else
						coord = alignmentTwo.getTagProfile(b).getCoords()[alignSectionTwo.get(z)];
					newAlignment.getTagProfile(numAligned+b).setCoord(antiZ, coord);					
					for (int c=0; c < numCond; c++){
						if (nw.isReverse()){	// reverse direction
							currWatsonVal=alignmentTwo.getTagProfile(b).getCrickProfile()[c][alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];
							currCrickVal=alignmentTwo.getTagProfile(b).getWatsonProfile()[c][alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];											
						}else{	// forward direction
							currWatsonVal=alignmentTwo.getTagProfile(b).getWatsonProfile()[c][alignSectionTwo.get(z)];
							currCrickVal=alignmentTwo.getTagProfile(b).getCrickProfile()[c][alignSectionTwo.get(z)];
						}
						newAlignment.getTagProfile(numAligned+b).setWatsonProfile(c, antiZ, currWatsonVal);
						newAlignment.getTagProfile(numAligned+b).setCrickProfile(c, antiZ, currCrickVal);
					}
					// control experiment
					if (hasControl){
						if (nw.isReverse()){
							currControlWatsonVal=alignmentTwo.getTagProfile(b).getControlCrickProfile()[alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];
							currControlCrickVal=alignmentTwo.getTagProfile(b).getControlWatsonProfile()[alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];
						}else{
							currControlWatsonVal=alignmentTwo.getTagProfile(b).getControlWatsonProfile()[alignSectionTwo.get(z)];
							currControlCrickVal=alignmentTwo.getTagProfile(b).getControlCrickProfile()[alignSectionTwo.get(z)];				
						}
						newAlignment.getTagProfile(numAligned+b).setControlWatsonProfile(antiZ, currControlWatsonVal);
						newAlignment.getTagProfile(numAligned+b).setControlCrickProfile(antiZ, currControlCrickVal);
					}
				}
			}else{
				// No gap; add in both alignments
				for (int a=0; a < numAligned; a++){
					newAlignment.getTagProfile(a).setCoord(antiZ, profileAlignment[a].getCoords()[alignSectionOne.get(z)]);
					for (int c=0; c < numCond; c++){
						newAlignment.getTagProfile(a).setWatsonProfile(c, antiZ, profileAlignment[a].getWatsonProfile()[c][alignSectionOne.get(z)]);
						newAlignment.getTagProfile(a).setCrickProfile(c, antiZ, profileAlignment[a].getCrickProfile()[c][alignSectionOne.get(z)]);
					}
					if (hasControl){
						newAlignment.getTagProfile(a).setControlWatsonProfile(antiZ, profileAlignment[a].getControlWatsonProfile()[alignSectionOne.get(z)]);
						newAlignment.getTagProfile(a).setControlCrickProfile(antiZ, profileAlignment[a].getControlCrickProfile()[alignSectionOne.get(z)]);
					}
				}
				for (int b=0; b < alignmentTwo.getNumAligned(); b++){
					if (nw.isReverse())
						coord = alignmentTwo.getTagProfile(b).getCoords()[alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];
					else
						coord = alignmentTwo.getTagProfile(b).getCoords()[alignSectionTwo.get(z)];
					newAlignment.getTagProfile(numAligned+b).setCoord(antiZ, coord);				
					for (int c=0; c < numCond; c++){
						if (nw.isReverse()){	// reverse direction
							currWatsonVal=alignmentTwo.getTagProfile(b).getCrickProfile()[c][alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];
							currCrickVal=alignmentTwo.getTagProfile(b).getWatsonProfile()[c][alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];											
						}else{	// forward direction
							currWatsonVal=alignmentTwo.getTagProfile(b).getWatsonProfile()[c][alignSectionTwo.get(z)];
							currCrickVal=alignmentTwo.getTagProfile(b).getCrickProfile()[c][alignSectionTwo.get(z)];
						}
						newAlignment.getTagProfile(numAligned+b).setWatsonProfile(c, antiZ, currWatsonVal);
						newAlignment.getTagProfile(numAligned+b).setCrickProfile(c, antiZ, currCrickVal);
					}
					// control experiment
					if (hasControl){
						if (nw.isReverse()){
							currControlWatsonVal=alignmentTwo.getTagProfile(b).getControlCrickProfile()[alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];
							currControlCrickVal=alignmentTwo.getTagProfile(b).getControlWatsonProfile()[alignmentTwo.getTagProfile(b).getAlignLength()-alignSectionTwo.get(z)-1];
						}else{
							currControlWatsonVal=alignmentTwo.getTagProfile(b).getControlWatsonProfile()[alignSectionTwo.get(z)];
							currControlCrickVal=alignmentTwo.getTagProfile(b).getControlCrickProfile()[alignSectionTwo.get(z)];				
						}
						newAlignment.getTagProfile(numAligned+b).setControlWatsonProfile(antiZ, currControlWatsonVal);
						newAlignment.getTagProfile(numAligned+b).setControlCrickProfile(antiZ, currControlCrickVal);
					}
				}
			}
			last1=alignSectionOne.get(z);
			last2=alignSectionTwo.get(z);
			antiZ++;
		}
		
		return newAlignment;
	}
	
	public Pair<double[][][], double[][][]> makePerPointCountsFromTagProfiles(){		
		double[][][] watsonProfile = new double[numAligned][numCond][alignL];
		double[][][] crickProfile = new double[numAligned][numCond][alignL];
		for (int x=0; x < numAligned; x++){
			for (int z=0; z < alignL; z++){
				if (profileAlignment[x].getWatsonProfile()[0][z]==-1){
					for (int c=0; c < numCond ; c++){
						watsonProfile[x][c][z]=0.0;
						crickProfile[x][c][z]=0.0;
					}
				}else{
					for (int c=0; c < numCond ; c++){
						watsonProfile[x][c][z]=profileAlignment[x].getWatsonProfile()[c][z];
						crickProfile[x][c][z]=profileAlignment[x].getCrickProfile()[c][z];
					}		
				}}}			
		
		return new Pair<double[][][], double[][][]>(watsonProfile,crickProfile);
	}
	
	public Pair<double[], double[]> makeControlCompositeCounts(){
		double[] ctrlWatsonComposite=null;
		double[] ctrlCrickComposite=null;
		if (hasControl){
			ctrlWatsonComposite = new double[alignL];
			ctrlCrickComposite = new double[alignL];
			for (int z=0; z < alignL; z++){
				ctrlWatsonComposite[z]=0;
				ctrlCrickComposite[z]=0;
				for (int x=0; x < numAligned; x++){
					if (profileAlignment[x].getWatsonProfile()[0][z]!=-1){
						ctrlWatsonComposite[z]+=profileAlignment[x].getControlWatsonProfile()[z];
						ctrlCrickComposite[z]+=profileAlignment[x].getControlCrickProfile()[z];
					}}}	
		}
		
		return new Pair<double[], double[]>(ctrlWatsonComposite, ctrlCrickComposite);
	}
	
	/**
	 * Testing only
	 */
	public void printOriginalRegionsToFile(String prefix, int win, boolean sort){
		//Sort by ids
		if (sort)
			Arrays.sort(profileAlignment);
		try {
			FileWriter fout = new FileWriter(prefix+"_original.regions");
			for (TagProfile profile : profileAlignment){
				StrandedPoint pt = profile.getStrandedPoint();
				fout.write("chr"+pt.getChrom()+":");
				if(pt.getStrand()=='+')
					for (int w=0; w < win; w++)
						fout.write(pt.getLocation()-win/2+w+",");
				else
					for (int w=0; w < win; w++)
						fout.write(pt.getLocation()+win/2-w-1+",");
				fout.write("\n");
			}fout.close();
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void printAlignedRegionsToFile(String prefix, boolean sort){
		//Sort by ids
		if (sort)
			Arrays.sort(profileAlignment);
		try {
			FileWriter fout = new FileWriter(prefix+"_aligned.regions");
			for (TagProfile profile : profileAlignment){
				fout.write("chr"+profile.getStrandedPoint().getChrom()+":");
				for (int w=0; w < alignL; w++)
					fout.write(profile.getCoords()[w]+",");
				fout.write("\n");
			}fout.close();
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void printAlignmentCenterToFile(String prefix, boolean sort){
		//sort by ids
		if (sort)
			Arrays.sort(profileAlignment);
		try {
			FileWriter fout = new FileWriter(prefix+"_alignmentcenter.bed");
			for (TagProfile profile: profileAlignment){
				String chr=profile.getStrandedPoint().getChrom();
				char strand ;
				int centercoord =1; int start=-1;int end=-1;
				// search for start
				for (int w=0; w < alignL; w++){
					if (profile.getCoords()[w]>-1){
						start= profile.getCoords()[w];
						break;
					}}	
				//search for end
				for (int w=alignL-1; w >=0; w--){
					if (profile.getCoords()[w]>-1){
						end= profile.getCoords()[w];
						break;
					}}
				if (start <end){strand='+';}
				else{strand='-';}
				//serach for center
				for (int w=0; w<alignL/2;w++){
					if (profile.getCoords()[alignL/2+w]>-1){
						centercoord=profile.getCoords()[alignL/2+w];
						break;
					}else if (profile.getCoords()[alignL/2-w]>-1){
						centercoord=profile.getCoords()[alignL/2-w];
						break;
					}
				}
				String outstring="chr"+chr+"\t"+centercoord+"\t"+	(centercoord+1)+"\t.\t.\t"+strand+"\n";
				fout.write(outstring);
				
			}fout.close();
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void printOriginalTagsToFile(ExperimentManager manager, CompositeTagDistribution maker, String prefix, boolean sort){
		//Sort by ids
		if (sort)
			Arrays.sort(profileAlignment);
		for (ExperimentCondition cond : manager.getConditions()){
			try {
				FileWriter fout = new FileWriter(prefix+"_original."+cond.getName()+".mat");
				for (int i=0; i < maker.getPoints().size(); i++){
					double[] watson = maker.getPointWatson(maker.getPoints().get(i), cond);
					double[] crick = maker.getPointCrick(maker.getPoints().get(i), cond);
					fout.write("index:"+i+",watson\t");
					for (int w=0; w < watson.length; w++)
						fout.write(watson[w]+",");
					fout.write("\nindex:"+i+",crick\t");
					for (int w=0; w < crick.length; w++)
						fout.write(crick[w]+",");
					fout.write("\n");	
				}fout.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}}
	}
	
	public void printAlignedTagsToFile(ExperimentManager manager, String prefix, boolean sort){			
		//Sort by ids
		if (sort)
			Arrays.sort(profileAlignment);
		for (ExperimentCondition cond: manager.getConditions()){
			try {
				FileWriter fout = new FileWriter(prefix+"_aligned."+cond.getName()+".mat");
				for (TagProfile profile : profileAlignment){
					double[] watson = profile.getWatsonProfile()[cond.getIndex()];
					double[] crick = profile.getCrickProfile()[cond.getIndex()];
					fout.write("index:"+profile.getID()+",watson\t");
					for (int w=0; w < alignL; w++)
						fout.write(watson[w]+",");			
					fout.write("\nindex:"+profile.getID()+",crick\t");
					for (int w=0; w < alignL; w++)
						fout.write(crick[w]+",");
					fout.write("\n");
				}fout.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}}
	}
	
	/**
     * Print composite plot to file
     * @param a set of alignment
     */
	public void printAlignedCompositeToFile(ExperimentManager manager, String prefix){	
		for (ExperimentCondition cond : manager.getConditions()){
			try {
				FileWriter fout = new FileWriter(prefix+"_aligned_composite."+cond.getName()+".txt");
				double[] watsonComposite = new double[alignL];
				double[] crickComposite = new double[alignL];
				for (int w=0; w < alignL; w++){
					watsonComposite[w]=0;
					crickComposite[w]=0;
				}
				for (int i=0; i < numAligned; i++){
					for (int w=0; w < alignL; w++){
						if (profileAlignment[i].getWatsonProfile()[cond.getIndex()][w]>=0){
							watsonComposite[w]+=profileAlignment[i].getWatsonProfile()[cond.getIndex()][w];
							crickComposite[w]+=profileAlignment[i].getCrickProfile()[cond.getIndex()][w];
						}}}	
				for(int w=0; w<alignL; w++){
					int pos = (w-alignL/2);
					fout.write(pos+"\t"+watsonComposite[w]+"\t"+crickComposite[w]+"\n");
				}
				fout.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public static List<Integer> removeZeroCounts(ExperimentManager manager, CompositeTagDistribution maker, boolean sort){						
		List<ArrayListSorter> filteredList = new ArrayList<ArrayListSorter>();
		//Iterate through points
		for (int p=0; p <maker.getPoints().size(); p++){
			StrandedPoint pt=maker.getPoints().get(p);
			float sum=0;
			ArrayListSorter elem = new ArrayListSorter();
			for(ExperimentCondition cond : manager.getConditions()){
				elem.setKey(p);
				double[] watsontags = maker.getPointWatson(pt,cond);
				double[] cricktags = maker.getPointCrick(pt, cond);
				for (int i=0; i<watsontags.length; i++)
					sum += (float) watsontags[i];
				for (int i=0; i<cricktags.length;i++)
					sum += (float) cricktags[i];
			}
			if (sum > 0){
				elem.setVal(sum);
				filteredList.add(elem);	
			}
		}
		// Sort per location tags based on occupancies
		if (sort)
			Collections.sort(filteredList);	
		
		List<Integer> filteredIndexes = new ArrayList<Integer>();
		for (ArrayListSorter sor : filteredList){
			System.out.println(sor.getKey().toString());
			filteredIndexes.add((Integer) sor.getKey());
		}
		
		return filteredIndexes;
	}
}
