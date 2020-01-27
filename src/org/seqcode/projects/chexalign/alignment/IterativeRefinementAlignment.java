package org.seqcode.projects.chexalign.alignment;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Pair;
import org.seqcode.projects.chexalign.multicompositemodel.CompositeXLFinderMultiCond;
import org.seqcode.projects.chexalign.xlanalysis.XLAnalysisConfig;

/**
 * IterativeAlignment implements iterative refinement algorithm from Barton-Sternberg multiple alignment.
 * 
 * i) Find the two sequences with the highest pairwise similarity and align them using standard pairwise dynamic programming alignment.
 * ii) Find the sequence that is most similar to a profile of the alignment of the first two, and align it to the first two by 
 *     profile-sequence alignment. Repeat until all sequences have been included in the multiple alignment.
 * iii) Remove sequence x1 and realign it to a profile of the other aligned sequences x2,...,xN by profile-sequence alignment. 
 * 		Repeat for sequences x2,...,xN.
 * iv) Repeat the previous realignment step a fixed number of times, or until the alignment score converges.
 * 
 * Depends on NeedlemanWunschAffine class.
 * 
 * @author naomi
 */
public class IterativeRefinementAlignment {	
	protected AlignmentConfig config;
	protected List<StrandedPoint> spts;
	protected int numPoints;
	protected int numCond;
	protected int win;
	protected ExperimentManager manager;
	protected CompositeTagDistribution signalComposite;
	protected CompositeTagDistribution controlComposite;
	protected MultipleAlignment alignmentRec;	//alignment records
	protected TagProfile[] tagProfiles;
	protected String filename;
	protected List<Integer> sortedPtIndexes;
	protected double gapScaling=0;
	protected CompositeXLFinderMultiCond xlFinder;
	
	public IterativeRefinementAlignment(AlignmentConfig config, ExperimentManager eMan, GenomeConfig gcon, ExptConfig econ, XLAnalysisConfig xlcon){
		this.config = config;
		spts = config.getStrandedPoints();
		manager = eMan;
		win = config.getWindowSize();
		filename = config.getFilename();
		numPoints = spts.size();
		numCond = manager.getNumConditions();
		signalComposite = new CompositeTagDistribution(spts, eMan, win, true);
		// check for control experiment
		if (manager.getReplicates().get(0).hasControl())
			controlComposite = new CompositeTagDistribution(spts, eMan, win, false);
		gapScaling = config.getGapExtScalingFactor();
		for(ExperimentCondition cond : eMan.getConditions()){
			String compositeFileName = filename+"_composite."+cond.getName()+".txt";
			signalComposite.printProbsToFile(cond, compositeFileName);
		}		
		
		// Make all initial tag profiles
		tagProfiles = new TagProfile[numPoints];
		for (int i=0; i < numPoints; i++){
			tagProfiles[i]= new TagProfile(win, numCond, signalComposite.getPointWatsons(i), signalComposite.getPointCricks(i),
					controlComposite==null? null: controlComposite.getPointWatsons(i)[0], controlComposite==null? null: controlComposite.getPointCricks(i)[0]);
			tagProfiles[i].setID(i);
			tagProfiles[i].setInitialCoords(spts.get(i));
		}
		
		xlFinder = new CompositeXLFinderMultiCond(gcon, econ, manager, xlcon);
	}
	
	
	public void alignByTag2Profile(){	
		// Remove regions with zero counts
		List<Integer> nonZerosReg = MultipleAlignment.removeZeroCounts(manager, signalComposite, false);
		
		// Keep track of which tags have been aligned
		List<Integer> orderedIndex = new ArrayList<Integer>();

		// 1) Determine tag pairs to align first
		double maxScore=config.MINIMUM_VALUE;
		NeedlemanWunschAffine currBestMatrix=null;
		int indexA=-1; int indexB=-1; 
		for (int i=0; i < nonZerosReg.size(); i++){
			// Convert i to a profile
			for (int j=i+1; j < nonZerosReg.size(); j++){
				NeedlemanWunschAffine dist = NeedlemanWunschAffine.AlignPair(manager, config, tagProfiles[i], tagProfiles[j], config.getGapPenalty(), gapScaling);
				double currScore = dist.getMaxScore();
				if (currScore > maxScore){
					maxScore=currScore;
					currBestMatrix=dist;
					indexA=i; indexB=j;
					}}}		
		alignmentRec = new MultipleAlignment(1,currBestMatrix.getAlignedLength(), numCond, controlComposite==null? false : true);
		alignmentRec.setTagProfile(0, tagProfiles[indexA]);
		alignmentRec = alignmentRec.SingleProfileAddition(currBestMatrix, tagProfiles[indexB]);
		orderedIndex.add(indexA);
		orderedIndex.add(indexB);
		
		//2) Determine the best profile to tag alignment
		while(orderedIndex.size() < nonZerosReg.size()){
			maxScore=config.MINIMUM_VALUE;
			NeedlemanWunschAffine bestMat=null;
			TagProfile currProfile = alignmentRec.Alignment2Profile();
			int bestJ=-1;
			for (int j=0; j < nonZerosReg.size(); j++){
				if (!orderedIndex.contains(nonZerosReg.get(j))){
					NeedlemanWunschAffine mat = NeedlemanWunschAffine.AlignPair(manager, config, currProfile, tagProfiles[j], config.getGapPenalty(), gapScaling);
					double currScore = mat.getMaxScore();
					if (currScore > maxScore){
						maxScore=currScore;
						bestMat=mat;
						bestJ=j;
						}}}	
			alignmentRec = alignmentRec.SingleProfileAddition(bestMat, tagProfiles[bestJ]);
			orderedIndex.add(bestJ);
		}		
		
//		multiAlign.printAlignedPointsToFile(filename+"c_0_aligned.points", orderedIndex);	

		//3) Remove element and re-align
		iterativeAlignment(orderedIndex);
	}
	
	public void alignHighOccupanciesFirst(){	
		List<Integer> sortedPtIndexes = MultipleAlignment.removeZeroCounts(manager, signalComposite, true);
		
		// First alignment
		NeedlemanWunschAffine firstdist = NeedlemanWunschAffine.AlignPair(manager, config, 
				tagProfiles[sortedPtIndexes.get(0)], tagProfiles[sortedPtIndexes.get(1)], config.getGapPenalty(), gapScaling);
		// Only make alignment profile for A, B will be made within SingleProfileAddition()
		alignmentRec = new MultipleAlignment(1,firstdist.getAlignedLength(), numCond, controlComposite==null? false : true);
		alignmentRec.setTagProfile(0, tagProfiles[sortedPtIndexes.get(0)]);
		alignmentRec = alignmentRec.SingleProfileAddition(firstdist, tagProfiles[sortedPtIndexes.get(1)]);
		// Align the rest
		for (int j=2; j < sortedPtIndexes.size(); j++){
			int currIndex = sortedPtIndexes.get(j);
			NeedlemanWunschAffine mat = NeedlemanWunschAffine.AlignPair(manager, config, alignmentRec.Alignment2Profile(), tagProfiles[currIndex], config.getGapPenalty(), gapScaling);
			alignmentRec = alignmentRec.SingleProfileAddition(mat, tagProfiles[currIndex]);
		}
		
//		multiAlign.printAlignedPointsToFile(filename+"c_0_aligned.points", sortedPtIndexes);
		
		//Iterate
		iterativeAlignment(sortedPtIndexes);
		
		// print alignment results
		alignmentRec.printOriginalRegionsToFile(filename, win, config.useSortForPrint());
		alignmentRec.printAlignedRegionsToFile(filename, config.useSortForPrint());
				
		alignmentRec.printOriginalTagsToFile(manager, signalComposite, filename, config.useSortForPrint());
		alignmentRec.printAlignedTagsToFile(manager, filename, config.useSortForPrint());
				
		alignmentRec.printAlignedCompositeToFile(manager, filename);
				
		if (config.doXLAnalysis()){
			//perform cross-linking analysis
			Pair<double[][][], double[][][]> perPointTags = alignmentRec.makePerPointCountsFromTagProfiles();
			Pair<double[], double[]> ctrlCompositeTags = alignmentRec.makeControlCompositeCounts();
			xlFinder.executeOnAlignedCompoiste(alignmentRec.getAlignedLength(), spts, perPointTags.car(), perPointTags.cdr(),ctrlCompositeTags.car(), ctrlCompositeTags.cdr());
		}
	}
	
	public void iterativeAlignment(List<Integer> orderedPtIndexes){
		for (int numItr=0; numItr < config.NUM_ITR; numItr++){			
			for (int j=0; j < numPoints; j++){
				int removeID = orderedPtIndexes.remove(0);
				alignmentRec.SingleProfileSubstraction(removeID);
				NeedlemanWunschAffine mat = NeedlemanWunschAffine.AlignPair(manager, config, alignmentRec.Alignment2Profile(), tagProfiles[removeID], config.getGapPenalty(), gapScaling);
				alignmentRec = alignmentRec.SingleProfileAddition(mat, tagProfiles[removeID]);
				orderedPtIndexes.add(removeID);
			}	
					
//			multiAlign.printAlignedPointsToFile(filename+"_c_"+(numItr+1)+"_aligned.points", orderedPtIndexes);
			
//			multiAlign.printAlignedTagsToFile(filename+"_c_"+(numItr+1),orderedPtIndexes);
		}				
	}
	
	
	// Main for testing
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);	
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		econ.setPerBaseReadFiltering(false);
		econ.setLoadRead2(false);
		if(args.length==0){
			System.err.println("IterativeRefinementAlignment:"+
					"\t--spoints <stranded point file> OR --points <point file>"+
					"\t--cwin <window around points>"+
					"\t--out <output file name>"+
					"Genome:" +
					"\t--species <Species;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"Experiment Design File:\n" +
					"\t--design <file name>\n");			
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
				
			AlignmentConfig config = new AlignmentConfig(gcon, args);
			
			XLAnalysisConfig xlcon = new XLAnalysisConfig(gcon, args);
			
			IterativeRefinementAlignment align = new IterativeRefinementAlignment(config, manager, gcon, econ, xlcon);

			if (ap.hasKey("hightagfirst"))
				align.alignHighOccupanciesFirst();
			else
				align.alignByTag2Profile(); // Align from the most similar pair
				
			manager.close();
		}	
	}
}
