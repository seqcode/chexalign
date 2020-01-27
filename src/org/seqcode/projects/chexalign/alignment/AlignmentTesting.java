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

public class AlignmentTesting {
	protected AlignmentConfig config;
	protected List<StrandedPoint> spts;
	protected int numPoints;
	protected int numCond;
	protected ExperimentManager manager;
	protected CompositeTagDistribution maker;
	protected String filename;
	protected static final double MINIMUM_VALUE = -10000;
	protected static final double NUM_ITR=2;
	protected int alignmentCenterIndex=0;
	protected List<Integer> sortedPtIndexes;
	protected MultipleAlignment[] multiAlign;
	protected TagProfile[] testTagProfile;
	protected int windowsize;	// test data windows size;
	
	public AlignmentTesting(AlignmentConfig config, ExperimentManager eMan, double [][][] testData){
		this.config = config;
		spts = config.getStrandedPoints();
		manager = eMan;
		filename = config.getFilename();
		numPoints = testData.length;
		numCond = manager.getNumConditions();
		maker = new CompositeTagDistribution(spts, eMan, config.getWindowSize(), true);
		for(ExperimentCondition cond : eMan.getConditions()){
			String compositeFileName = filename+"_composite."+cond.getName()+".txt";
			maker.printProbsToFile(cond, compositeFileName);
		}	
		
		// make test data
		windowsize = testData[0][0].length;
		testTagProfile = new TagProfile[4];
		for (int i=0; i < 4; i++){
			double[][] wData = new double[1][windowsize];
			double[][] cData = new double[1][windowsize];
			wData[0] = testData[i][0];
			cData[0] = testData[i][1];
			testTagProfile[i]= new TagProfile(windowsize, 1, wData, cData, null, null);	
			testTagProfile[i].setID(i);
		}
		multiAlign = new MultipleAlignment[3];
	}
	
	public void iterativeAlignmentTesting(double gap_penalty, double gapScalingFactor){
		System.out.println("iterative alignment testing");
		
		List<Integer> orderedIndex = new ArrayList<Integer>();
		orderedIndex.add(0);
		NeedlemanWunschAffine dist = NeedlemanWunschAffine.AlignPair(manager, config, testTagProfile[3], testTagProfile[0], gap_penalty,gapScalingFactor);
		multiAlign[0] = new MultipleAlignment(1, dist.getAlignedLength(), numCond, false);		
		dist.setIndexB(1);
		multiAlign[0].setTagProfile(0, testTagProfile[3]);
		multiAlign[0] = multiAlign[0].SingleProfileAddition(dist, testTagProfile[0]);
		multiAlign[0].printAlignedTagsToFile(manager, filename+"-intermediate1", config.useSortForPrint());
		orderedIndex.add(1);
		
		NeedlemanWunschAffine dist2 = NeedlemanWunschAffine.AlignPair(manager, config, multiAlign[0].Alignment2Profile(), testTagProfile[2], gap_penalty, gapScalingFactor);		
		dist2.setIndexB(2);
		multiAlign[1] = multiAlign[0].SingleProfileAddition(dist2, testTagProfile[2]);
		multiAlign[1].printAlignedTagsToFile(manager, filename+"-intermediate2", config.useSortForPrint());
		orderedIndex.add(2);
		
		NeedlemanWunschAffine dist3 = NeedlemanWunschAffine.AlignPair(manager, config, multiAlign[1].Alignment2Profile(), testTagProfile[1], gap_penalty, gapScalingFactor);
		dist3.setIndexB(3);
		multiAlign[2] = multiAlign[1].SingleProfileAddition(dist3, testTagProfile[1]);
		orderedIndex.add(3);;
		
		multiAlign[2].printOriginalTagsToFile(manager, maker, filename,config.useSortForPrint());
		multiAlign[2].printAlignedTagsToFile(manager, filename, config.useSortForPrint());
	}
	
	public void multipleAlignmentTesting(double gap_penalty, double gapScalingFactor){
		System.out.println("multiple alignment testing");
		
		List<Integer> orderedIndex1 = new ArrayList<Integer>();
		List<Integer> orderedIndex2 = new ArrayList<Integer>();
		
		// First set of multiple alignment
		orderedIndex1.add(0);
		NeedlemanWunschAffine dist = NeedlemanWunschAffine.AlignPair(manager, config, testTagProfile[1], testTagProfile[3], gap_penalty, gapScalingFactor);		
		dist.setIndexB(1);
		multiAlign[0] = new MultipleAlignment(1, dist.getAlignedLength(), numCond, false);
		multiAlign[0].setTagProfile(0, testTagProfile[1]);
		multiAlign[0] = multiAlign[0].SingleProfileAddition(dist, testTagProfile[3]);
		orderedIndex1.add(1);
		
		multiAlign[0].printAlignedTagsToFile(manager, filename+"-intermediate1", config.useSortForPrint());
	
		// second set of multiple alignment
		orderedIndex2.add(2);
		NeedlemanWunschAffine dist2 = NeedlemanWunschAffine.AlignPair(manager, config, testTagProfile[2], testTagProfile[0], gap_penalty, gapScalingFactor);				
		dist2.setIndexB(3);
		multiAlign[1] = new MultipleAlignment(1, dist.getAlignedLength(), numCond, false);
		multiAlign[1].setTagProfile(0, testTagProfile[2]);
		multiAlign[1] = multiAlign[1].SingleProfileAddition(dist2, testTagProfile[0]);
		orderedIndex2.add(3);
		
		multiAlign[1].printAlignedTagsToFile(manager, filename+"-intermediate2", config.useSortForPrint());
		
		NeedlemanWunschAffine nw = NeedlemanWunschAffine.AlignPair(manager, config, multiAlign[0].Alignment2Profile(), multiAlign[1].Alignment2Profile(), gap_penalty, gapScalingFactor);
		multiAlign[2] = multiAlign[0].MultipleProfileAddition(nw, multiAlign[1]);
		
		orderedIndex1.addAll(orderedIndex2);	
		System.out.print(orderedIndex1.toString());
		
		multiAlign[2].printAlignedTagsToFile(manager, filename, config.useSortForPrint());
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
			
			// hard coded test data
			double[][][] testdata = new double[4][2][];				
			double[] t0_watson = { 26.0,32.0,16.0,58.0,52.0,37.0,48.0,30.0,38.0,57.0,32.0,30.0,18.0,26.0,15.0 };
			double[] t0_crick = { 1.0,2.0,10.0,2.0,14.0,12.0,5.0,3.0,9.0,22.0,26.0,13.0,10.0,9.0,23.0 };		
			// gapped example
			double[] t1_watson = { 26.0,32.0,16.0,0.0,0.0,58.0,52.0,37.0,48.0,30.0,38.0,57.0,32.0,30.0,18.0 };
			double[] t1_crick = { 1.0,2.0,10.0,0.0,0.0,2.0,14.0,12.0,5.0,3.0,9.0,22.0,26.0,13.0,10.0 };
			// gapped and reverse
			double[] t2_watson = { 10.0,13.0,0.0,0.0,0.0,26.0,22.0,9.0,3.0,5.0,12.0,14.0,2.0,10.0,2.0 };
			double[] t2_crick = {18.0,30.0,0.0,0.0,0.0,32.0,57.0,38.0,30.0,48.0,37.0,52.0,58.0,16.0,32.0 };
			// gapped and reverse (2)
			double[] t3_watson = { 9.0,10.0,13.0,0.0,26.0,22.0,9.0,3.0,5.0,12.0,14.0,2.0,10.0,2.0,1.0 };
			double[] t3_crick = {26.0,18.0,30.0,0.0,32.0,57.0,38.0,30.0,48.0,37.0,52.0,58.0,16.0,32.0,26.0 };
			
			testdata[0][0]=t0_watson; testdata[0][1]=t0_crick; 
			testdata[1][0]=t1_watson; testdata[1][1]=t1_crick; 
			testdata[2][0]=t2_watson; testdata[2][1]=t2_crick; 
			testdata[3][0]=t3_watson; testdata[3][1]=t3_crick; 
			
			
			ExperimentManager manager = new ExperimentManager(econ);
					
			AlignmentConfig config = new AlignmentConfig(gcon, args);
					
			AlignmentTesting test = new AlignmentTesting(config, manager, testdata);
				
			if (ap.hasKey("iterative")){
				// iterative alignment testing
				test.iterativeAlignmentTesting(config.getGapPenalty(), config.getGapExtScalingFactor());
			}else{
				//multiple alignment testing
				test.multipleAlignmentTesting(config.getGapPenalty(), config.getGapExtScalingFactor());
			}
					
			manager.close();
		}	
	}			
}
