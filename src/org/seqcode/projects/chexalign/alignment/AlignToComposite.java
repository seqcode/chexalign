package org.seqcode.projects.chexalign.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Args;
import org.seqcode.projects.chexalign.multicompositemodel.CompositeXLFinderMultiCond;
import org.seqcode.projects.chexalign.xlanalysis.XLAnalysisConfig;

public class AlignToComposite {
	protected Genome genome;
	protected AlignmentConfig config;
	protected XLAnalysisConfig xlconfig;
	protected List<StrandedPoint> spts;
	protected int numPoints;
	protected int numCond;
	protected int win;
	protected List<String> fileListName;
	protected TagProfile[] tagProfiles;	// initial tag profiles
	protected double gapScaling=0;
	protected ExperimentManager manager;
	protected CompositeTagDistribution signalComposite;
	protected CompositeTagDistribution controlComposite=null;
	protected MultipleAlignment alignmentRec;	//alignment records
	protected String filename; 
	protected List<Integer> anchorIndices = new ArrayList<Integer>();
	protected CompositeXLFinderMultiCond xlFinder;
	protected boolean normScore=true;
	
	public AlignToComposite(AlignmentConfig config, XLAnalysisConfig xlcon, Genome gen, ExptConfig econ, ExperimentManager eMan, List<String> flistName){
		this.config=config;
		genome=gen;
		xlconfig = xlcon;
		spts = config.getStrandedPoints();
		manager = eMan;
		win = config.getWindowSize();
		filename = config.getFilename();
		numPoints = spts.size();
		numCond = manager.getNumConditions();
		fileListName = flistName;
		gapScaling = config.getGapExtScalingFactor();
		signalComposite = new CompositeTagDistribution(spts, eMan, win, true);
		// check for control experiment
		if (manager.getReplicates().get(0).hasControl())
			controlComposite = new CompositeTagDistribution(spts, eMan, win, false);
		for(ExperimentCondition cond : eMan.getConditions()){
			String compositeFileName = this.filename+"_composite."+cond.getName()+".txt";
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
	}
	
	public TagProfile makeTagProfileFromCompositeFiles() throws IOException{
		// check number of lines from first file
		BufferedReader reader = new BufferedReader(new FileReader(fileListName.get(0)));
		int numLine=0; 
		String line;
		while ((line = reader.readLine())!=null)
			numLine++;

		double[][] watsonProfile = new double[numCond][numLine];
		double[][] crickProfile = new double[numCond][numLine];	
		for (String fname : fileListName){
			String base = new File(fname).getName();
			String factorName = base.split("\\.")[1];
			System.out.println(factorName);
			int condIndex=-1;
			for (ExperimentCondition cond : manager.getConditions()){
				System.out.println("cond name : "+cond.getName());
				if (factorName.equals(cond.getName()))
					condIndex = cond.getIndex();
			}
			if (condIndex==-1){
				System.err.println("condition name does not match");
				System.exit(-1);
			}
			
			reader = new BufferedReader(new FileReader(fname));
			int j=0;
			while ((line = reader.readLine())!=null){
				line = line.trim();
				String[] words = line.split("\t");
				watsonProfile[condIndex][j]=Double.valueOf(words[1]);
				crickProfile[condIndex][j]=Double.valueOf(words[2]);		
				j++;
			}		
		}
		
		TagProfile profile = new TagProfile(numLine, numCond, watsonProfile, crickProfile, null, null);	
		return profile;	
	}
	
	public void execute() throws IOException{
		TagProfile compositeProfile = makeTagProfileFromCompositeFiles();
		compositeProfile.setID(numPoints);
		compositeProfile.setInitialCoords(new StrandedPoint(genome, "1", 1000, '+')); // some location
		for (int i=0; i < numPoints; i++){
			NeedlemanWunschAffine dist = NeedlemanWunschAffine.AlignPair(manager, config, compositeProfile, tagProfiles[i], config.getGapPenalty(), gapScaling);
			alignmentRec = new MultipleAlignment(1,dist.getAlignedLength(), numCond, controlComposite==null? false : true);
			alignmentRec.setTagProfile(0, compositeProfile);
			alignmentRec = alignmentRec.SingleProfileAddition(dist, tagProfiles[i]);
		}
		
		// print alignment results
		alignmentRec.printOriginalRegionsToFile(filename, win, config.useSortForPrint());  // this is making errors
		alignmentRec.printAlignedRegionsToFile(filename, config.useSortForPrint());
						
		alignmentRec.printOriginalTagsToFile(manager, signalComposite, filename, config.useSortForPrint());
		alignmentRec.printAlignedTagsToFile(manager, filename, config.useSortForPrint());
						
		alignmentRec.printAlignedCompositeToFile(manager, filename);	
	}
	
	
	/**
	 * Main driver method
	 * @param args
	 * @throws IOException 
	 * @throws Exception 
	 */
	public static void main(String[] args) throws IOException{
		
		GenomeConfig gcon = new GenomeConfig(args);
		Genome gen = gcon.getGenome();
		String fileName = Args.parseString(args, "complist", null);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);						
		AlignmentConfig config = new AlignmentConfig(gcon, args);			
		if (!config.useReadFilter())
			econ.setPerBaseReadFiltering(false);
		econ.setLoadRead2(false);
			
		if(gen==null || fileName==null){printError();}
			
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		ArrayList<String> fileList = new ArrayList<String>();
		String line;
		while ((line = reader.readLine())!=null){
			fileList.add(line.trim());		
		}

			
		ExperimentManager manager = new ExperimentManager(econ);
			
		XLAnalysisConfig ccon = new XLAnalysisConfig(gcon, args);					
		AlignToComposite alignment = new AlignToComposite(config,ccon, gen, econ, manager, fileList);	
		alignment.execute();
	}
	
	private static void printError(){
		System.err.println("Usage: AlignToComposite --species <organism;genome> --seq <sequence path name>\n" +
				"--complist <list of composite file name>\n"+
				"--motif <motif names> --mback <background model name> --mthres <threshold>\n" +
				"--posf <coordinate file name> --out <output root name> \n" +
				"--color <red/green/blue> \n");
		System.exit(1);
	}

}
