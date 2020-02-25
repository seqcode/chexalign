package org.seqcode.projects.chexalign.alignment;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.projects.chexalign.ChExAlign;
import org.seqcode.projects.chexalign.xlanalysis.XLAnalysisConfig;


public class RandomRegionAlignment {
	
	protected AlignmentConfig config;
	protected XLAnalysisConfig xlconfig;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
	
	public RandomRegionAlignment(AlignmentConfig config, XLAnalysisConfig xlcon, GenomeConfig gcon,ExptConfig econ, ExperimentManager eMan){
		this.config = config;
		xlconfig = xlcon;
		gconfig = gcon;
		econfig = econ;
		manager = eMan;
	}
	
	protected List<StrandedPoint> randomPointPick(int numSamples){
		List<StrandedPoint> points = new ArrayList<StrandedPoint>();
		Random rand = new Random();
		
		//First see how big the genome is:
		int numChroms=0;
		long genomeSize=0;
		long [] chromoSize = new long[gconfig.getGenome().getChromList().size()];
		String [] chromoNames = new String[gconfig.getGenome().getChromList().size()];
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gconfig.getGenome());
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			genomeSize += (double)currentChrom.getWidth();
			chromoSize[numChroms]=currentChrom.getWidth();
			chromoNames[numChroms]=currentChrom.getChrom();
			numChroms++;				
		}

		//Now, iteratively generate random positions and check if they are valid and not overlapping repeats. 
		for (int i=0; i < numSamples; i++){
			StrandedPoint p;				
			long randPos = (long)(1+(rand.nextDouble()*genomeSize));
			//find the chr
			long total=0;
			for(int c=0; c<numChroms; c++){
				if(randPos<total+chromoSize[c]){
					int pos = (int)(randPos-total+1);				
					p = new StrandedPoint(gconfig.getGenome(), chromoNames[c], pos, '+');
					points.add(p);
				}total+=chromoSize[c];
			}
		}
		return(points);
	}
	
	public void run() throws IOException{
		
		List<StrandedPoint> points = randomPointPick(config.getStrandedPoints().size());
		config.setStrandedPoints(points);
		
		ChExAlign nodes = new ChExAlign(config, xlconfig, gconfig, econfig, manager);		
		nodes.buildTree();	
	}
	
	/**
	 * Main driver method for ChExAlign
	 * @param args
	 * @throws IOException 
	 * @throws Exception 
	 */
	public static void main(String[] args) throws IOException{
		System.setProperty("java.awt.headless", "true");
		System.out.println("ChExAlign version "+AlignmentConfig.version+"\n\n");
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);						
		AlignmentConfig config = new AlignmentConfig(gcon, args);
		
		// generate random regions
		
		
		if (!config.useReadFilter())
			econ.setPerBaseReadFiltering(false);
		econ.setLoadRead2(false);
		
		if(args.length==0 || config.helpWanted()){
			System.err.println(ChExAlign.getChExAlignArgsList());	
		}else{
			
			ExperimentManager manager = new ExperimentManager(econ);
			
			XLAnalysisConfig ccon = new XLAnalysisConfig(gcon, args);
			ccon.makeXLAnalysisOutputDirs(true);
			
			RandomRegionAlignment alignment = new RandomRegionAlignment(config,ccon, gcon, econ, manager);
			alignment.run();
			manager.close();

		}
	}

}
