package org.seqcode.projects.chexalign.xlanalysis;

import java.io.File;
import java.util.List;

import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.StrandedPoint;


public class XLScanner {

	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected XLAnalysisConfig cconfig;
	protected ProteinDNAInteractionModel model;
	protected WindowedTagDistributions tagsDistribs;
	protected CompositeModelScan scanner;
	protected List<Point> scanPoints;
	protected int scanWindow=50;
	
	
	public XLScanner(GenomeConfig gcon, ExptConfig econ, XLAnalysisConfig ccon){
		gconfig = gcon;
		econfig = econ;
		cconfig = ccon;
		cconfig.makeXLAnalysisOutputDirs(true);
		manager = new ExperimentManager(econfig);
	}
	
	public void execute(){
		//Load appropriate data
		model = ProteinDNAInteractionModel.loadFromFile(cconfig, new File(cconfig.getModelFilename()));
		int winSize = model.getWidth();
		scanPoints = cconfig.getScanPoints();
		
		tagsDistribs = new WindowedTagDistributions(scanPoints, manager, winSize+scanWindow, true);
		
		System.err.println("Scanning "+scanPoints.size()+" windows of length "+scanWindow+"bp.");
		scanner = new CompositeModelScan(model, scanWindow, gconfig, econfig, cconfig, manager);
		List<StrandedPoint> scanMaxPoints = scanner.scan(tagsDistribs);
		
		//Report (TODO: replace with file output)
		for(StrandedPoint pt : scanMaxPoints){
			System.out.println(pt.toString());
		}
	}
	
	
	
	public static void main(String[] args){
		System.setProperty("java.awt.headless", "true");
		System.err.println("XLAnalysis version: "+XLAnalysisConfig.version);
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		//econ.setLoadRead2(false);//Enforce for chip-exo
		XLAnalysisConfig ccon = new XLAnalysisConfig(gcon, args);
		if(ccon.helpWanted()){
			System.err.println(ccon.getArgsList());
		}else if(ccon.getModelFilename()==null){
			System.err.println("Error: no ChExMix model provided. Use --model");
		}else{
			XLScanner xlScan = new XLScanner(gcon, econ, ccon);
			xlScan.execute();
		}
		
	}
}
