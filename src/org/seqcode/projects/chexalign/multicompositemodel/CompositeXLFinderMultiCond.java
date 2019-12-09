package org.seqcode.projects.chexalign.multicompositemodel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.IOUtil;
import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.projects.chexalign.xlanalysis.XLAnalysisConfig;



public class CompositeXLFinderMultiCond {

	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected XLAnalysisConfig cconfig;
	protected CompositeModelMixtureMultiCond mixtureModel;
	protected CompositeTagDistribution signalComposite;
	protected CompositeTagDistribution controlComposite;
	protected List<StrandedPoint> compositePoints;
	
	public CompositeXLFinderMultiCond(GenomeConfig gcon, ExptConfig econ, XLAnalysisConfig ccon){
		gconfig = gcon;
		econfig = econ;
		cconfig = ccon;
		cconfig.makeXLAnalysisOutputDirs(true);
		manager = new ExperimentManager(econfig);
	}
	
	public CompositeXLFinderMultiCond(GenomeConfig gcon, ExptConfig econ, ExperimentManager manager, XLAnalysisConfig ccon){
		gconfig = gcon;
		econfig = econ;
		cconfig = ccon;
//		cconfig.makeXLAnalysisOutputDirs(true);
		this.manager = manager;
	}
	
	public void execute(){
		//Load appropriate options
		compositePoints = cconfig.getCompositePoints();
		int winSize = cconfig.getCompositeWinSize();
		
		//Build the composite distribution(s)
		signalComposite = new CompositeTagDistribution(compositePoints, manager, winSize, true);
		
		//Check for control experiment
		boolean hasControl=true;
		for (ControlledExperiment rep : manager.getReplicates()){
			if (!rep.hasControl()){
				hasControl=false;
				break;
			}
		}
		controlComposite =null;
		if (hasControl)
			controlComposite = new CompositeTagDistribution(compositePoints, manager, winSize, false);
		
		//Initialize the mixture model 
		mixtureModel = new CompositeModelMixtureMultiCond(signalComposite, controlComposite, gconfig, econfig, cconfig, manager);
		
		//Train EM
		System.err.println("EM training");
		mixtureModel.trainEM();
		
		
		//ML assignment
		System.err.println("ML assignment");
		mixtureModel.assignML(false);
		
		//Report
		for(ExperimentCondition cond : manager.getConditions()){
			String compositeFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_composite."+cond.getName()+".txt";
			signalComposite.printProbsToFile(cond, compositeFileName);
			
			String perSiteRespFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_site-component-ML."+cond.getName()+".txt";
			mixtureModel.printPerSiteComponentResponsibilitiesToFile(cond, perSiteRespFileName);
			
			//Print component responsibility profiles
			String compProfileFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
			+"_component-profile."+cond.getName()+".txt";
			mixtureModel.printComponentProfilesToFile(cond, compProfileFileName);
		}
		mixtureModel.saveCompositePlots();
		
		//Save the model
		for(ExperimentCondition cond : manager.getConditions()){
			String modelFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()+"."+cond.getName()+".chexmix";
			mixtureModel.getModel().saveToFile(modelFileName, cond.getIndex());
		}
	}
	
	/**
	 * Run XL Finder on composite profiles from alignment methods
	 * @param alignL
	 * @param perPointWatson
	 * @param perPointCrick
	 */
	public void executeOnAlignedCompoiste(int alignL, List<StrandedPoint> spts, double[][][] perPointWatson, double[][][] perPointCrick, double[] ctrlWatsonComposite, double[] ctrlCrickComposite){
		
		signalComposite = new AlignedCompositeTagDistribution();
		((AlignedCompositeTagDistribution) signalComposite).initializeClass(spts, manager, alignL, false);
		((AlignedCompositeTagDistribution) signalComposite).setPerPointWatson(perPointWatson, ctrlWatsonComposite);
		((AlignedCompositeTagDistribution) signalComposite).setPerPointCrick(perPointCrick, ctrlCrickComposite);
		((AlignedCompositeTagDistribution) signalComposite).normalizeComposite();
		controlComposite =null;
		
		//Initialize the mixture model 
		mixtureModel = new CompositeModelMixtureMultiCond(signalComposite, controlComposite, gconfig, econfig, cconfig, manager);
		
		//Train EM
		System.err.println("EM training");
		mixtureModel.trainEM();
		
		
		//ML assignment
		System.err.println("ML assignment");
		mixtureModel.assignML(false);
		
		//Report
		for(ExperimentCondition cond : manager.getConditions()){
			String compositeFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_composite."+cond.getName()+".txt";
			signalComposite.printProbsToFile(cond, compositeFileName);
			
			String perSiteRespFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_site-component-ML."+cond.getName()+".txt";
			mixtureModel.printPerSiteComponentResponsibilitiesToFile(cond, perSiteRespFileName);
			
			//Print component responsibility profiles
			String compProfileFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
			+"_component-profile."+cond.getName()+".txt";
			mixtureModel.printComponentProfilesToFile(cond, compProfileFileName);
		}
		
		//Save the model
		for(ExperimentCondition cond : manager.getConditions()){
			String modelFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()+"."+cond.getName()+".chexmix";
			mixtureModel.getModel().saveToFile(modelFileName, cond.getIndex());
		}
		
		mixtureModel.saveCompositePlots();	//this step seems to give problems to some datasets
	}
	
	// Read aligned profiles 
	public void excecuteOnAlignedCompositeFile(String[] lines){
		
		int alignL=-1;
		double[][][] perPointWatson = null; double[][][] perPointCrick = null;	
		// match with experiments specified in manager.
		for (ExperimentCondition cond : manager.getConditions()){
			// search for this condition name
			for (String fname : lines){
				if (fname.contains(cond.getName())){
					System.out.println("file found with condition name : "+cond.getName()+" in line "+fname);
					List<String> data = new ArrayList<String>();
					// open file and read from aligned composite plots
					File file = new File(fname);
					if(!file.isFile()){System.err.println("Invalid file name: "+fname);System.exit(1);}
			        BufferedReader reader;
					try {
						reader = new BufferedReader(new FileReader(file));
						String line;
				        while ((line = reader.readLine()) != null) {
				            line = line.trim();
				            String[] words = line.split("\\s+");
				            data.add(words[1]);
				        }
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					if (perPointWatson == null){
						alignL = data.get(0).split(",").length;
						perPointWatson = new double[data.size()][manager.getNumConditions()][alignL];
						perPointCrick = new double[data.size()][manager.getNumConditions()][alignL];
					}			
					int index=0;
					for (int row=0; row < data.size(); row++){
						String[] results = data.get(row).split(",");
						if (row%2==0){ // watson strand
							int c=0; 
							for (String item : results){
								perPointWatson[index][cond.getIndex()][c]=Double.parseDouble(item);
								c++;
							}
						}else{ // crick strand
							int c=0; 
							for (String item : results){
								perPointCrick[index][cond.getIndex()][c]=Double.parseDouble(item);
								c++;
							}
							index++;
						}}}}}			
		
		executeOnAlignedCompoiste(alignL, cconfig.getCompositePoints(), perPointWatson, perPointCrick, null, null);
	}
	
	//Main
	public static void main(String[] args){
		System.setProperty("java.awt.headless", "true");
		System.err.println("XL analysis version: "+XLAnalysisConfig.version);
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		econ.setPerBaseReadFiltering(false);	
		econ.setLoadRead2(false);//Enforce for chip-exo
		XLAnalysisConfig ccon = new XLAnalysisConfig(gcon, args);
		
		ArgParser ap = new ArgParser(args);
		
		if(ccon.helpWanted()){
			System.err.println(ccon.getArgsList());
		}else{
			CompositeXLFinderMultiCond xlFinder = new CompositeXLFinderMultiCond(gcon, econ, ccon);
			if (ap.hasKey("flist")){
				// perform xl analysis on aligned files
				String fname=null;
				if (ap.hasKey("flist")){
					fname=Args.parseString(args, "flist", null);
				}
				if (fname!=null){
					String[] lines= IOUtil.readFile2Array(fname);
					System.out.println("printing the file names in flist");
					for (int i=0; i <lines.length; i++)
						System.out.println(lines[i]);
					xlFinder.excecuteOnAlignedCompositeFile(lines);
				}else{
					System.err.println("'flist' is specified but no file is given.");
				}				
			}else{
				xlFinder.execute();
			}
		}
	}
}
