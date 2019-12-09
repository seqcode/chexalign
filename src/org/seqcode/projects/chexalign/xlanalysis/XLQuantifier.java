package org.seqcode.projects.chexalign.xlanalysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;


public class XLQuantifier {

	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected XLAnalysisConfig cconfig;
	protected ProteinDNAInteractionModel model;
	protected CompositeModelMixture mixtureModel;
	protected CompositeModelEM EMtrainer;
	protected CompositeTagDistribution signalComposite;
	protected CompositeTagDistribution controlComposite;
	protected List<StrandedPoint> compositePoints;
	private final char[] LETTERS = {'A','C','G','T'};
	
	
	public XLQuantifier(GenomeConfig gcon, ExptConfig econ, XLAnalysisConfig ccon){
		gconfig = gcon;
		econfig = econ;
		cconfig = ccon;
		cconfig.makeXLAnalysisOutputDirs(true);
		manager = new ExperimentManager(econfig);
	}
	
	public void execute(){
		//Load appropriate data
		model = ProteinDNAInteractionModel.loadFromFile(cconfig, new File(cconfig.getModelFilename()));
		compositePoints = cconfig.getCompositePoints();
		int winSize = cconfig.getCompositeWinSize();
		
		//Build the composite distribution(s)
		signalComposite = new CompositeTagDistribution(compositePoints, manager, winSize, true);
		//controlComposite = new CompositeTagDistribution(points, manager, winSize, false);
		controlComposite =null;
		
		//Initialize the mixture model 
		mixtureModel = new CompositeModelMixture(signalComposite, controlComposite, gconfig, econfig, cconfig, manager);
		
		//Set the loaded model
		mixtureModel.setModel(model);
		
		//Set composite-level responsibilities in model
		EMtrainer = new CompositeModelEM(signalComposite, cconfig, manager);
		try {
			EMtrainer.train(model, 0, false);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		//ML assignment at per-site level
		System.err.println("ML assignment");
		mixtureModel.assignML(false);
		//Get the per-site assignments
		List<CompositeModelSiteAssignment> assignments = mixtureModel.getSiteAssignments();
		Map<StrandedPoint, String> pointSeqs = retrieveSequences(assignments, signalComposite.getWinSize());
		 
		//Build per-XL-point weighted PWMs
		Map<CompositeModelComponent, WeightMatrix> componentMotifs = getPerComponentWeightMatrices(assignments, pointSeqs, 20, false);
		Map<CompositeModelComponent, WeightMatrix> componentMotifsMaxPeaksOnly = getPerComponentWeightMatrices(assignments, pointSeqs, 20, true);
		WeightMatrix allSiteMotif = getAllPointWeightMatrix(pointSeqs, 100);
		
		//Estimate corrected XL tag counts
		for(ExperimentCondition cond : manager.getConditions()){
			String correctionFilename = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_XL-bias-correction."+cond.getName()+".txt";
			correctXLBias(cond, assignments, new File(correctionFilename));
		}
		
		//Print XL bias factors
		for(ExperimentCondition cond : manager.getConditions()){
			String biasFilename = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_XL-bias-factors."+cond.getName()+".txt";
			printXLBias(cond, assignments, new File(biasFilename));
		}
		
		//Print SW base indicators
		String swFilename = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()+"_SWbases.txt";
		printSWBases(new File(swFilename), pointSeqs, signalComposite.getWinSize());
		
		//Report
		for(ExperimentCondition cond : manager.getConditions()){
			String compositeFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_composite."+cond.getName()+".txt";
			signalComposite.printProbsToFile(cond, compositeFileName);
			
			String perSiteRespFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_site-component-ML."+cond.getName()+".txt";
			mixtureModel.printPerSiteComponentResponsibilitiesToFile(cond, perSiteRespFileName);
		}
		mixtureModel.saveCompositePlots();
		//Print motif logos
		for(CompositeModelComponent comp : componentMotifs.keySet()){
			String motifFileName = cconfig.getOutputImagesDir()+File.separator+cconfig.getOutBase()
					+"_xl-bias-motif."+(comp.getPosition()-signalComposite.getCenterOffset())+".png";
			WeightMatrix motif = componentMotifs.get(comp);
			String motifLabel = cconfig.getOutBase()+" "+(comp.getPosition()-signalComposite.getCenterOffset())+" XL bias motif";
			org.seqcode.motifs.DrawMotifs.printMotifLogo(motif, new File(motifFileName), 150, motifLabel, true);
			String motifFileNameMaxOnly = cconfig.getOutputImagesDir()+File.separator+cconfig.getOutBase()
					+"_xl-bias-motif-max-only."+(comp.getPosition()-signalComposite.getCenterOffset())+".png";
			motif = componentMotifsMaxPeaksOnly.get(comp);
			motifLabel = cconfig.getOutBase()+" "+(comp.getPosition()-signalComposite.getCenterOffset())+" XL bias motif (max XL only)";
			org.seqcode.motifs.DrawMotifs.printMotifLogo(motif, new File(motifFileNameMaxOnly), 150, motifLabel, true);
		}
		String allMotifFileName = cconfig.getOutputImagesDir()+File.separator+cconfig.getOutBase()
				+"_all-site-motif.png";
		String motifLabel = cconfig.getOutBase()+" all site motif";
		org.seqcode.motifs.DrawMotifs.printMotifLogo(allSiteMotif, new File(allMotifFileName), 710, motifLabel, true);
		
	}
	
	/**
	 * Get sequences corresponding to the stranded points in the site assignments
	 * @param siteAssignments
	 * @param seqWin
	 * @return
	 */
	protected Map<StrandedPoint, String> retrieveSequences(List<CompositeModelSiteAssignment> siteAssignments, int seqWin){
		Map<StrandedPoint, String> pointSequences = new HashMap<StrandedPoint, String>();
		SequenceGenerator<Region> seqgen = gconfig.getSequenceGenerator();
		
		for(CompositeModelSiteAssignment sa : siteAssignments){
			Region reg = sa.getPoint().expand(seqWin/2);
			String seq = seqgen.execute(reg);
			if(sa.getPoint().getStrand()=='-')
				seq = SequenceUtils.reverseComplement(seq);
			pointSequences.put(sa.getPoint(), seq);
		}
		return pointSequences;
	}
	
	
	/**
	 * Calculate weighted motifs for each non-zero XL component
	 * @return
	 */
	protected WeightMatrix getAllPointWeightMatrix(Map<StrandedPoint, String> pointSeqs, int motifWidth){
		int motifStartPos = signalComposite.getCenterOffset() - (motifWidth/2)+1;
		float[][] freq =  new float[motifWidth][WeightMatrix.MAXLETTERVAL];
		float[][] pwm = new float[motifWidth][WeightMatrix.MAXLETTERVAL];
		float[] sum =  new float[motifWidth];
		float pseudoCount=(float) 0.01;
		for(int x=0; x<motifWidth; x++){ 
			sum[x]=0; 
			for (int b=0;b<LETTERS.length;b++){
				freq[x][LETTERS[b]]=0;
		}}
		
		for(StrandedPoint pt : pointSeqs.keySet()){
			String currSeq = pointSeqs.get(pt).toUpperCase();
			for(int x=0; x<motifWidth; x++){ 
				char letter = currSeq.charAt(motifStartPos+x);
				freq[x][letter]++;
				sum[x]++;
			}
		}
		//Normalize freq matrix
		//Log-odds over expected motif
		for(int x=0; x<motifWidth; x++){ 
			for (int b=0;b<LETTERS.length;b++){
				freq[x][LETTERS[b]]/=sum[x];
			}
			for (int b=0;b<LETTERS.length;b++){
				freq[x][LETTERS[b]] = (freq[x][LETTERS[b]]+pseudoCount)/(1+(4*pseudoCount));
				pwm[x][LETTERS[b]] = (float) (Math.log((double)freq[x][LETTERS[b]] / 0.25) / Math.log(2.0));
		}}
		
		
		//Save as WeightMatrix
		WeightMatrix motif = new WeightMatrix(pwm);
		return motif;
	}
	/**
	 * Calculate weighted motifs for each non-zero XL component
	 * @return
	 */
	protected Map<CompositeModelComponent, WeightMatrix> getPerComponentWeightMatrices(List<CompositeModelSiteAssignment> siteAssignments, Map<StrandedPoint, String> pointSeqs, int motifWidth, boolean motifWeightedByMaxCompsOnly){
		Map<CompositeModelComponent, WeightMatrix> xlComponentMotifs = new HashMap<CompositeModelComponent, WeightMatrix>();
		double maxWeightProp = 0.05; //cap per-site motif weights as 5% of total component tag count
		
		for(CompositeModelComponent xlComp : model.getXLComponents()){ if(xlComp.isNonZero()){
				int compOffset = xlComp.getPosition();
				int motifStartPos = compOffset - (motifWidth/2)+1;
				float[][] freq =  new float[motifWidth][WeightMatrix.MAXLETTERVAL];
				float[][] freqWeight =  new float[motifWidth][WeightMatrix.MAXLETTERVAL];
				float[][] pwm = new float[motifWidth][WeightMatrix.MAXLETTERVAL];
				float[] sum =  new float[motifWidth];
				float[] sumWeight =  new float[motifWidth];
				float pseudoCount=(float) 0.01;
				for(int x=0; x<motifWidth; x++){ 
					sum[x]=0; sumWeight[x]=0; 
					for (int b=0;b<LETTERS.length;b++){
						freq[x][LETTERS[b]]=0; freqWeight[x][LETTERS[b]]=0;
				}}
				int numCountedSites=0;
				
				//Sum of weights (for capping weights)
				float sumWeights=0;
				for(CompositeModelSiteAssignment sa : siteAssignments){
					for(ExperimentCondition cond : manager.getConditions()){
						sumWeights+= sa.getCompResponsibility(cond, xlComp.getIndex());
					}
				}
				
				//Weighted freq matrix
				for(CompositeModelSiteAssignment sa : siteAssignments){
					StrandedPoint pt = sa.getPoint();
					String currSeq = pointSeqs.get(pt).toUpperCase();
					
					for(ExperimentCondition cond : manager.getConditions()){
						//This block weights XL sites by proportion of total XL tags
						/*
						double currWeight =0;
						double currResp = sa.getCompResponsibility(cond, xlComp.getIndex());
						double totXLResp =0;
						if(currResp>10){
							for(CompositeModelComponent comp : model.getXLComponents()){ if(comp.isNonZero()){
								totXLResp+=sa.getCompResponsibility(cond, comp.getIndex());
							}}
							currWeight=currResp/totXLResp;
						}*/
						
						double currWeight =0;
						if(motifWeightedByMaxCompsOnly){
							//This block weights XL sites as component tags iff they are strongest XL site
							double currResp = sa.getCompResponsibility(cond, xlComp.getIndex());
							double maxXLResp =0;
							for(CompositeModelComponent comp : model.getXLComponents()){ if(comp.isNonZero()){
								if(sa.getCompResponsibility(cond, comp.getIndex())>maxXLResp)
									maxXLResp =sa.getCompResponsibility(cond, comp.getIndex()); 
							}}
							//currWeight=(currResp==maxXLResp && currResp>1) ? currResp : 0;
							currWeight=(currResp==maxXLResp && currResp>1) ? 1 : 0;
						}else{						 
							//Weight = XL tags at this component
							double currResp = sa.getCompResponsibility(cond, xlComp.getIndex());
							currWeight = currResp>1 ? sa.getCompResponsibility(cond, xlComp.getIndex()) : 0;
							currWeight/=sumWeights;
							if(currWeight>maxWeightProp)
								currWeight=maxWeightProp;
						}
						
						if(currWeight>0)
							numCountedSites++;
						
						for(int x=0; x<motifWidth; x++){ 
							char letter = currSeq.charAt(motifStartPos+x);
							freqWeight[x][letter]+=currWeight;
							sumWeight[x]+=currWeight;
							freq[x][letter]++;
							sum[x]++;
						}
					}
				}
				//Normalize freq matrix
				//Log-odds over expected motif
				for(int x=0; x<motifWidth; x++){ 
					for (int b=0;b<LETTERS.length;b++){
						freq[x][LETTERS[b]]/=sum[x];
						freqWeight[x][LETTERS[b]]/=sumWeight[x];
					}
					for (int b=0;b<LETTERS.length;b++){
						freq[x][LETTERS[b]] = (freq[x][LETTERS[b]]+pseudoCount)/(1+(4*pseudoCount));
						freqWeight[x][LETTERS[b]] = (freqWeight[x][LETTERS[b]]+pseudoCount)/(1+(4*pseudoCount));
						pwm[x][LETTERS[b]] = (float) (Math.log((double)freqWeight[x][LETTERS[b]] / (double)freq[x][LETTERS[b]]) / Math.log(2.0));
				}}
				
				
				//Save as WeightMatrix
				WeightMatrix motif = new WeightMatrix(pwm);
				xlComponentMotifs.put(xlComp, motif);
				
				//Temp print
				System.out.println("XLComponent motif: "+(xlComp.getPosition()-signalComposite.getCenterOffset()));
				for(int x=0; x<motifWidth; x++){
					System.out.println(x+"\t"+pwm[x]['A']+"\t"+pwm[x]['C']+"\t"+pwm[x]['G']+"\t"+pwm[x]['T']);
					//System.out.println(x+"\t"+freqWeight[x]['A']+"\t"+freqWeight[x]['C']+"\t"+freqWeight[x]['G']+"\t"+freqWeight[x]['T']);
					//System.out.println(x+"\t"+freq[x]['A']+"\t"+freq[x]['C']+"\t"+freq[x]['G']+"\t"+freq[x]['T']);
				}System.out.println("");
				
				System.out.println("Component "+(xlComp.getPosition()-signalComposite.getCenterOffset())+": motif created using "+numCountedSites);
			}
		}		
		
		return xlComponentMotifs;
	}
	
	/**
	 * Attempt to correct XL tag counts for XL bias.
	 *  The intuition here is if an XL component is taking proportionally more responsibility for tags than
	 *  it does on average, it must be compensating for inefficient crosslinking at other components.
	 *  Therefore, this method infers total XL tag counts from the over-strength XL components.
	 *    
	 * @param siteAssignments
	 */
	protected void correctXLBias(ExperimentCondition cond, List<CompositeModelSiteAssignment> siteAssignments, File correctedOutFile){
		try{
			FileWriter fout = new FileWriter(correctedOutFile);
			fout.write("#Point\tObsXLTags\tEstXLTags\n");
			
			//Find the expected/average XL proportions.
			double totalXLpi = 0;
			double[] expectedXLProps = new double[model.getXLComponents().size()];
			Map<CompositeModelComponent, Integer> component2Index = new HashMap<CompositeModelComponent, Integer>();
			int i=0;
			for(CompositeModelComponent xlComp : model.getXLComponents()){
				totalXLpi+=xlComp.getPi();
				expectedXLProps[i]=xlComp.getPi();
				component2Index.put(xlComp, i);
				i++;
			}
			for(int x=0; x<expectedXLProps.length; x++)
				expectedXLProps[x]/=totalXLpi;
			
			//Examine actual XL proportions for each ML assigned site  
			for(CompositeModelSiteAssignment sa : siteAssignments){
				StrandedPoint pt = sa.getPoint();
				
				//Calculate observed XL proportions & counts for this site
				double siteXLtotal = 0;
				double[] obsXLProps = new double[model.getXLComponents().size()];
				double[] obsXLCounts = new double[model.getXLComponents().size()];
				for(CompositeModelComponent xlComp : model.getXLComponents()){ if(xlComp.isNonZero()){
					siteXLtotal+=sa.getCompResponsibility(cond, xlComp.getIndex());
					obsXLProps[component2Index.get(xlComp)]=sa.getCompResponsibility(cond, xlComp.getIndex());
					obsXLCounts[component2Index.get(xlComp)]=sa.getCompResponsibility(cond, xlComp.getIndex());
				}}
				for(CompositeModelComponent xlComp : model.getXLComponents()){ if(xlComp.isNonZero()){
					obsXLProps[component2Index.get(xlComp)]/=siteXLtotal;
				}}
				
				//For over-strength XL components, estimate totals
				double sumEstTot = 0, count=0;
				for(CompositeModelComponent xlComp : model.getXLComponents()){ if(xlComp.isNonZero()){
					if(obsXLProps[component2Index.get(xlComp)] > expectedXLProps[component2Index.get(xlComp)]){
						sumEstTot+=obsXLCounts[component2Index.get(xlComp)]/expectedXLProps[component2Index.get(xlComp)];;
						count++;
					}
				}}
				
				double estimatedXLtotal = count>0 ? sumEstTot/count : siteXLtotal;
				fout.write(pt.toString()+"\t"+siteXLtotal+"\t"+estimatedXLtotal+"\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	

	/**
	 * Print per-XL-component bias factors at each site
	 *    
	 * @param siteAssignments
	 */
	protected void printXLBias(ExperimentCondition cond, List<CompositeModelSiteAssignment> siteAssignments, File outFile){
		try{
			FileWriter fout = new FileWriter(outFile);
			fout.write("Point\tTotalXLTags");
			for(CompositeModelComponent nz : model.getXLComponents()){ if(nz.isNonZero()){
				int pos = nz.getPosition() - model.getCenterOffset();
				fout.write("\t"+"Obs_"+nz.getLabel()+"_"+pos);
				fout.write("\t"+"Exp_"+nz.getLabel()+"_"+pos);
				fout.write("\t"+"Log2Ratio_"+nz.getLabel()+"_"+pos);
			}}fout.write("\n");
			
			//Find the expected/average XL proportions.
			double totalXLpi = 0;
			double[] expectedXLProps = new double[model.getXLComponents().size()];
			Map<CompositeModelComponent, Integer> component2Index = new HashMap<CompositeModelComponent, Integer>();
			int i=0;
			for(CompositeModelComponent xlComp : model.getXLComponents()){
				totalXLpi+=xlComp.getPi();
				expectedXLProps[i]=xlComp.getPi();
				component2Index.put(xlComp, i);
				i++;
			}
			for(int x=0; x<expectedXLProps.length; x++)
				expectedXLProps[x]/=totalXLpi;
			
			//Examine actual XL proportions for each ML assigned site  
			for(CompositeModelSiteAssignment sa : siteAssignments){
				StrandedPoint pt = sa.getPoint();
				
				//Calculate observed XL proportions & counts for this site
				double siteXLtotal = 0;
				double[] obsXLProps = new double[model.getXLComponents().size()];
				double[] obsXLCounts = new double[model.getXLComponents().size()];
				for(CompositeModelComponent xlComp : model.getXLComponents()){ if(xlComp.isNonZero()){
					siteXLtotal+=sa.getCompResponsibility(cond, xlComp.getIndex());
					obsXLProps[component2Index.get(xlComp)]=sa.getCompResponsibility(cond, xlComp.getIndex());
					obsXLCounts[component2Index.get(xlComp)]=sa.getCompResponsibility(cond, xlComp.getIndex());
				}}
				for(CompositeModelComponent xlComp : model.getXLComponents()){ if(xlComp.isNonZero()){
					obsXLProps[component2Index.get(xlComp)]/=siteXLtotal;
				}}
				
				fout.write(pt.toString()+"\t"+siteXLtotal);
				for(CompositeModelComponent xlComp : model.getXLComponents()){ if(xlComp.isNonZero()){
					double ratio = Math.log((obsXLProps[component2Index.get(xlComp)]+0.01) / (expectedXLProps[component2Index.get(xlComp)]+0.01))/Math.log(2);
					fout.write("\t"+obsXLProps[component2Index.get(xlComp)]+"\t"+expectedXLProps[component2Index.get(xlComp)]+"\t"+ratio);
				}}
				fout.write("\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Print indicators of S/W bases in the windows covered by XL points
	 * 	S=1, W=0, no-data or N = -1
	 *    
	 * @param 
	 */
	protected void printSWBases(File outFile, Map<StrandedPoint, String> pointSeqs, int seqWin){
		try{
			FileWriter fout = new FileWriter(outFile);
			fout.write("Point");
			for(int pos = -1*(seqWin/2); pos<(seqWin/2); pos++)
				fout.write("\tpos_"+pos);
			fout.write("\n");
			
			for(StrandedPoint pt : pointSeqs.keySet()){
				String currSeq = pointSeqs.get(pt).toUpperCase();
				
				fout.write(pt.toString());
				
				//Account for any offsets (i.e. seqWin goes beyond chromosome bounds) 
				int startOffset = (pt.getLocation()-(seqWin/2) < 0) ? (seqWin/2)-pt.getLocation() : 0;
				int endOffset = (pt.getLocation()+(seqWin/2) > gconfig.getGenome().getChromLength(pt.getChrom())) ? (pt.getLocation()+(seqWin/2)-gconfig.getGenome().getChromLength(pt.getChrom())) : 0;
				int i=0;
				for(i=0; i<startOffset; i++)
					fout.write("\t-1");
				for(i=i; i<seqWin-endOffset; i++){
					char letter = currSeq.charAt(i-startOffset);
					if(letter=='C' || letter=='G')
						fout.write("\t1");
					else if(letter=='A' || letter=='T')
					fout.write("\t0");
					else
						fout.write("\t-1");
				}
				for(i=0; i<endOffset; i++)
					fout.write("\t-1");
				
				fout.write("\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
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
			XLQuantifier xlQuant = new XLQuantifier(gcon, econ, ccon);
			xlQuant.execute();
		}
		
	}
}
