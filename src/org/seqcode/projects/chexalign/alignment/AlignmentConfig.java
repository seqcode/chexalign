package org.seqcode.projects.chexalign.alignment;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

/**
 * AlignmentConfig: 
 * 		Maintains all constants needed by ChIP-exo tag alignment. 
 *     
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class AlignmentConfig {
	
	public static String version = "0.1";
	protected boolean printHelp=false;
	protected GenomeConfig gconfig;
	protected Genome gen=null;
	List<StrandedPoint> spts;
	protected int w;
	protected double gap_penalty;	//gap penalty of NeedlemanWunsch affine alignment
	protected double minGapScaling;	//Scaling factor to determine minimum gap penalty when using gradient gap penalty
	protected double gapExtensionScalingFactor;	//Scaling factor to determine gap extension penalty relative to gap opening
	protected double smoothingGaussianVariance;
	protected boolean upgma=false;	// Use UPGMA guide tree to align
	protected boolean gradient=true;	//Use gradient gap penalty
	protected boolean smooth=false; 	//Smooth profiles before alignment
	protected boolean normalizeBySignal =false;	//Normalize profile based on signal proportions calculated by NCIS scaling ratio
	protected boolean ihstrans =false; // Transform data using inverse hyperbolic sine function
	protected boolean doXLAnalysis = true;
	protected boolean sortforprint = true; // sort tag profiles before printing
	protected boolean normalize = true;
	protected boolean doReadFilter=false;	// Turn on per base read filter in case of highly duplicated experiment
	protected String filename;
	
	
	/** Similarity metrics for NeedlemanWunsch Alignment    **/
	public boolean linear = false;
	public boolean euclidean = true; 	// 1. Minkowski family (default)
	public boolean sorensen = false; 	// 2. L1 family
	public boolean soergel = false;		// 2. L1 family
	public boolean lorentzian = false;	// 2. L1 family
	public boolean cosine = false;		// 3. Inner product family
	public boolean pce = false;			// 3. Inner product family
	public boolean chisquare = false;	// 4. Squared L2 family or Chi-squre family
	public boolean divergence = false;	// 4. Squared L2 family or Chi-squre family
	public boolean clark = false;		// 4. Squared L2 family or Chi-squre family
	public boolean kl = false;			// 5. Shannon's entropy family
	
	public boolean pearson = false;		// Pearson correlation
	
	public boolean sum_similarity = false;		// Calculate single similarity score across expt rather than sum of individual experiments
	public boolean debugMode=false;		// print extra outputs for debugging 
	
	// Constants
	public final double MIN_OVERLAP=0.5;	// minimum overlap required for Needleman-Wunsch alignment
	public final double GAPP_FRAC=0.5;	// minimum overlap required for Needleman-Wunsch alignment
	public final double MINIMUM_VALUE = -10000;
	public final double NUM_ITR=2;
	
	protected String[] args;
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	
	public AlignmentConfig(GenomeConfig gcon, String[] arguments){
		gconfig = gcon;
		gen = gconfig.getGenome();
		this.args=arguments; 
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h")){
			printHelp=true;			
		}else{
				
			// parse command line arguments
			w = Args.parseInteger(args, "cwin", 201);
			// convert to odd number to avoid 1 bp off
			if (w % 2==0){ w+=1;}
			gap_penalty = Args.parseDouble(args, "gap", 5.0);
			//If gradient gap is used, specify the minimum gap penalty
			minGapScaling = Args.parseDouble(args, "mingapscaling", 0.1);
			if (minGapScaling > 1){
				System.err.println("Minimum gap penalty should be smaller than the gap open penalty. Setting the mingapscaling to be 1");
				minGapScaling=1.0;
			}
			//Fraction of gap penalty used for gap extension
			gapExtensionScalingFactor = Args.parseDouble(args, "extscaling", 0.1);
			if (gapExtensionScalingFactor > 1){
				System.err.println("Gap extension penalty should be smaller than the gap open penalty. Setting the extscaling to be 1");
				gapExtensionScalingFactor=1.0;
			}
			
			smoothingGaussianVariance = Args.parseDouble(args, "gausssmoothparam", 1.0);
			//Use UPGMA guide tree for multiple alignment
			upgma = Args.parseFlags(args).contains("upgma") ? true : false;
			// Use constant gap penalty
			gradient = Args.parseFlags(args).contains("constantgap") ? false : true;
			// Smooth profiles before alignment
			smooth = Args.parseFlags(args).contains("smooth") ? true : false;
			// Normalize tags by maximum values within windows
			normalizeBySignal = Args.parseFlags(args).contains("normalizebysig") ? true : false;
			// Transform data using inverse hyperbolic sine function
			ihstrans = Args.parseFlags(args).contains("ihstrans")? true : false;
			// Print extra outputs for debugging
			debugMode = Args.parseFlags(args).contains("debug") ? true : false;
			// Sort tag profile based on the order of input region order, this option will turn off  sorting
			sortforprint = Args.parseFlags(args).contains("nosort") ? false : true;
			
			spts = new ArrayList<StrandedPoint>();
			if (ap.hasKey("cpoints")){
				spts = RegionFileUtilities.loadStrandedPointsFromFile(gen, Args.parseString(args, "cpoints", null));
			}else if (ap.hasKey("points")){
				List<Point> pts = RegionFileUtilities.loadPointsFromFile(Args.parseString(args, "points", null), gen);
				for (Point p : pts)
					spts.add(new StrandedPoint(p,'+'));
			}else{
				System.err.println("Please provide --cpoints <stranded point file> OR --points <point file>");
			}
					
			filename = Args.parseString(args, "out", "out");
			
			/** Similarity metrics for NeedlemanWunsch Alignment    **/
			linear = Args.parseFlags(args).contains("linear") ? true : false;
			// 2. L1 family
			sorensen = Args.parseFlags(args).contains("sorensen") ? true : false;
			soergel = Args.parseFlags(args).contains("soergel") ? true : false;
			lorentzian = Args.parseFlags(args).contains("lorentzian") ? true : false;
			// 3. Inner product family
			cosine = Args.parseFlags(args).contains("cosine") ? true : false;
			pce = Args.parseFlags(args).contains("pce") ? true : false;
			// 4. Squared L2 family or Chi-squre family
			chisquare = Args.parseFlags(args).contains("chisquare") ? true : false;
			divergence = Args.parseFlags(args).contains("divergence") ? true : false;
			clark = Args.parseFlags(args).contains("clark") ? true : false;
			// 5. Shannon's entropy family
			kl = Args.parseFlags(args).contains("kl") ? true : false;
			// Pearson correlation coefficient
			pearson = Args.parseFlags(args).contains("pearson") ? true : false;
			
			//calculate similarity score and sum across experiments, rather than calculate a single score across experiments
			sum_similarity = Args.parseFlags(args).contains("sumsimilarity") ? true : false;
			
			//Do not perform cross-linking analysis after tag alignment
			doXLAnalysis = Args.parseFlags(args).contains("noxlanalysis") ? false : true;
			
			// Turn on per base read filtering
			doReadFilter = Args.parseFlags(args).contains("readfilter") ? true : false;	
		}
	}
	
	//Accessors
	public boolean helpWanted(){return printHelp;}
	public int getWindowSize(){return w;}
	public double getGapPenalty(){return gap_penalty;}
	public double getGapExtScalingFactor(){return gapExtensionScalingFactor;}
	public double getMinGapScalingFactor(){return minGapScaling;}
	public String getFilename(){return filename;}
	public List<StrandedPoint> getStrandedPoints(){return spts;}
	public boolean doXLAnalysis(){return doXLAnalysis;}
	public boolean useUPGMA(){return upgma;}
	public boolean useGradientGap(){return gradient;}
	public boolean useSmoothing(){return smooth;}
	public boolean useSigNormalize(){return normalizeBySignal;}
	public boolean useInverseHyperbolicSineTrans(){return ihstrans;}
	public boolean useSortForPrint(){return sortforprint;}
	public boolean useNormalize(){return normalize;}
	public double getSmoothingGaussianVariance(){return smoothingGaussianVariance;}
	public boolean useReadFilter(){return doReadFilter;}

}
