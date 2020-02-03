package org.seqcode.projects.chexalign.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.math.probability.NormalDistribution;

/**
 * NeedlemanWunschAffine: A global alignment algorithm to align two sets of profiles with an affine gap.
 * 
 * @author naomi yamada
 *
 */
public class NeedlemanWunschAffine {
	protected SimilarityScore sim;
	protected ExperimentManager manager;
	protected AlignmentConfig config;
	protected int numCond;			// number of experiments
	protected int indexA;
	protected int indexB;
	protected double[] gaussian;
	protected double[][] watsonA;	// watson tags from region A is indexed by genomic positions and experiments
	protected double[][] crickA;	// crick tags from region A is indexed by genomic positions and experiments
	protected double[][] watsonB;	// watson tags from region B is indexed by genomic positions and experiments
	protected double[][] crickB;	// crick tags from region B is indexed by genomic positions and experiments
	protected double[][] nWatsonA;	// normalized watson tags from region A is indexed by genomic positions and experiments
	protected double[][] nCrickA;	// normalized crick tags from region A is indexed by genomic positions and experiments
	protected double[][] nWatsonB;	// normalized watson tags from region B is indexed by genomic positions and experiments
	protected double[][] nCrickB;	// normalized crick tags from region B is indexed by genomic positions and experiments
	protected int[] gapsA;			// gap counts from region A is indexed by genomic positions
	protected int[] gapsB;			// gap counts from region B is indexed by genomic positions
	protected double[][] M;	// Ungapped alignment
	protected double[][] Ix;	// Gapped alignment with respect to X
	protected double[][] Iy; 	// Gapped alignment with respect to Y
	protected int[][] point_i;	// i pointer
	protected int[][] point_j;	// j pointer
	protected int x;	//Length of profile A
	protected int y;	//Length of profile B
	protected ArrayList<Integer> alignSectionX = new ArrayList<Integer>();
	protected ArrayList<Integer> alignSectionY = new ArrayList<Integer>();
	protected double maxScore;
	protected int alignLen=0;
	protected static double GAP_OPEN;
	protected static double GAP_EXT;	// GAP_EXT = GAP_OPEN * alpha
	protected boolean reverse=false;	//Alignment is performed using a reverse compliment of profile B
	
	public NeedlemanWunschAffine(ExperimentManager manager, AlignmentConfig config, int ia, int ib, double[][] watsonA, double[][] crickA, int[] gapsA, 
			double[][] watsonB, double[][] crickB, int[] gapsB, int aLen, int bLen, boolean reverse, double gap_penalty, double gapExtScalingFactor){
		this.config = config;
		this.manager = manager;
		this.watsonA= watsonA;
		this.crickA = crickA;
		this.watsonB = watsonB;
		this.crickB = crickB;
		this.gapsA = gapsA;
		this.gapsB = gapsB;
		this.numCond = manager.getNumConditions();
		this.reverse=reverse;
		indexA=ia;
		indexB=ib;
		x=aLen;
		y=bLen;
		initGaussianKernel(config.getSmoothingGaussianVariance());
		// initialize matrices and smooth profiles
		initMatrices(config.useSmoothing());
		sim = new SimilarityScore(config, numCond);
		
		//Normalize or transform data
		if (config.useSigNormalize())
			normalizeTagsBySignals();
		else if (config.useInverseHyperbolicSineTrans())
			inverseHyperbolicSineTransformation();
		else
			normalizeTagsByMax();
		
		// Adjust gap open penalty based on number of condition
		if (config.sum_similarity)
			GAP_OPEN=gap_penalty*(double) numCond;
		else
			GAP_OPEN=gap_penalty;
		
		// Compute gap extension penalty based on gap open penalty
		GAP_EXT=GAP_OPEN*gapExtScalingFactor;
	}
	
	// accessors
	public double getMaxScore(){return maxScore;}
	public int getAlignedLength(){return alignLen;}
	public boolean isReverse(){return reverse;}
	public int getIndexA(){return indexA;}
	public int getIndexB(){return indexB;}
	public double[][] getWatsonA(){return watsonA;}
	public double[][] getCrickA(){return crickA;}
	public double[][] getWatsonB(){return watsonB;}
	public double[][] getCrickB(){return crickB;}
	public int getALen(){return x;}
	public int getBLen(){return y;}
	public ArrayList<Integer> getAlignSectionX(){return alignSectionX;}
	public ArrayList<Integer> getAlignSectionY(){return alignSectionY;}
	
	// setters for testing
	public void setIndexB(int ib){indexB=ib;}
	
	public void initMatrices(boolean smoothProfile){
		//initialization of M, Ix, and Iy matrix
		M = new double [x+1][y+1];
		Ix = new double[x+1][y+1];
		Iy = new double[x+1][y+1];
		point_i = new int[x+1][y+1];
		point_j = new int[x+1][y+1];
		for (int i = 0 ; i <= x; i++){
			for (int j = 0 ; j <= y; j++){
				M[i][j]=0;
				Ix[i][j]=0;
				Iy[i][j]=0;
				point_i[i][j]=0;
				point_j[i][j]=0;
			}
		}
		// Edge initialization
		for (int i=1; i <=x; i++){
			point_i[i][0]=i-1;
			point_j[i][0]=0;
		}
		for (int j=1; j <=y; j++){
			point_i[0][j]=0;
			point_j[0][j]=j-1;
		}			
		
		//initialize and copy watson & crick profiles
		nWatsonA = new double[numCond][x];
		nCrickA = new double[numCond][x];
		nWatsonB = new double[numCond][y];
		nCrickB = new double[numCond][y];		
		for (int c=0; c < numCond; c++){
			for (int i=0; i < x; i++){
				nWatsonA[c][i] = watsonA[c][i];
				nCrickA[c][i] = crickA[c][i];
			}
			for (int i=0; i < y; i++){
				nWatsonB[c][i] = watsonB[c][i];
				nCrickB[c][i] = crickB[c][i];
			}
		}
		
		if (smoothProfile){
			for (int c=0; c < numCond ; c++){
				nWatsonA[c] = symmetricKernelSmoother(nWatsonA[c], gaussian);
				nCrickA[c] = symmetricKernelSmoother(nCrickA[c], gaussian);
				nWatsonB[c] = symmetricKernelSmoother(nWatsonB[c], gaussian);
				nCrickB[c] = symmetricKernelSmoother(nCrickB[c], gaussian);
			}	
		}
	}
	
	// normalize tags by maximum height
	public void normalizeTagsByMax(){
		// If cosine distance is used subtract mean 
		
		/**
		if (config.cosine){
			for (int c=0; c < numCond; c++){
				double sumA=0.0; double sumB=0.0;
				for (int i=0; i < al; i++){
					sumA+=nWatsonA[c][i];
					sumA+=nCrickA[c][i];
				}
				for (int i=0; i < al; i++){
					nWatsonA[c][i]-=(sumA/(double) al);
					nCrickA[c][i]-=(sumA/(double) al);
				}
				for (int i=0; i < bl; i++){
					sumB+=nWatsonB[c][i];
					sumB+=nCrickB[c][i];
				}
				for (int i=0; i < bl; i++){
					nWatsonB[c][i]-=(sumB/(double) bl);
					nCrickB[c][i]-=(sumB/(double) bl);
				}			
			}
		}
		**/
		
		// normalize
		for (int c=0; c < numCond; c++){
			double maxA=0, maxB=0;
			for (int i=0; i < x; i++){
				if (nWatsonA[c][i] > maxA)
					maxA = nWatsonA[c][i];
				if (nCrickA[c][i] > maxA)
					maxA = nCrickA[c][i];
			}
			if (maxA>0){
				for (int i=0; i < x; i++){
					nWatsonA[c][i] /= maxA;
					nCrickA[c][i] /= maxA;				
				}
			}
			for (int i=0; i < y; i++){
				if (nWatsonB[c][i] > maxB)
					maxB = nWatsonB[c][i];
				if (nCrickB[c][i] > maxB)
					maxB = nCrickB[c][i];
			}
			if (maxB>0){
				for (int i=0; i < y; i++){
					nWatsonB[c][i] /= maxB;
					nCrickB[c][i] /= maxB;				
				}
			}
		}
	}
	
	// normalize tags by signal strengths estimated by NCIS scaling
	public void normalizeTagsBySignals(){
		// 1. normalize each factor to a control experiment
		double maxA=0, maxB=0;
		for (ExperimentCondition cond: manager.getConditions()){
			int c = cond.getIndex();
			double norm = cond.getPooledSampleControlScaling();
			for (int i=0; i < x; i++){
				nWatsonA[c][i] /= norm;
				nCrickA[c][i] /= norm;	
				if (nWatsonA[c][i] > maxA)
					maxA = nWatsonA[c][i];
				if (nCrickA[c][i] > maxA)
					maxA = nCrickA[c][i];
			}
			for (int i=0; i < y; i++){
				nWatsonB[c][i] /= norm;
				nCrickB[c][i] /= norm;	
				if (nWatsonB[c][i] > maxB)
					maxB = nWatsonB[c][i];
				if (nCrickB[c][i] > maxB)
					maxB = nCrickB[c][i];
			}
		}		
		// 2. normalize using a maximum point across all experiments
		for (int c=0; c < numCond; c++){
			for (int i=0; i < x; i++){
				nWatsonA[c][i]/=maxA;
				nCrickA[c][i]/=maxA;
			}
			for (int i=0; i < y; i++){
				nWatsonB[c][i]/=maxB;
				nCrickB[c][i]/=maxB;
			}
		}
	}
	
	public void inverseHyperbolicSineTransformation(){
		for (int c=0; c < numCond; c++){
			for (int i=0; i < x; i++){
				nWatsonA[c][i]= inverseHyperbolicSine(nWatsonA[c][i]);
				nCrickA[c][i]=inverseHyperbolicSine(nCrickA[c][i]);
			}
			for (int i=0; i< y;i++){
				nWatsonB[c][i]=inverseHyperbolicSine(nWatsonB[c][i]);
				nCrickB[c][i]=inverseHyperbolicSine(nWatsonB[c][i]);
			}
		}
	}
	
	public static double inverseHyperbolicSine(double x){
		double y = Math.log(x+Math.sqrt(Math.pow(x, 2)+1));
		return y;
	}
	
	//pre-calculate and store the Guassian kernel prob., for efficiency
	private void initGaussianKernel(double var){
		gaussian = new double[250]; 
		NormalDistribution gaussianDist = new NormalDistribution(0, var*var);
		for (int i=0;i<gaussian.length;i++)
			gaussian[i]=gaussianDist.calcProbability((double)i);
	}
	
	// take a prob. density dist., smooth it using specified kernel density
	// assume kernel density is symmetrical (e.g. Gaussian).
	// assume X is evenly spaced.
	public static double[] symmetricKernelSmoother( double[]Y, double[]kernel){
		int length = Y.length;
		int kernel_length = kernel.length;
		double[] yy = new double[length];
		double total=0;
		for (int i=0;i<length;i++){
			double v=kernel[0]*Y[i] + 2.0E-300;		// init with very small number
            double weight=kernel[0];
            for (int j = 1; j < kernel_length && i+j < length; j++) {
                v+=Y[i+j]*kernel[j];
                weight += kernel[j];                
            }
            for (int j = 1; j < kernel_length && i-j >= 0; j++) {
                v+=Y[i-j]*kernel[j];
                weight += kernel[j];                
            }
			v = v / weight;
            yy[i] = v;
			total+=v;
		}
//		for (int i=0;i<length;i++){ // Do not normalize to be probabilities
//			yy[i]=yy[i]/total;
//		}
		return yy;
	}	
	
	// Free memories used for alignment matrices
	public void clear(){
		M=null;
		Ix=null;
		Iy=null;
		point_i=null;
		point_j=null;
	}

	public void buildMatrix(){	
		maxScore = config.MINIMUM_VALUE;
		int max_i=0; int max_j=0;
		double tmp=0; double currMax=0;	
		for (int i = 1 ; i <= x; i++){
			for (int j = 1 ; j <= y ; j++){	
				// M update
				M[i][j] = Math.max(M[i-1][j-1], Math.max(Ix[i-1][j-1], Iy[i-1][j-1]))+sim.computeScore(i-1,j-1, nWatsonA, nWatsonB, nCrickA, nCrickB);			
				point_i[i][j]=i-1;
				point_j[i][j]=j-1;
				currMax=M[i][j];
								
				// Ix update
				if (gapsA[i-1]!=0 || gapsB[j-1]!=0){ //Special case for partial gaps
					Ix[i][j] = M[i-1][j]-GAP_OPEN*config.GAPP_FRAC;
					tmp = Ix[i-1][j]-GAP_EXT*config.GAPP_FRAC;				
				}else{
					Ix[i][j] = M[i-1][j]-GAP_OPEN;
					tmp = Ix[i-1][j]-GAP_EXT;
				}
				if(Ix[i][j] < tmp && i>1 ){ Ix[i][j]=tmp;}
				if (Ix[i][j]> currMax){
					currMax = Ix[i][j];
					point_i[i][j]=i-1;
					point_j[i][j]=j;
				}
				
				//Iy update
				if (gapsA[i-1]!=0 || gapsB[j-1]!=0){ //Special case for partial gaps
					Iy[i][j] = M[i][j-1]-GAP_OPEN*config.GAPP_FRAC;
					tmp = Iy[i][j-1]-GAP_EXT*config.GAPP_FRAC;				
				}else{
					Iy[i][j] = M[i][j-1]-GAP_OPEN;
					tmp = Iy[i][j-1]-GAP_EXT;	
				}
				if (Iy[i][j] < tmp && j>1){Iy[i][j]=tmp;}
				if (Iy[i][j]> currMax){
					currMax=Iy[i][j];
					point_i[i][j]=i;
					point_j[i][j]=j-1;
				}	
				// Minimum overlap length enforced
				if (currMax > maxScore && (i==x && j > (int) Math.floor(config.MIN_OVERLAP*x) || j==y && j > (int) Math.floor(config.MIN_OVERLAP*y))){
					maxScore = currMax;
					max_i=i;
					max_j=j;
				}
			}
		}
		
		// Edge
		int start_i=max_i;
		int start_j=max_j;
		if (start_i==x){
			int i=start_i;
			for (int j=start_j+1; j<=y;j++){
				M[i][j]=maxScore;
				point_i[i][j]=i;
				point_j[i][j]=j-1;
			}
		}else{
			int j=start_j;
			for (int i=start_i+1; i <= x;i++){
				M[i][j]=maxScore;
				point_i[i][j]=i-1;
				point_j[i][j]=j;		
			}
		}
		start_i=x;
		start_j=y;	
		
		// Traceback
		tmp = maxScore;
		int tmp_si, tmp_sj;
		int z=0;
		while (start_i!=0 || start_j!=0){
			alignSectionX.add(start_i-1);
			alignSectionY.add(start_j-1);
			tmp_si=start_i; tmp_sj=start_j;
			start_i=point_i[tmp_si][tmp_sj];
			start_j=point_j[tmp_si][tmp_sj];
			z++;
		}	
		alignLen=z;
		
		// Print out content
		if (config.debugMode){
			System.out.println("max score: "+maxScore+"alignedLen: "+alignLen+",");
			System.out.println("alignSectionX: "+alignSectionX.toString());
			System.out.println("alignSectionY: "+alignSectionY.toString());
			printMatrix();
		}
		
		// Free memory used in arrays
		clear();
	}
	
	/**
	 * Perform a pairwise alignment
	 */	
	public static NeedlemanWunschAffine AlignPair(ExperimentManager manager, AlignmentConfig config, TagProfile profileA, TagProfile profileB, double gap_penalty, double gapExtFactor){
		
		int numCond = manager.getNumConditions();		
		double[][] watsonA = profileA.getWatsonProfile();
		double[][] crickA = profileA.getCrickProfile();
		int[] gapsA = profileA.getGaps();
		int aL = profileA.getAlignLength();
		int pA = profileA.getID();
		
		double[][] watsonB = profileB.getWatsonProfile();
		double[][] crickB = profileB.getCrickProfile();
		int[] gapsB = profileB.getGaps();
		int bL=profileB.getAlignLength();
		int pB = profileB.getID();
		double[][] rWatsonB = new double[numCond][bL];
		double[][] rCrickB = new double[numCond][bL];
		int[] rGapsB = new int[bL];
		
		//Reverse B tags
		for (int i=0; i < bL; i++){
			rGapsB[bL-i-1] = gapsB[i];
			for (int c=0; c < numCond; c++){
				rWatsonB[c][bL-i-1] = crickB[c][i];
				rCrickB[c][bL-i-1] = watsonB[c][i];
			}
		}	
		/**
		System.out.println("lengths "+aL+","+bL);		
		for (int c=0; c < numCond; c++){
			for (int i=0; i < watsonA[c].length; i++)
				System.out.print(watsonA[c][i]+",");
			System.out.println();
			for (int i=0; i < crickA[c].length; i++)
				System.out.print(crickA[c][i]+",");
			System.out.println();
			for (int i=0; i < watsonB[c].length; i++)
				System.out.print(watsonB[c][i]+",");
			System.out.println();
			for (int i=0; i < crickB[c].length; i++)
				System.out.print(crickB[c][i]+",");
			System.out.println();
			for (int i=0; i < rWatsonB[c].length; i++)
				System.out.print(rWatsonB[c][i]+",");
			System.out.println();
			for (int i=0; i < rCrickB[c].length; i++)
				System.out.print(rCrickB[c][i]+",");
			System.out.println();
		}
		**/
				
		// align using two possible ways
		// use estimated signals
		//System.out.println("pA: "+pA+"\tpB: "+pB);
		NeedlemanWunschAffine alignFor = new NeedlemanWunschAffine(manager,config, pA, pB, watsonA, crickA, gapsA, watsonB, crickB, gapsB, aL, bL, false, gap_penalty, gapExtFactor);
		NeedlemanWunschAffine alignRev = new NeedlemanWunschAffine(manager,config, pA, pB, watsonA, crickA, gapsA, rWatsonB, rCrickB, rGapsB, aL, bL, true, gap_penalty, gapExtFactor);
		alignFor.buildMatrix();
		alignRev.buildMatrix();
			
		if (alignFor.getMaxScore() > alignRev.getMaxScore()){
			alignRev=null;
			return alignFor;
		}else{
			alignFor=null;
			return alignRev;
		}
	}	
	
	// for testing only
	public void printMatrix(){
		System.out.println("printing M matrix");
		for (int j = 0; j <=y ; j++){
			for (int i = 0; i <=y; i++)
				System.out.print(M[i][j]+",");
			System.out.println();
		}
		System.out.println("printing Ix matrix");
		for (int j = 0; j <=y ; j++){
			for (int i = 0; i <=y; i++)
				System.out.print(Ix[i][j]+",");
			System.out.println();
		}
		System.out.println("printing Iy matrix");
		for (int j = 0; j <=y ; j++){
			for (int i = 0; i <=y; i++)
				System.out.print(Iy[i][j]+",");
			System.out.println();
		}
		System.out.println("printing pointer_i");
		for (int j = 0; j <=y ; j++){
			for (int i = 0; i <=y; i++)
				System.out.print(point_i[i][j]+",");
			System.out.println();
		}	
		System.out.println("printing pointer_j");
		for (int j = 0; j <=y ; j++){
			for (int i = 0; i <=y; i++)
				System.out.print(point_j[i][j]+",");
			System.out.println();
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
			System.err.println("NeedlemanWunschAffine:"+
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
			
			// parse command line arguments
			int w = Args.parseInteger(args, "cwin", 200);
			double gap = Args.parseDouble(args, "gap", 10);
			double gapScalingFactor = Args.parseDouble(args, "gapScaling", 0.5);
			List<StrandedPoint> spts = new ArrayList<StrandedPoint>();
			if (ap.hasKey("spoints")){
				spts = RegionFileUtilities.loadStrandedPointsFromFile(gcon.getGenome(), Args.parseString(args, "spoints", null));
			}else if (ap.hasKey("points")){
				List<Point> pts = RegionFileUtilities.loadPointsFromFile(Args.parseString(args, "points", null), gcon.getGenome());
				for (Point p : pts)
					spts.add(new StrandedPoint(p,'+'));
			}else{
				System.err.println("Please provide --spoints <stranded point file> OR --points <point file>");
			}
						
			CompositeTagDistribution maker = new CompositeTagDistribution(spts, manager, w, true);
			
			for (int i=0; i < maker.getPoints().size(); i++){
				for (int j=i+1; j < maker.getPoints().size(); j++){
					double[][] n_watson_a = maker.getPointWatsons(i);
					double[][] n_crick_a = maker.getPointCricks(i);
					double[][] n_watson_b = maker.getPointWatsons(j);
					double[][] n_crick_b = maker.getPointCricks(j);
					int[] gaps = new int[w];
					for (int l=0; l< w; l++)
						gaps[l]=0;
					for (int c=0; c < manager.getNumConditions(); c++){
						double maxA=0; double maxB=0;
						for (int l=0; l< w; l++){
							if (n_watson_a[c][l] > maxA){ maxA=n_watson_a[c][l];}
							if (n_crick_a[c][l] > maxA){ maxA=n_crick_a[c][l];}		
							if (n_watson_b[c][l] > maxB){ maxB=n_watson_b[c][l];}
							if (n_crick_b[c][l] > maxB) { maxB=n_crick_b[c][l];}
						}
						for (int l=0; l< w; l++){
							n_watson_a[c][l]/=maxA;
							n_crick_a[c][l]/=maxA;
							n_watson_b[c][l]/=maxB;
							n_crick_b[c][l]/=maxB;
						}
					}				
					
					if (config.debugMode){					
						for (int c=0; c < manager.getNumConditions(); c++){
							n_watson_a[c]=RandomizeArray(n_watson_a[c]);
							n_crick_a[c]=RandomizeArray(n_crick_a[c]);
							n_watson_b[c]=RandomizeArray(n_watson_b[c]);
							n_crick_b[c]=RandomizeArray(n_crick_b[c]);
						}						
					}						
										
					NeedlemanWunschAffine alignment = new NeedlemanWunschAffine(manager, config, i, j, 
							n_watson_a, n_crick_a, gaps, n_watson_b, n_crick_b,gaps, w,w, false, gap, gapScalingFactor);	
					alignment.buildMatrix();
				}
			}
			manager.close();
		}
	}
	
	public static double[] RandomizeArray(double[] array){
		Random rgen = new Random();  // Random number generator			
 
		for (int i=0; i<array.length; i++) {
		    int randomPosition = rgen.nextInt(array.length);
		    double temp = array[i];
		    array[i] = array[randomPosition];
		    array[randomPosition] = temp;
		}
 
		return array;
	}
}