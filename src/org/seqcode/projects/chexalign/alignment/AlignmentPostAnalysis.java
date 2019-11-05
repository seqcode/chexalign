package org.seqcode.projects.chexalign.alignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.BackgroundModelIO;
import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.RandomSequenceGenerator;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScorer;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.motifs.FreqMatrixImport;
import org.seqcode.motifs.MarkovMotifThresholdFinder;

/**
 * AlignmentPostAnalysis : a class to output matrix of 1) motif scores, 2) genomic annotation positions, 3) aligned sequences from tag alignments.
 * 
 * @author nuy11
 *
 */

public class AlignmentPostAnalysis {
	protected Genome gen=null;
	protected WeightMatrix motif;
	protected WeightMatrixScorer scorer;
	protected MarkovBackgroundModel backMod;
	protected static int numTest=100000;
	protected int win=300;
	protected SequenceGenerator seqgen;
	protected double minThreshold=0.1;
	protected String outFileName;
	protected ArrayList<String[]> coords;
	protected boolean revComplement=false;
	protected boolean strandedScore=true;
	protected boolean singlePosAnalysis=false;
	protected boolean makesummary=false;
	protected int numRegs;
	protected int len;
	
	public AlignmentPostAnalysis(Genome genome, ArrayList<String[]> coords, boolean useCache, boolean rc, String seqPath, String outName, boolean makesummary, boolean singlepos){
		gen=genome;
		this.coords = coords;
		outFileName = outName;
		seqgen = new SequenceGenerator();
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.setGenomePath(seqPath);
		}
		revComplement=rc;
		numRegs = coords.size();
		len = coords.get(0).length;
		this.makesummary = makesummary;
		singlePosAnalysis = singlepos;
	}
	
	public int[] getPositions(int i){
		int[] pos = new int[len];
		pos[0]=Integer.parseInt(coords.get(i)[0].split(":")[1]);
		for (int j=1; j < len; j++)
			pos[j]=Integer.parseInt(coords.get(i)[j]);
		return pos;
	}
	
	public StrandedRegion getStrandedRegion(int i, int[] pos){
		String chrm=coords.get(i)[0].split(":")[0].replace("chr", "");
		// get the start and end positions
		int start=-1; int end=-1;
		for (int j=0; j < len ; j++){if (pos[j]!=-1){ start=pos[j];break;}}
		for (int j=len-1; j >=0; j--){if (pos[j]!=-1){end=pos[j];break;}}
		// check direction			
		StrandedRegion reg = null;
		char strand;
		if (start<end){ reg = new StrandedRegion(gen,chrm,Math.max(0, start), Math.min(end, gen.getChromLength(chrm)),'+');}
		else{ reg = new StrandedRegion(gen, chrm, Math.max(0, end), Math.min(start, gen.getChromLength(chrm)),'-');}
		return reg;
	}
	
	private void maxToArray(int startPos, int stopPos, int[] pos, double[] array, double value){		
		for (int k=startPos; k < stopPos; k++)
			for (int i=0; i< pos.length; i++)
				if (pos[i]==k)
					array[i]= Math.max(array[i], value);
	}
	
	private void printScores(double[][] scores){
		try{
    		FileWriter fout = new FileWriter(outFileName+".out");
    		for (int i=0; i < scores.length; i++){
    			for (int j=0; j < scores[i].length; j++)
    				fout.write(scores[i][j]+",");
    			fout.write("\n");
    		}
    		fout.close();
    	}catch (IOException e) {
    		// TODO Auto-generated catch block
    		e.printStackTrace();
    	}	
	}
	
	private void printSummary(double[][] scores){
		float[] blankcounts=new float[len];
		float[] scoresum = new float[len];
		for (int j=0; j < len; j++){
			blankcounts[j]=0;
			scoresum[j]=0;
			for (int i=0; i < scores.length; i++){
				if (scores[i][j] < 0)
					blankcounts[j]++;
				else if (scores[i][j] > 0)
					scoresum[j]+=scores[i][j];
			}
		}
		// normalize
		for (int j=0; j < len; j ++){
			blankcounts[j]/=((float) numRegs);
			scoresum[j]/=((float) numRegs);
		}
		try{
			FileWriter fout = new FileWriter(outFileName+".summary.out");
			for (int j=0; j < len ; j++)
				fout.write((1-blankcounts[j])+",");
			fout.write("\n");
			for (int j=0; j < len ; j++)
				fout.write(scoresum[j]+",");
			fout.write("\n");
			fout.close();
		}catch (IOException e) {
    		// TODO Auto-generated catch block
    		e.printStackTrace();
    	}	
	}
	
	public void executeMotifMaker(WeightMatrix wm, MarkovBackgroundModel back, double minthres, boolean stranded){	
		motif = wm;
		minThreshold=minthres;
		backMod = back;
		if (revComplement){
			motif = WeightMatrix.reverseComplement(motif);
		}
		scorer = new WeightMatrixScorer(motif);
		strandedScore = stranded;
		
		//Pre-load the random sequences
        ArrayList<String> randSeq = new ArrayList<String>();
        RandomSequenceGenerator gen = new RandomSequenceGenerator(backMod);
		for(int i=0; i<numTest; i++){
			randSeq.add(gen.execute(win));
		}
		
		MarkovMotifThresholdFinder finder = new MarkovMotifThresholdFinder(motif, backMod, numTest);
		finder.setWin(win);
    	finder.setRandomSeq(randSeq);
    	double thres = finder.execute(minThreshold);
    	System.out.println("motif threshold : "+thres);
		
    	// get motif score
    	double[][] scores = new double[numRegs][len];
		for (int i=0; i < numRegs; i++)
			for (int j=0 ; j < len; j ++)
				scores[i][j]=-1;	
		for(int i=0; i < numRegs; i++){
			int[] pos = getPositions(i);
			StrandedRegion reg = getStrandedRegion(i, pos);
			
			String seq = seqgen.execute(reg);
			WeightMatrixScoreProfile profiler = scorer.execute(seq);
			
			for (int z=reg.getStart(); z < reg.getEnd(); z++){
				int currIndex = z-reg.getStart();
				double currScore =-1;
				if (strandedScore){
					currScore = reg.getStrand()  == '+' ? profiler.getForwardScores()[currIndex] : profiler.getReverseScores()[currIndex];				
				}else{
					currScore = profiler.getMaxScore(currIndex);				
				}
				if (!singlePosAnalysis){
					int startPos = z;
					int stopPos = z+motif.length()-1;
					maxToArray(startPos, stopPos, pos, scores[i], currScore>= thres ? currScore : 0);
				}else{	// record a score at single position
					int midPos = z+motif.length()/2;
					for (int l=0; l < pos.length; l++)
						if (pos[l] == midPos)
							scores[i][l]=currScore>= thres ? currScore : 0;		
				}}}		
		
		printScores(scores);
		if (makesummary){printSummary(scores);}
	}
	
	public void PrintSequencesFromAlignedRegions(){
		char[][] mappedSequences = new char[numRegs][len];		
		for (int i=0; i < numRegs; i++)
			for (int j=0 ; j < len; j ++)
				mappedSequences[i][j]='-';	
		for(int i=0; i < numRegs; i++){			
			int[] pos = getPositions(i);
			StrandedRegion reg = getStrandedRegion(i, pos);		
			String seq = seqgen.execute(reg);
			for (int z=reg.getStart(); z < reg.getEnd(); z++){
				int currIndex = z-reg.getStart();
				for (int l=0; l < pos.length; l++){
					if (pos[l] == z){
						char base = seq.charAt(currIndex);		
						mappedSequences[i][l] = reg.getStrand() == '+' ? base : SequenceUtils.complementChar(base);
					}}}}
		
		if (revComplement){
			// copy all sequences to temp and delete sequences from mappedSequences
			for (int i=0; i <numRegs; i++){
				char [] temp = new char[len];
				for (int j=0; j < len; j ++){
					temp[j]=mappedSequences[i][j];
					mappedSequences[i][j]='-';
				}		
				// reverse order and get reverse complement
				for (int j=0; j < len ; j ++){
					char currChar = temp[len-j];	
					if (currChar != '-')
						mappedSequences[i][j] = SequenceUtils.complementChar(currChar);
				}}}	
			
		try{
	    	FileWriter fout = new FileWriter(outFileName+".fa");
	    	for (int k=0; k < mappedSequences.length; k++){
	    		fout.write(">seq"+k+"\n");
	    		for (int j=0; j < mappedSequences[k].length; j++)
	    			fout.write(mappedSequences[k][j]);
	    		fout.write("\n");
	    	}
	    	fout.close();
	    }catch (IOException e) {
	    	// TODO Auto-generated catch block
	    	e.printStackTrace();
	    }
	}
	
	public void mapAnnotation(List<StrandedPoint> spts, int width){
		double[][] scores = new double[numRegs][len];
		int sameStrandCount=0; int oppositeStrandCount =0;
			
		for(int i=0; i < numRegs; i++){
			int[] pos = getPositions(i);
			for (int l=0; l < pos.length ; l++){
				if (pos[l] >=0)
					scores[i][l]=0;
				else
					scores[i][l]=-1;
			}
			
			StrandedRegion reg = getStrandedRegion(i, pos);			

			for (StrandedPoint p : spts){
				if (reg.contains(p)){	
					if (width == -1 || singlePosAnalysis){
						for (int l=0; l < pos.length; l++)
							if (pos[l] == p.getLocation())					
								scores[i][l]+=1;					
					}else{
						int startPos = p.getLocation()-((int) width/2);
						int endPos = p.getLocation()+ ((int) width/2);
						maxToArray(startPos, endPos, pos, scores[i],1);			
					}
					if (reg.getStrand() == p.getStrand())
						sameStrandCount +=1;
					else
						oppositeStrandCount +=1;
				}}		
		}
		
		System.out.println("same strand: "+sameStrandCount+"\n"+"opposite strand: "+oppositeStrandCount);
					
		printScores(scores);
		if (makesummary){printSummary(scores);}
	}
	
	public static void main(String[] args) throws IOException, ParseException {
		ArgParser ap = new ArgParser(args);	
		GenomeConfig gcon = new GenomeConfig(args);
		Genome gen = gcon.getGenome();
		String posFileName = Args.parseString(args, "posf", null);
		String motifName = Args.parseString(args,"motif", null);
		String backName = Args.parseString(args,"mback", null);
		double minthres = Args.parseDouble(args, "mthres", 0);
		String outName = Args.parseString(args, "out", "motifout");
		boolean useCache = Args.parseFlags(args).contains("cache") ? true : false;
		boolean rc = Args.parseFlags(args).contains("rc") ? true : false;
		boolean strandedScore = Args.parseFlags(args).contains("notstranded") ? false : true;
		boolean singlepos = Args.parseFlags(args).contains("singlepos") ? true: false;
		boolean makesummary = Args.parseFlags(args).contains("summary") ? true: false;
		String seqPathName="";
		if(useCache){
			seqPathName = Args.parseString(args, "seq", "");
		}
		List<StrandedPoint> spts = null;
		if (ap.hasKey("points"))
			spts = RegionFileUtilities.loadStrandedPointsFromFile(gen, Args.parseString(args, "points", null));
		int width = Args.parseInteger(args, "w", -1);
		
		if(gen==null || posFileName==null){printError();}
						
		// motif maker		
		ArrayList<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
		MarkovBackgroundModel back = null;
		if (motifName != null){		
			//Load the background
			back = BackgroundModelIO.parseMarkovBackgroundModel(backName, gen);
			//Load the motifs
			FreqMatrixImport motifImport = new FreqMatrixImport();
			motifImport.setBackground(back);
			for(WeightMatrix wm : motifImport.readTransfacMatrices(motifName)){
				motifs.add(wm);
			}
		}	
		
		BufferedReader reader = new BufferedReader(new FileReader(posFileName));
		ArrayList<String[]> coords = new ArrayList<String[]>();
		String line;
		while ((line = reader.readLine())!=null){
			line = line.trim();
			String[] words = line.split(",");
			if (words.length > 1)
				coords.add(words);
		}	
		
		AlignmentPostAnalysis analysis = new AlignmentPostAnalysis(gen, coords, useCache, rc, seqPathName, outName, makesummary, singlepos);
		if (!motifs.isEmpty())
			analysis.executeMotifMaker(motifs.get(0), back, minthres, strandedScore);
		else if (spts!=null)
			analysis.mapAnnotation(spts, width);			
		else		
			analysis.PrintSequencesFromAlignedRegions();
	}

	private static void printError(){
		System.err.println("Usage: AlignmentPostAnalysis --species <organism;genome> --seq <sequence path path>\n" +
				"--motif <motif names> --mback <background model name> --mthres <threshold>\n" +
				"--posf <coordinate file name> --out <output root name> \n" +
				"--color <red/green/blue> \n");
		System.exit(1);
	}
}
