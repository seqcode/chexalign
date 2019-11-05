package org.seqcode.projects.chexalign;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Pair;
import org.seqcode.projects.chexalign.alignment.AlignmentConfig;
import org.seqcode.projects.chexalign.alignment.MultipleAlignment;
import org.seqcode.projects.chexalign.alignment.NeedlemanWunschAffine;
import org.seqcode.projects.chexalign.alignment.TagProfile;
import org.seqcode.projects.chexalign.multicompositemodel.CompositeXLFinderMultiCond;
import org.seqcode.projects.chexalign.xlanalysis.XLAnalysisConfig;

/**
 * Tree implements progressive profile alignment and UPGMA.
 * 
 * @author naomi
 */
public class ChExAlign {
	protected AlignmentConfig config;
	protected XLAnalysisConfig xlconfig;
	protected List<StrandedPoint> spts;
	protected int numPoints;
	protected int numCond;
	protected int win;
	protected double[][] pairwiseDist;
	protected double minScore;
	protected double maxScore;
	protected double gapScaling;
	protected NeedlemanWunschAffine[][] nwaffine;
	protected ExperimentManager manager;
	protected CompositeTagDistribution signalComposite;
	protected CompositeTagDistribution controlComposite=null;
	protected TreeNode[] nodes;	// Tree nodes
	protected MultipleAlignment[] alignmentRec;	//alignment records
	protected TagProfile[] tagProfiles;	// initial tag profiles
	protected TreeNode root;
	protected String filename; 
	protected List<Integer> anchorIndices = new ArrayList<Integer>();
	protected CompositeXLFinderMultiCond xlFinder;
	protected boolean normScore=true;
		
	public ChExAlign(AlignmentConfig config, XLAnalysisConfig xlcon, GenomeConfig gcon,ExptConfig econ, ExperimentManager eMan){
		this.config=config;
		xlconfig = xlcon;
		spts = config.getStrandedPoints();
		manager = eMan;
		win = config.getWindowSize();
		filename = config.getFilename();
		numPoints = spts.size();
		numCond = manager.getNumConditions();
		gapScaling = config.getGapExtScalingFactor();
		signalComposite = new CompositeTagDistribution(spts, eMan, win, true);
		// check for control experiment
		if (manager.getReplicates().get(0).hasControl())
			controlComposite = new CompositeTagDistribution(spts, eMan, win, false);
		for(ExperimentCondition cond : eMan.getConditions()){
			String compositeFileName = this.filename+"_composite."+cond.getName()+".txt";
			signalComposite.printProbsToFile(cond, compositeFileName);
		}	
		nodes = new TreeNode[numPoints];
		alignmentRec = new MultipleAlignment[numPoints];
		
		// Make all initial tag profiles
		tagProfiles = new TagProfile[numPoints];
		for (int i=0; i < numPoints; i++){
			tagProfiles[i]= new TagProfile(win, numCond, signalComposite.getPointWatsons(i), signalComposite.getPointCricks(i), 
					controlComposite==null? null: controlComposite.getPointWatsons(i)[0], controlComposite==null? null: controlComposite.getPointCricks(i)[0]);	
			tagProfiles[i].setID(i);
			tagProfiles[i].setInitialCoords(spts.get(i));
		}
		
		// Do all pairwise alignment
		alignAllPairwise();
		
		xlFinder = new CompositeXLFinderMultiCond(gcon, econ, manager, xlconfig);
	}
	
	public void alignAllPairwise(){
		// Initialize matrices
		double[][] pairwiseSimilarities = new double[numPoints][numPoints];
		pairwiseDist = new double[numPoints][numPoints];
		nwaffine = new NeedlemanWunschAffine[numPoints][numPoints];
		for (int i=0; i < numPoints; i++){
			for (int j=0; j < numPoints; j++){
				pairwiseSimilarities[i][j] = 0;
				pairwiseDist[i][j] = 0;
			}
		}
		
		// Compute pairwise similarities
		minScore= Double.MAX_VALUE; maxScore = config.MINIMUM_VALUE;
		for (int i=0; i<numPoints; i++){
			for (int j=i+1; j< numPoints; j++){
				nwaffine[i][j]=NeedlemanWunschAffine.AlignPair(manager, config, tagProfiles[i], tagProfiles[j], config.getGapPenalty(), gapScaling);
				double score = nwaffine[i][j].getMaxScore();
				pairwiseSimilarities[i][j] = score;
				pairwiseSimilarities[j][i] = score;
				if (score < minScore)
					minScore=score;
				if (score > maxScore)
					maxScore = score;
			}
		}
		
		System.out.println("raw similarity scores");
		for (int i=0; i < numPoints; i++){
			for (int j=i+1; j< numPoints; j++)
				System.out.print(pairwiseSimilarities[i][j]+",");
			System.out.println();
		}
		
		// Normalize similarity matrices
		for (int i = 0; i <numPoints;i++){
			for (int j=i+1; j< numPoints; j++){
				pairwiseSimilarities[i][j] = (pairwiseSimilarities[i][j]-minScore)/(maxScore-minScore);
				pairwiseSimilarities[j][i] = (pairwiseSimilarities[j][i]-minScore)/(maxScore-minScore);
		}}		
		
		System.out.println("normalized similarity scores");
		for (int i=0; i < numPoints; i++){
			for (int j=i+1; j< numPoints; j++)
				System.out.print(pairwiseSimilarities[i][j]+",");
			System.out.println();
		}	
		
		//Convert similarities to distances		
		double generalMean=0; double count=0;
		for (int i=0; i < numPoints; i++){
			for (int j=i+1; j< numPoints; j++){
				double currScore = 1-pairwiseSimilarities[i][j];
				pairwiseDist[i][j] = currScore;
				pairwiseDist[j][i] = currScore;
				generalMean+=currScore;
				count++;			
			}
		}
		
		System.out.println("normalized distance");
		for (int i=0; i < numPoints; i++){
			for (int j=i+1; j< numPoints; j++)
				System.out.print(pairwiseDist[i][j]+",");
			System.out.println();
		}
		
		generalMean/=count;	
		System.out.println(generalMean);		
	}
	
	/**
     * Build UPGMA tree based on pairwise distances/similarities 
     */	
	public void buildTree(){			
		// Keep track of which node is active
		boolean[] active = new boolean[numPoints];
		//Put all profiles to TreeNode
		for (int i=0; i < numPoints; i++){
			nodes[i] = new TreeNode(true, i, 1);
			active[i] = true;
		}	
		
		TreeNode tmpNode = null;	
		MultipleAlignment multialign = null; // multiple alignment results
		// Iterate until all nodes are on the tree
		for (int z=0; z < (numPoints-1); z++){
			double minDist = Double.MAX_VALUE;
			int minNodeA=-1; int minNodeB=-1;
			//Step 1: Find the two nodes with the minimum dij
			for (int i=0; i< numPoints; i++){
				if (active[i]){
					for (int j=i+1; j < numPoints; j++){
						if (active[j] && (pairwiseDist[i][j] < minDist) ){
							minDist = pairwiseDist[i][j];
							minNodeA = i; minNodeB=j;
						}}}}	
			System.out.println("minimum nodes are "+minNodeA+" , "+minNodeB);
			
			//Step 2: put minNodes A&B into new Node
			double h = minDist/2.0;
			tmpNode = new TreeNode(z, nodes[minNodeA], nodes[minNodeB], h, nodes[minNodeA].getNumMembers()+nodes[minNodeB].getNumMembers());
			nodes[minNodeB].setParent(tmpNode);		
			nodes[minNodeA].setEdge(h - nodes[minNodeA].getHeight());
			nodes[minNodeB].setEdge(h - nodes[minNodeB].getHeight());
			
			// Make alignment
			ArrayList<Integer> leaflist = new ArrayList<Integer>();
			NeedlemanWunschAffine nw;
			// Both are leaves
			if (nodes[minNodeA].isLeaf() && nodes[minNodeB].isLeaf()){
				nw = nwaffine[minNodeA][minNodeB];
				alignmentRec[minNodeA] = new MultipleAlignment(1,nw.getAlignedLength(), numCond, controlComposite==null? false : true);
				alignmentRec[minNodeA].setTagProfile(0, tagProfiles[minNodeA]);
				alignmentRec[minNodeA] = alignmentRec[minNodeA].SingleProfileAddition(nw, tagProfiles[minNodeB]);				
				leaflist.add(nodes[minNodeA].getLeafID());
				leaflist.add(nodes[minNodeB].getLeafID());
				anchorIndices.add(minNodeA);
			//A is a leaf; B is a profile
			}else if (nodes[minNodeA].isLeaf()){
				leaflist = nodes[minNodeB].getLeafList();
				TagProfile profileB = alignmentRec[minNodeB].Alignment2Profile();
				nw = NeedlemanWunschAffine.AlignPair(manager, config, profileB, tagProfiles[minNodeA], computeGapPenalty(1), gapScaling);
				alignmentRec[minNodeA] = alignmentRec[minNodeB].SingleProfileAddition(nw, tagProfiles[minNodeA]);
				leaflist.add(nodes[minNodeA].getLeafID());		
			// A is a profile; B is a leaf
			}else if (nodes[minNodeB].isLeaf()){
				leaflist = nodes[minNodeA].getLeafList();
				TagProfile profileA = alignmentRec[minNodeA].Alignment2Profile();
				nw = NeedlemanWunschAffine.AlignPair(manager, config, profileA, tagProfiles[minNodeB], computeGapPenalty(1), gapScaling);
				alignmentRec[minNodeA] = alignmentRec[minNodeA].SingleProfileAddition(nw, tagProfiles[minNodeB]);
				leaflist.add(nodes[minNodeB].getLeafID());		
			// Both are profiles
			}else{	
				leaflist = nodes[minNodeA].getLeafList();
				ArrayList<Integer> leaflistB = nodes[minNodeB].getLeafList();
				TagProfile profileA = alignmentRec[minNodeA].Alignment2Profile();
				TagProfile profileB = alignmentRec[minNodeB].Alignment2Profile();
				nw = NeedlemanWunschAffine.AlignPair(manager, config, profileA, profileB, computeGapPenalty(Math.min(leaflist.size(), leaflistB.size())), gapScaling);
				alignmentRec[minNodeA] = alignmentRec[minNodeA].MultipleProfileAddition(nw, alignmentRec[minNodeB]);
				leaflist.addAll(leaflistB);			
			}
			tmpNode.setLeafList(leaflist);		
						
			//Step 3: Add tmpNode into minNodeA's spot, inactivate minNodeB's spot, update dists for new node
			active[minNodeB] = false;
			for (int j=0; j < numPoints; j++){
				if (minNodeA!=j && active[j]){				
					if (config.useUPGMA())
						updateDistanceUsingUPGMA(j, minNodeA, minNodeB);
					else
						updateDistanceByAlignment(j, minNodeA, leaflist);
				}
			}
			nodes[minNodeA] = tmpNode;
			multialign = alignmentRec[minNodeA];
			///////////Calinski & Harabasz///////////////////////////
	        //Internal similarity of the nodes & between clusters
			/**
			System.out.println("distance matrix for iteration : "+z);
			for (int i=0; i < pairwiseDist.length; i++){
				for (int j=0; j < pairwiseDist[i].length; j++)
					System.out.print(pairwiseDist[i][j]+",");
				System.out.println();
			}
			**/
			
		}				
		
		//Tree building finished, the last node (i.e. tmpNode), should be the root
		root = tmpNode;
		
		System.out.println("tree building complete");
		
		// print alignment results
		multialign.printOriginalRegionsToFile(filename, win, config.useSortForPrint());
		multialign.printAlignedRegionsToFile(filename, config.useSortForPrint());
		
		MultipleAlignment.printOriginalTagsToFile(manager, signalComposite, filename);
		multialign.printAlignedTagsToFile(manager, filename, config.useSortForPrint());
		
		multialign.printAlignedCompositeToFile(manager, filename);
		
		if (config.doXLAnalysis()){
			//perform cross-linking analysis
			Pair<double[][][], double[][][]> perPointTags = multialign.makePerPointCountsFromTagProfiles();
			Pair<double[], double[]> ctrlCompositeTags = multialign.makeControlCompositeCounts();
			xlFinder.executeOnAlignedCompoiste(multialign.getAlignedLength(), spts, perPointTags.car(), perPointTags.cdr(),ctrlCompositeTags.car(), ctrlCompositeTags.cdr());
		}
		
		/**
		// this is not working
		System.out.println(root.getNodeID()+"\t"+root.getLeafList().toString()+"\t"+root.getNumMembers());
		System.out.println("right");
		System.out.println(root.getRight().getNodeID()+"\t"+root.getLeafID());
		System.out.println(root.getRight().getLeafList().toString());
		System.out.println("right of right");
		System.out.println(root.getRight().getRight().getLeafList().toString());
		System.out.println("left");
		System.out.println(root.getLeft().getLeafList().toString());
		System.out.println("left of left");
		System.out.println(root.getLeft().getLeft().getLeafList().toString());
		System.out.println("right of left");
		System.out.println(root.getLeft().getRight().getLeafList().toString());
		System.out.println("right");
		System.out.println(root.getRight().getLeafList().toString());
		**/
	
	}
	
	// Update distances using UPGMA calculations
	public void updateDistanceUsingUPGMA(int j, int minNodeA, int minNodeB){
		double aMem = nodes[minNodeA].getNumMembers();
		double bMem = nodes[minNodeB].getNumMembers();
		System.out.println("pairwiseDist[minNodeA][j] "+pairwiseDist[minNodeA][j]+", aMem "+aMem+", pairwiseDist[minNodeB][j] "+pairwiseDist[minNodeB][j]+" ,bMem "+bMem);
		pairwiseDist[minNodeA][j] = (pairwiseDist[minNodeA][j]*aMem + pairwiseDist[minNodeB][j]*bMem)/(aMem +bMem);
		pairwiseDist[j][minNodeA] = (pairwiseDist[j][minNodeA]*aMem + pairwiseDist[j][minNodeB]*bMem)/(aMem +bMem);	
		System.out.println("new distance for index "+minNodeA+" and "+j+ "is : "+pairwiseDist[minNodeA][j]);		
	}
	
	// Calculate distance by aligning a new profile to the rest
	public void updateDistanceByAlignment(int j, int minNodeA, ArrayList<Integer> leaflist){
		NeedlemanWunschAffine nw;
		if (nodes[j].isLeaf()){
			nw = NeedlemanWunschAffine.AlignPair(manager, config, alignmentRec[minNodeA].Alignment2Profile(), tagProfiles[j], computeGapPenalty(1), gapScaling);
		}else{
			double minNumLeaves = Math.min(leaflist.size(), nodes[j].getLeafList().size());
			// Gradient gap penality
			System.out.println("currGapPenalty is: "+computeGapPenalty(minNumLeaves));
			nw = NeedlemanWunschAffine.AlignPair(manager, config, 
					alignmentRec[minNodeA].Alignment2Profile(), alignmentRec[j].Alignment2Profile(), computeGapPenalty(minNumLeaves), gapScaling);
		}
		// convert similarity to distance
		double nSimilarity = (nw.getMaxScore()-minScore)/(maxScore-minScore);
		double distance = 1-nSimilarity;
		pairwiseDist[minNodeA][j] = distance;
		pairwiseDist[j][minNodeA] = distance;
		System.out.println("new distance for index "+minNodeA+" and "+j+ "is : "+pairwiseDist[minNodeA][j]);		
	}
	
	public double computeGapPenalty(double minNumLeaves){
		double currGapPenalty=config.getGapPenalty();
		if (config.useGradientGap() && minNumLeaves>0)
			currGapPenalty = config.getGapPenalty()*(1-(1-config.getMinGapScalingFactor())*minNumLeaves/((double) numPoints*0.5));			
		return currGapPenalty;	
	}
	
	public int getAnchorIndex(List<Integer> indices){
		int anchorIndex=-1;
		for (int index : indices){
			if (anchorIndices.contains(index)){
				anchorIndex=index; break;
			}
		}
		return anchorIndex;
	}
	
	public void PostorderListChildren(TreeNode n, TreeNode start){
		if (n.getLeft()!=null)
			PostorderListChildren(n.getLeft(), start);
		if (n.getRight()!=null)
			PostorderListChildren(n.getRight(), start);
		if (n.isLeaf()){
			
			
			
		}
		
	}
	
	private class TreeNode{
		protected boolean leaf=false;
		protected int leafID=-1;
		protected int nodeID=-1;
		protected TreeNode left;
		protected TreeNode right;
		protected TreeNode parent;
		protected double members=0;
		protected ArrayList<Integer> leafList;
		protected double height=0;
		protected double edge;
		
		/**
	     * Constructor for leaf. 
	     * @param leaf 
	     * @param leafID
	     * @param number of progeny members
	     */
		public TreeNode(boolean leaf, int leafID, double members){
			this.leaf=leaf;
			this.leafID=leafID;
			this.members=members;
		}
		
		/**
	     * Constructor for non-leaf. 
	     * @param nodeID 
	     * @param TreeNode left
	     * @param TreeNode right
	     * @param height of tree
	     * @param number of progeny members
	     */
		public TreeNode(int nodeID, TreeNode left, TreeNode right, double height, double members){
			this.nodeID = nodeID;
			this.left=left;
			this.right=right;
			this.height=height;
			this.members=members;		
			// Perform deep copy of TreeNode
			this.left=deepCopyNode(left);
			this.right=deepCopyNode(right);
		}
		
		// accessors
		public boolean isLeaf(){return leaf;}
		public int getLeafID(){return leafID;}
		public int getNodeID(){return nodeID;}
		public TreeNode getLeft(){return left;}
		public TreeNode getRight(){return right;}
		public TreeNode getParent(){return parent;}
		public double getHeight(){return height;}
		public double getNumMembers(){return members;}	
		public ArrayList<Integer> getLeafList(){return leafList;}
		
		// setters
		public void setParent(TreeNode parent){this.parent = parent;}
		public void setEdge(double edge){this.edge = edge;}
		public void setLeafList(ArrayList<Integer> leafList){this.leafList=leafList;}
		
		// Deep copy of node class
		public TreeNode deepCopyNode(TreeNode node){
			TreeNode copiedNode=null;
			if (node.isLeaf())
				copiedNode = new TreeNode(true,node.getLeafID(),node.getNumMembers());
			else
				copiedNode = new TreeNode(node.getNodeID(), node.getLeft(), node.getRight(),node.getHeight(),node.getNumMembers());
			return copiedNode;
		}
		
	}
	
	private class Child{
		
		public Child(){
			
		}
		
	}
	
	
	/**
	 * Main driver method
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args){
		
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);						
		AlignmentConfig config = new AlignmentConfig(gcon, args);			
		if (!config.useReadFilter())
			econ.setPerBaseReadFiltering(false);
		econ.setLoadRead2(false);
		
		if(args.length==0){
			System.err.println("ProgressiveProfileAlignment:"+
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
			
			XLAnalysisConfig ccon = new XLAnalysisConfig(gcon, args);					
			ChExAlign nodes = new ChExAlign(config,ccon, gcon, econ, manager);	
			
//			Tree nodes = new Tree(config,gcon, econ, manager);
			nodes.buildTree();
					
			manager.close();
		}
	}
}