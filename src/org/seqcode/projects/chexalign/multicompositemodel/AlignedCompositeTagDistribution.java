package org.seqcode.projects.chexalign.multicompositemodel;


import java.util.List;

import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.location.StrandedPoint;

public class AlignedCompositeTagDistribution extends CompositeTagDistribution {
	
	public AlignedCompositeTagDistribution(){
		super();
	}
	
	public double[][] makeComposite(double[][][] perPointCounts, double[] ctrlComposite){
		double[][] composite= new double[numConditions][win];
		for (int p=0; p < numPoints; p++)
			for (int c=0; c < numConditions; c++)
				for (int w=0; w < win; w++)
					composite[c][w]+=perPointCounts[p][c][w];
		
		//test start//
		for (int c=0; c < numConditions; c++){
			if (ctrlComposite!=null){
				System.out.println("scaled control composite");
				for (int w=0; w < win; w++)
					System.out.print(ctrlComposite[w]*exptMan.getIndexedCondition(c).getPooledSampleControlScaling()+",");
				System.out.println();
			}
		
			System.out.println("before subtracting");
			for (int w=0; w < win; w++)
				System.out.print(composite[c][w]+",");
			System.out.println();
		}// test end
		
		
		if (ctrlComposite!=null){
			for (int c=0; c < numConditions; c++){
				double scaling = exptMan.getIndexedCondition(c).getPooledSampleControlScaling();
				for (int w=0; w < win; w++){
					double currVal=composite[c][w]-ctrlComposite[w]*scaling;
					composite[c][w] = currVal>0? currVal: 0;
				}}}
		
		//test start
		for (int c=0; c < numConditions; c++){
			System.out.println("after subtracting");
			for (int w=0; w < win; w++)
				System.out.print(composite[c][w]+",");
			System.out.println();
			
		}
		//test end
		
		return composite;	
	}
	
	public void normalizeComposite(){
		//Normalize
		for (int c=0; c < numConditions;c++){
			double sum=0;
			for(int w=0; w<win; w++){
				sum+=watson[c][w]; sum+=crick[c][w];
			}for(int w=0; w<win; w++){
				watson[c][w]/=sum; crick[c][w]/=sum;
			}
		}
	}
	
	// Setter
	public void setPerPointWatson(double[][][] perPointWatson, double[] ctrlWatsonComposite){
		this.perPointWatson = perPointWatson; 
		watson = makeComposite(perPointWatson, ctrlWatsonComposite);
	}	
	public void setPerPointCrick(double[][][] perPointCrick, double[] ctrlCrickComposite){
		this.perPointCrick = perPointCrick; 
		crick = makeComposite(perPointCrick, ctrlCrickComposite);
	}
	
	public void initializeClass(List<StrandedPoint> points, ExperimentManager eMan, int win, boolean loadSignal){
		exptMan = eMan;
		this.win = win;
		centerOffset = win/2;
		this.numConditions=exptMan.getNumConditions();
		this.points = points;
		numPoints = points.size();
		isSignal = loadSignal;
		
		for(int p=0; p<numPoints; p++)
			pointIndex.put(points.get(p), p);
	}
}
