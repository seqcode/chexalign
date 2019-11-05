package org.seqcode.projects.chexalign.alignment;

/**
 * SimilarityScore: Computes similarity scores based on distance metrics described in the following paper.
 * Comprehensive Survey on Distance/Similarity Measures between Probability Density Functions: Sung-Hyuk Cha
 * 
 * @author naomi yamada
 *
 */

public class SimilarityScore {	
	protected AlignmentConfig config;
	protected int numCond;
	
	public SimilarityScore(AlignmentConfig config, int numC){
		numCond = numC;
		this.config = config;
	}
	
	public double computeScore(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double score = 0.0;
		double dist=0.0;
		double penalty=0.0;
		boolean euclidean=false;
		
		if (config.pearson){
			score =(double) numCond*pearsonCC(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			if (config.debugMode)
				System.out.print("score: "+score+",");
		}else{
			if (config.sorensen)
				dist=sorensen(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			else if (config.soergel)
				dist=soergel(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			else if (config.cosine)
				dist=cosine(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			else if (config.pce)
				dist=PCE(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			else if (config.divergence)
				dist=divergence(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			else if (config.clark)
				dist=clark(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);		
			else if (config.linear)
				dist=linear(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			else if (config.lorentzian)
				dist=lorentzian(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			else if (config.chisquare)
				dist=chiSquare(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			else if (config.kl)
				dist=KL(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
			else{
				dist=euclidean(aIndex, bIndex, nWatsonA, nWatsonB, nCrickA, nCrickB);
				euclidean=true;
			}
			
			if (euclidean){
				for (int e=0; e < numCond; e++)
					penalty+=(Math.abs(nWatsonA[e][aIndex]-nWatsonB[e][bIndex]) + Math.abs(nCrickA[e][aIndex]-nCrickB[e][bIndex]));
				score = (double) numCond -dist-penalty;
			}else{
				score = (double) numCond*(1 - dist); // penalty is zero if similarity is calculated using vectorized profiles
			}
		}
		if (config.debugMode)
			System.out.print("dist: "+dist+",penalty: "+penalty+","+"score: "+score+",");
		if (Double.isNaN(dist) || Double.isNaN(score)){
			System.out.println();
			System.out.println("NAN value found");
			System.out.println("aIndex,"+aIndex+ " bIndex,"+bIndex);
			System.out.println("nWatsonA");
			for (int e=0; e < numCond; e++)
				System.out.print(nWatsonA[e][aIndex]+",");
			System.out.println();
			System.out.println("nWatsonB[e][bIndex]");
			for (int e=0; e < numCond; e++)
				System.out.print(nWatsonB[e][bIndex]+",");
			System.out.println("nCrickA[e][aIndex]");
			for (int e=0; e < numCond; e++)
				System.out.print(nCrickA[e][aIndex]+",");
			System.out.println("nCrickB[e][bIndex]");
			for (int e=0; e < numCond; e++)
				System.out.print(nCrickB[e][bIndex]+",");
		}
		
		/**	
		double score = 0.0;
		for (int e=0; e< numCond; e++){
			double dist=0.0;
			// score distributions indicate that not subtracting penalties works better for the following distance metrics
			if (config.sorensen||config.soergel||config.cosine||config.pce||config.divergence||config.clark){
				if (config.sorensen)
					dist=sorensen(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				if (config.soergel)
					dist=soergel(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				if (config.cosine)
					dist=cosine(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				if (config.pce)
					dist=PCE(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				if (config.divergence)
					dist=divergence(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				if (config.clark)
					dist=clark(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				score += (1 - dist);
				if (config.debugMode)
					System.out.print("score: "+score+",");
			}else{			
				if (config.linear)
					dist=linear(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				else if (config.lorentzian)
					dist=lorentzian(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				else if (config.pearson)
					dist=pearson(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				else if (config.chisquare)
					dist=chiSquare(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				else if (config.kl)
					dist=KL(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);
				else
					dist=euclidean(nWatsonA[e][aIndex], nWatsonB[e][bIndex], nCrickA[e][aIndex], nCrickB[e][bIndex]);

				double penalty = Math.abs(nWatsonA[e][aIndex]-nWatsonB[e][bIndex]) + Math.abs(nCrickA[e][aIndex]-nCrickB[e][bIndex]);
				score += (1 - dist - penalty);
				if (config.debugMode)
					System.out.print("dist: "+dist+",penalty: "+penalty+",");
			}
		}
	**/
		return score;
	}
	
	// Linear
	protected double linear(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double doubledist=0.0;
		for (int e=0; e < numCond; e++)
			doubledist += (nWatsonA[e][aIndex] + nWatsonB[e][bIndex] + nCrickA[e][aIndex] + nCrickB[e][bIndex]);
		return (doubledist/2);
	}
	
	// 1. Minkowski family
	protected double euclidean(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){	
		double dist=0.0;
		if (config.sum_similarity){
			for (int e=0; e< numCond; e++)
				dist+=Math.sqrt(Math.pow(nWatsonA[e][aIndex] - nWatsonB[e][bIndex], 2) + Math.pow(nCrickA[e][aIndex] - nCrickB[e][bIndex], 2));
		}else{
			double sqdist=0.0;
			for (int e=0; e < numCond; e++)
				sqdist += (Math.pow(nWatsonA[e][aIndex] - nWatsonB[e][bIndex], 2) + Math.pow(nCrickA[e][aIndex] - nCrickB[e][bIndex], 2));
			dist=Math.sqrt(sqdist);
		}
		return dist;
	}
	
	// 2. L1 family
	protected double sorensen(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){	
		double numerator=0.0;
		double denom = 0.0;				
		for (int e=0; e < numCond; e++){
			numerator += (Math.abs(nWatsonA[e][aIndex] - nWatsonB[e][bIndex]) + Math.abs(nCrickA[e][aIndex] - nCrickB[e][bIndex]));
			denom += (nWatsonA[e][aIndex] + nWatsonB[e][bIndex] + nCrickA[e][aIndex] + nCrickB[e][bIndex]);
		}
		if (denom == 0)
			return 1;
		else
			return (numerator/denom);
	}	
	protected double soergel(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double numerator=0.0;
		double denom = 0.0;				
		for (int e=0; e < numCond; e++){
			numerator += (Math.abs(nWatsonA[e][aIndex] - nWatsonB[e][bIndex]) + Math.abs(nCrickA[e][aIndex] - nCrickB[e][bIndex]));
			denom += (Math.max(nWatsonA[e][aIndex], nWatsonB[e][bIndex])+ Math.max(nCrickA[e][aIndex], nCrickB[e][bIndex]));
		}
		if (denom == 0)
			return 1;
		else
			return (numerator/denom);
	}
	protected double lorentzian(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double dist=0.0;
		for (int e=0; e < numCond; e++)
			dist += (Math.log(1+Math.abs(nWatsonA[e][aIndex] - nWatsonB[e][bIndex])) + Math.log(1+Math.abs(nCrickA[e][aIndex] - nCrickB[e][bIndex])));
		return dist;		
	}
	
	// 3. Inner product family
	protected double cosine(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double numerator=0.0;
		double sqdenomA=0.0;
		double sqdenomB=0.0;
		for (int e=0; e < numCond; e++){
			numerator += (nWatsonA[e][aIndex]*nWatsonB[e][bIndex] + nCrickA[e][aIndex]*nCrickB[e][bIndex]);
			sqdenomA += (Math.pow(nWatsonA[e][aIndex], 2) + Math.pow(nCrickA[e][aIndex] , 2));
			sqdenomB += (Math.pow(nWatsonB[e][bIndex], 2) + Math.pow(nCrickB[e][bIndex], 2));
		}
		if (sqdenomA != 0 && sqdenomB !=0)
			return (numerator/(Math.sqrt(sqdenomA)*Math.sqrt(sqdenomB)));
		else
			return 1;
	}
	// Kumar-Hassebrook (PCE)
	protected double PCE(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double numerator = 0.0;
		double denom = 0.0;
		for (int e=0; e < numCond; e++){
			numerator += (nWatsonA[e][aIndex]*nWatsonB[e][bIndex] + nCrickA[e][aIndex]*nCrickB[e][bIndex]);
			denom += (Math.pow(nWatsonA[e][aIndex], 2) + Math.pow(nCrickA[e][aIndex], 2) + Math.pow(nWatsonB[e][bIndex], 2) + Math.pow(nCrickB[e][bIndex], 2)
				- (nWatsonA[e][aIndex]*nWatsonB[e][bIndex]+nCrickA[e][aIndex]*nCrickB[e][bIndex]));
		}
		if (denom != 0)
			return (numerator/denom);
		else
			return 1;
	}
	
	// 4. Squared L2 family or Chi-squre family
	protected double chiSquare(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double dist=0.0;
		for (int e=0; e<numCond; e++){
			if (nWatsonA[e][aIndex]+nWatsonB[e][bIndex] !=0)
				dist += (Math.pow(nWatsonA[e][aIndex]-nWatsonB[e][bIndex], 2)/(nWatsonA[e][aIndex]+nWatsonB[e][bIndex]));
			if (nCrickA[e][aIndex]+nCrickB[e][bIndex] != 0)
				dist += (Math.pow(nCrickA[e][aIndex]-nCrickB[e][bIndex], 2)/(nCrickA[e][aIndex]+nCrickB[e][bIndex]));
		}
		return dist;
	}
	protected double divergence(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double halfdist=0.0;
		for (int e=0; e< numCond; e++){
			if (nWatsonA[e][aIndex]+nWatsonB[e][bIndex] !=0)
				halfdist += (Math.pow(nWatsonA[e][aIndex]-nWatsonB[e][bIndex], 2)/Math.pow(nWatsonA[e][aIndex]+nWatsonB[e][bIndex], 2));
			if (nCrickA[e][aIndex]+nCrickB[e][bIndex] !=0)
				halfdist += (Math.pow(nCrickA[e][aIndex]-nCrickB[e][bIndex], 2)/Math.pow(nCrickA[e][aIndex]+nCrickB[e][bIndex], 2));
		}
		return halfdist*2;
	}
	protected double clark(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double sqdist=0.0;
		for (int e=0; e< numCond; e++){
			if (nWatsonA[e][aIndex]+nWatsonB[e][bIndex] !=0)
				sqdist += (Math.pow(Math.abs(nWatsonA[e][aIndex]-nWatsonB[e][bIndex])/(nWatsonA[e][aIndex]+nWatsonB[e][bIndex]), 2));
			if (nCrickA[e][aIndex]+nCrickB[e][bIndex] !=0)
				sqdist += (Math.pow(Math.abs(nCrickA[e][aIndex]-nCrickB[e][bIndex])/(nCrickA[e][aIndex]+nCrickB[e][bIndex]), 2));
		}
		return Math.sqrt(sqdist);	
	}
	
	// 5. Shannon's entropy family
	// Kullback-Leibler
	protected double KL(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double dist=0.0;
		for (int e=0; e< numCond; e++){
			if (nWatsonA[e][aIndex]!=0 && nWatsonB[e][bIndex]!=0){
				dist += (nWatsonA[e][aIndex]*Math.log(nWatsonA[e][aIndex]/nWatsonB[e][bIndex]));
				dist += (nWatsonB[e][bIndex]*Math.log(nWatsonB[e][bIndex]/nWatsonA[e][aIndex]));
			}
			if (nCrickA[e][aIndex]!=0 && nCrickB[e][bIndex]!=0){
				dist += (nCrickA[e][aIndex]*Math.log(nCrickA[e][aIndex]/nCrickB[e][bIndex]));
				dist += (nCrickB[e][bIndex]*Math.log(nCrickB[e][bIndex]/nCrickA[e][aIndex]));
			}
		}
		return dist;	
	}
	
	// Pearson's correlation coefficient
	protected double pearsonCC(int aIndex, int bIndex, double[][] nWatsonA, double[][] nWatsonB, double[][] nCrickA, double[][] nCrickB){
		double sumA=0.0, aveA=0.0, sumB=0.0, aveB=0.0;
		for (int e=0; e< numCond; e++){
			sumA += (nWatsonA[e][aIndex]+nCrickA[e][aIndex]);
			sumB += (nWatsonB[e][bIndex]+nCrickB[e][bIndex]);
		}
		aveA=sumA/((double) numCond*2);
		aveB=sumB/((double) numCond*2);
		
		double cov =0.0, varA = 0.0, varB = 0.0;
		for (int e=0; e< numCond; e++){
			double wa = nWatsonA[e][aIndex] - aveA;
			double wb = nWatsonB[e][bIndex] - aveB;
			cov += wa*wb;
			varA += wa*wa;
			varB += wb*wb;
			double ca = nCrickA[e][aIndex] - aveA;
			double cb = nCrickB[e][bIndex] - aveB;
			cov += ca*cb;
			varA += ca*ca;
			varB += cb*cb;
		}
		double corr=0;
		if (Math.sqrt(varA)*Math.sqrt(varB)!=0)
			corr=cov/(Math.sqrt(varA)*Math.sqrt(varB));
		return corr;
	}
	
	/**
	// Linear
	protected double linear(double x1, double x2, double y1, double y2){
		return ((x1 + x2)/2 + (y1 + y2)/2);	
	}
	
	// 1. Minkowski family
	protected double euclidean(double x1, double x2, double y1, double y2){	
		return (Math.sqrt(Math.pow(x1-x2, 2) + Math.pow(y1-y2, 2)));
	}
	
	// 2. L1 family
	protected double sorensen(double x1, double x2, double y1, double y2){	
		if (x1+x2+y1+y2 > 0)
			return ((Math.abs(x1-x2) + Math.abs(y1-y2))/(x1+x2+y1+y2));	
		else
			return 1;
	}	
	protected double soergel(double x1, double x2, double y1, double y2){	
		if (x1+x2+y1+y2 > 0)
			return ((Math.abs(x1-x2) + Math.abs(y1-y2))/(Math.max(x1, x2) + Math.max(y1, y2)));	
		else
			return 1;
	}	
	protected double lorentzian(double x1, double x2, double y1, double y2){	
		return ((Math.log(1+Math.abs(x1-x2))+Math.log(1+Math.abs(y1-y2))));		
	}
	
	// 3. Inner product family
	protected double cosine(double x1, double x2, double y1, double y2){
		if (x1+y1>0 && x2+y2 > 0)
			return ((x1*x2 + y1*y2)/(Math.sqrt(Math.pow(x1, 2)+Math.pow(y1, 2))*Math.sqrt(Math.pow(x2, 2)+Math.pow(y2, 2))));
		else
			return 1;
	}	
	// Kumar-Hassebrook (PCE)
	protected double PCE(double x1, double x2, double y1, double y2){	
		if (x1+x2+y1+y2 > 0)
			return ((x1*x2 + y1*y2)/(Math.pow(x1, 2) + Math.pow(y1, 2) + Math.pow(x2, 2) + Math.pow(y2, 2) - (x1*x2 + y1*y2)));
		else
			return 1;
	}
	
	// 4. Squared L2 family or Chi-squre family
	protected double pearson(double x1, double x2, double y1, double y2){
		double dist=0.0;
		if (x1+x2>0)
			dist=(Math.pow(x1-x2, 2)*2)/(x1+x2);
		if (y1+y2>0)
			dist+=(Math.pow(y1-y2, 2)*2)/(y1+y2);
		return dist;				
	}
	protected double chiSquare(double x1, double x2, double y1, double y2){
		double dist=0.0;
		if (x1+x2>0)
			dist= Math.pow(x1-x2, 2)/(x1+x2);
		if (y1+y2>0)
			dist += Math.pow(y1-y2, 2)/(y1+y2);
		return dist;
	}
	protected double divergence(double x1, double x2, double y1, double y2){		
		double dist = 0.0;
		if (x1+x2> 0)
			dist = 2*Math.pow(x1-x2, 2)/Math.pow(x1+x2,2);
		if (y1+y2> 0)
			dist += 2*Math.pow(y1-y2, 2)/Math.pow(y1+y2,2);
		return dist*2;
	}
	protected double clark(double x1, double x2, double y1, double y2){		
		double dist = 0.0;
		if (x1+x2 >0)
			dist = Math.pow(Math.abs(x1-x2)/(x1+x2),2);
		if (y1+y2 >0)
			dist+=Math.pow(Math.abs(y1-y2)/(y1+y2),2);
		return Math.sqrt(dist);
	}	
	
	// 5. Shannon's entropy family
	// Kullback-Leibler
	protected double KL(double x1, double x2, double y1, double y2){
		double dist=0.0;
		if (x1>0 && x2>0){
			dist+=x1*Math.log(x1/x2);
			dist+=x2*Math.log(x2/x1);
		}
		if (y1>0 && y2>0){
			dist+=y1*Math.log(y1/y2);
			dist+=y2*Math.log(y2/y1);	
		}
		return dist;
	}
	**/
}
