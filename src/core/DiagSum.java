/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.util.Arrays;

/**
 *
 * @author ben
 */
public class DiagSum {
    
    //sum(PP*)
    private float[] diagsum=null;
    //#words effectivily summed in this diagsum
    private int[] effectiveWordCount=null;
    //corresponding position in reference
    private int[] referencePosition=null;

    
    private int queryLength=-1;
    private int referenceLength=-1;
    private int minOverlap=-1;
    private int k=-1;
    private int s=-1;
    

    
    /**
     * 
     * @param queryLength
     * @param referenceLength
     * @param minOverlap minimum overlap required between query and reference to build a diagSum score
     * @param k word size
     * @param s step used in query words cut
     */
    public DiagSum(int queryLength, int referenceLength, int minOverlap, int k, int s) {
        this.queryLength=queryLength;
        this.referenceLength=referenceLength;
        this.minOverlap=minOverlap;
        this.k=k;
        this.s=s;
        this.diagsum=new float[queryLength-minOverlap+referenceLength-minOverlap];
        this.effectiveWordCount=new int[queryLength-minOverlap+referenceLength-minOverlap];
        this.referencePosition=new int[queryLength-minOverlap+referenceLength-minOverlap];
    }

    public void sum(int diagSumPosition, float val) {
        this.diagsum[diagSumPosition]+=val;
        this.effectiveWordCount[diagSumPosition]+=1;
    }


    public float getSum(int diagSumPosition) {
        return this.diagsum[diagSumPosition];
    }

    public int getEffectiveWordCount(int diagSumPosition) {
        return this.effectiveWordCount[diagSumPosition];
    }
    
    public int getPositionInReference(int diagSumPosition) {
        return this.referencePosition[diagSumPosition];
    }
    
    
    
    /**
     * init all positions of the DiagSum, as the #words*PPStratThreshold
     * @param PPStarThreshold 
     */
    public void init(float PPStarThreshold) {
        for (int i = 0; i < diagsum.length; i++) {
            diagsum[i]=getMaxWordCount(i)*PPStarThreshold;
            referencePosition[i]=i-(queryLength-minOverlap);
        }
    }
    
    
    
    /**
     * return the maximum number of words that could have been summed in 
     * a particular DiagSum position
     * @param diagSumPosition
     * @return 
     */
    public int getMaxWordCount(int diagSumPosition) {
        if ( diagSumPosition<(queryLength-minOverlap) ) {
            return((queryLength-(queryLength-minOverlap)+diagSumPosition-k+1)/s);
        }
        if ( diagSumPosition>((queryLength-minOverlap)+referenceLength-queryLength-1) ) {
            return (referenceLength+1-diagSumPosition+(queryLength-minOverlap)-(k+1))/s;
        }            
        return (queryLength-k+1)/s;      
    }
    
    
    public int getSize() {
        return this.diagsum.length;
    }

    @Override
    public String toString() {
        
        return "PP*:"+Arrays.toString(diagsum)+"\n#Words:"+Arrays.toString(effectiveWordCount);
        
    }
    
    
    
    
}
