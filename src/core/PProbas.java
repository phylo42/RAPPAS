/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.io.Serializable;
import java.util.Arrays;

/**
 *
 * @author ben
 */
public class PProbas implements Serializable {
    
    private static final long serialVersionUID = 3000L;
    
    double[][][] pp=null; //[nodeid][site][state]

    public PProbas(int nodeCount, int siteCount, int stateCount) {
        pp=new double[nodeCount][siteCount][stateCount];
    }
    
    public void setStates(int nodeId, int site, double[] states) {
        for (int i = 0; i < pp[0][0].length; i++) {
            pp[nodeId][site][i]=states[i];
        }
    }
    public void setState(int nodeId, int site, int stateIndex, double ppValue ) {
        pp[nodeId][site][stateIndex]=ppValue;
    }    
    public double getPP(int nodeId, int site, int state) {
        return pp[nodeId][site][state];
    }
    
    public double[][] getPPSet(int nodeId, int siteStart, int siteEnd) {
        if (siteStart>siteEnd) {return null;}
        return Arrays.copyOfRange(pp[nodeId], siteStart, siteEnd+1);
    }
    
    public int getNodeCount() {
        return pp.length;
    }
    
    public int getSiteCount() {
        return pp[0].length;
    }
    
    
    /**
     * flatten the probas by applying a rounding to 1e-(rank)
     * @param rank 
     */
    public void simplifyByRounding (int rank) {
        
    }
    
}
