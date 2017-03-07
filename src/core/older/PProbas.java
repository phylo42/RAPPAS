/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.older;

import java.io.Serializable;
import java.util.Arrays;

/**
 *
 * @author ben
 */
public class PProbas implements Serializable {
    
    private static final long serialVersionUID = 3000L;
    
    float[][][] pp=null; //[nodeIndex][site][state]
    //note that the index of a node is defined through its DFS pre_order

    public PProbas(int nodeCount, int siteCount, int stateCount) {
        pp=new float[nodeCount][siteCount][stateCount];
    }
    
    public void setStates(int nodeId, int site, float[] states) {
        for (int i = 0; i < pp[0][0].length; i++) {
            pp[nodeId][site][i]=states[i];
        }
    }
    public void setState(int nodeId, int site, int stateIndex, float ppValue ) {
        pp[nodeId][site][stateIndex]=ppValue;
    }    
    public double getPP(int nodeId, int site, int state) {
        return pp[nodeId][site][state];
    }
    
    public float[][] getPPSet(int nodeId, int siteStart, int siteEnd) {
        if (siteStart>siteEnd) {return null;}
        return Arrays.copyOfRange(pp[nodeId], siteStart, siteEnd+1);
    }
    
    public int getNodeCount() {
        return pp.length;
    }
    
    public int getSiteCount() {
        return pp[0].length;
    }
    
    public float[][][] accessTable() {
        return pp;
    }
}
