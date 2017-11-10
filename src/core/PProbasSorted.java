/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author ben
 */
public class PProbasSorted implements Serializable {
    
    
    float[][][] pp=null; //[nodeIndex][site][index]=PPStar
    byte[][][] states=null;//[nodeIndex][site][index]= State associated to PPStar
    //example on pplacer benchmark:
    //pp[O][0]=[-0.51749474, -0.5934179, -0.64727, -0.665648]
    //states[0][0]=[3, 0, 2, 1]
    //statesIndex[0][0]=[1, 3, 2, 0]
    //note that the index of a node is defined through its DFS pre_order

    public PProbasSorted(int nodeCount, int siteCount, int stateCount) {
        pp=new float[nodeCount][siteCount][stateCount];
        states=new byte[nodeCount][siteCount][stateCount];
    }
    
    public void setStates(int nodeId, int site, ArrayList<SiteProba> probas) {

        assert probas.size()==pp[0][0].length;
        
        byte i=0;
        for (SiteProba sp:probas) {
            pp[nodeId][site][i]=sp.proba;
            states[nodeId][site][i]=sp.state;
            i++;
        }
    }
      
    
    public double getPP(int nodeId, int site, int index) {
        return pp[nodeId][site][index];
    }
    
    public byte getState(int nodeId, int site, int index) {
        return states[nodeId][site][index];
    }
    
    public byte getStateIndex(int nodeId, int site, byte state) {
        return states[nodeId][site][state];
    }
    

    public float[][] getPPSet(int nodeId, int siteStart, int siteEnd) throws IndexOutOfBoundsException {
        if (siteEnd>=pp[0].length) {
            throw new IndexOutOfBoundsException("siteEnd="+siteEnd+" > to last site position");
        } else {
            return Arrays.copyOfRange(pp[nodeId], siteStart, siteEnd+1);
        }
    }
    
    public byte[][] getStateSet(int nodeId, int siteStart, int siteEnd) {
        if (siteEnd>=pp[0].length) {
            throw new IndexOutOfBoundsException("siteEnd="+siteEnd+" > to last site position");
        } else {
            return Arrays.copyOfRange(states[nodeId], siteStart, siteEnd+1);
        }
    }
    
    public byte[][] getStateIndexSet(int nodeId, int siteStart, int siteEnd) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public int getNodeCount() {
        return pp.length;
    }
    
    public int getSiteCount() {
        return pp[0].length;
    }
    
    public int getStateCount() {
        return pp[0][0].length;
    }
    
    
    
}
