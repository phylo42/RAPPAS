/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author ben
 */
public interface HashPointer {
    
    /**
     * register a tuple into the hash
     * @param nodeId
     * @param refPosition
     * @param PPStar 
     */
    public void registerTuple(int nodeId, int refPosition, float PPStar);
    
    /**
     * ref alignment positions associated to this node, order by underlying max(PP*)
     * @return 
     */
    public int[] getPositions();

    /**
     * get all Pairs associated to a particular reference position.
     * @param refPosition
     * @return 
     */
    public List<Pair> getPairList(int refPosition);
    
    /**
     * get number of nodes associated to best position
     * @return 
     */
    public int getPairCountInTopPosition();
    
    /**
     * get best (nodeId;PP*).
     * @return 
     */
    public Pair getBestPair();
    
    /**
     * get best reference position.
     * @return 
     */
    public int getBestPosition();
    
    /**
     * sort positions pointer (through their 1st associated Pair_16_32_bit PP*), then
 sort the Pair_16_32_bit list itself
     */
    public void sort();

}
