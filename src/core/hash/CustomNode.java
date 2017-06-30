/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;


/**
 * Intermediate table used to select the (nodeId,PP*) pairs through 
 * associated reference position.
 * @author ben
 */
public class CustomNode implements Serializable {
    
    private static final long serialVersionUID = 7200L;


    //small init capacity because if k around 8 to 12, few occurences are
    //expected in the ref alignment
    private ArrayList<PositionPointer> positionsPointers=new ArrayList<>(3);

    public CustomNode() {}

    public void registerTuple(int nodeId, int refPosition, float PPStar) {
        boolean refPositionAlreadyRegistered=false;
        for (PositionPointer p:positionsPointers) {
            if (p.getRefPosition()==refPosition) {
                refPositionAlreadyRegistered=true;
                break;
            }
        }
        if (!refPositionAlreadyRegistered) {
            PositionPointer pl=new PositionPointer(refPosition);
            pl.add(new Pair(nodeId, PPStar));
            positionsPointers.add(pl);
        } else {
            getPairList(refPosition).add(new Pair(nodeId, PPStar));
        }
    }

    /**
     * ref alignment positions associated to this node, order by underlying max(PP*)
     * @return 
     */
    public int[] getPositions() {
        return positionsPointers.stream().mapToInt(pp->pp.refPosition).toArray();
    }

    /**
     * get all Pairs associated to a particular reference position.
     * @param refPosition
     * @return 
     */
    public ArrayList<Pair> getPairList(int refPosition) {
        for (int i=0;i<positionsPointers.size();i++) {
            if (positionsPointers.get(i).getRefPosition()==refPosition)
                return positionsPointers.get(i);
        }
        return null;
    }

    /**
     * get best (nodeId;PP*).
     * @return 
     */
    public Pair getBestPair() {
        return positionsPointers.get(0).get(0);
    }

    /**
     * get best reference position.
     * @return 
     */
    public Integer getBestPosition() {
        return positionsPointers.get(0).getRefPosition();
    }

    public void sort() {
        //sort each list of Pair objects
        for (Iterator<PositionPointer> iterator = positionsPointers.iterator(); iterator.hasNext();) {
            PositionPointer ppl = iterator.next();
            Collections.sort(ppl);
        }
        //then sort these list (= sort the position by the PP* at 1st 
        //position of the list
        Collections.sort(positionsPointers);            
    }  

    
    /**
     * from the "full" hash, remove all (nodeId,PP*) pairs associated to positions
     * which are not the position holding the best PP* --> "best position" hash
     */
    public void clearPairsOfWorsePositions() {
        for (int i = 1; i < positionsPointers.size(); i++) {
            positionsPointers.get(i).clear();
        }
    }
    
    /**
     * from the "full" or "best position" hash, remove all (nodeId,PP*) 
     * pairs associated to positions which are not the position holding
     * the best PP*, then keeps only X (nodeid,PP*) pairs
     */
    public void clearPairsOfWorsePositionsAndLimitToXPairs(int X) {
        clearPairsOfWorsePositions();
        if (positionsPointers.get(0).size()>X) {
            for (int i = X; i < positionsPointers.get(0).size(); i++) {
                positionsPointers.get(0).remove(i);
            }
        }
    }
        
    

    /**
     * internal class just to link a reference position to a LinkedList
     */
    private class PositionPointer extends ArrayList<Pair> implements Comparable<PositionPointer>,Serializable {
        
        private static final long serialVersionUID = 7210L;
        
        private int refPosition=-1;

        public PositionPointer(int refPosition) {
            this.refPosition=refPosition;
        }

        public int getRefPosition() {
            return refPosition;
        }

        /**
         * Comparable retrieving inversed order to set. Work considering 
         * that the list of Pair object is already sorted (best PP* as
         * index 0).
         * @param o
         * @return 
         */
        public int compareTo(PositionPointer o) {
            if ((this.get(0).getPPStar()-o.get(0).getPPStar())<0.0) {
                return 1;
            } else if ((this.get(0).getPPStar()-o.get(0).getPPStar())>0.0){
                return -1;
            } else {
                return 0;
            }
        }


    }


}
