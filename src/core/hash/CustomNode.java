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
    ArrayList<PairList> positionsPointers=new ArrayList<>(3);

    public CustomNode() {}

    public void registerTuple(int nodeId, int refPosition, float PPStar) {
        boolean refPositionAlreadyRegistered=false;
        for (PairList p:positionsPointers) {
            if (p.getRefPosition()==refPosition) {
                refPositionAlreadyRegistered=true;
                break;
            }
        }
        if (!refPositionAlreadyRegistered) {
            PairList pl=new PairList(refPosition);
            pl.add(new Pair(nodeId, PPStar));
            positionsPointers.add(pl);
        } else {
            getPairList(refPosition).add(new Pair(nodeId, PPStar));
        }
    }

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
        for (Iterator<PairList> iterator = positionsPointers.iterator(); iterator.hasNext();) {
            PairList ppl = iterator.next();
            Collections.sort(ppl);
        }
        //then sort these list (= sort the position by the PP* at 1st 
        //position of the list
        Collections.sort(positionsPointers);            
    }  


    /**
     * internal class just to link a reference position to a LinkedList
     */
    private class PairList extends ArrayList<Pair> implements Comparable<PairList>,Serializable {
        
        private static final long serialVersionUID = 7210L;
        
        int refPosition=-1;

        public PairList(int refPosition) {
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
        public int compareTo(PairList o) {
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
