/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;


/**
 * Intermediate table used to select the (nodeId,PP*) pairs through 
 * associated reference position.
 * @author ben
 */
public class PositionPointer implements HashPointer,Serializable {
    
    private static final long serialVersionUID = 7200L;

    //small init capacity because if k around 8 to 12, few occurences are
    //expected in the ref alignment
    private ArrayList<PositionList> positionPointerList=new ArrayList<>(3);

    /**
     *
     */
    public PositionPointer() {}

    @Override
    public void registerTuple(int nodeId, int refPosition, float PPStar) {
        boolean refPositionAlreadyRegistered=false;
        for (PositionList p:positionPointerList) {
            if (p.getRefPosition()==refPosition) {
                refPositionAlreadyRegistered=true;
                break;
            }
        }
        if (!refPositionAlreadyRegistered) {
            //create new position otherwise
            PositionList pl=new PositionList(refPosition);
            pl.add(new Pair_16_32_bit(nodeId, PPStar));
            positionPointerList.add(pl);
        } else {
            getPairList(refPosition).add(new Pair_16_32_bit(nodeId, PPStar));
        }
    }

    /**
     * ref alignment positions associated to this node, order by underlying max(PP*)
     * @return 
     */
    @Override
    public int[] getPositions() {
        int[] pos=new int[positionPointerList.size()];
        for (int i=0;i<positionPointerList.size();i++) {
            pos[i]=positionPointerList.get(i).getRefPosition();
        }
        return pos;
    }

    /**
     * get all Pairs associated to a particular reference position.
     * @param refPosition
     * @return 
     */
    @Override
    public ArrayList<Pair> getPairList(int refPosition) {
        //return positionPointerList.stream().filter(p->p.getRefPosition()==refPosition).findFirst().orElse(null);
        for (PositionList p:positionPointerList) {
            if (p.getRefPosition()==refPosition) {
                return p;
            }
        }
        return null;
    }
    
    /**
     * get number of nodes associated to best position
     * @return 
     */
    @Override
    public int getPairCountInTopPosition() {
        
        return this.positionPointerList.get(0).size();
    }
    
    /**
     * get best (nodeId;PP*).
     * @return 
     */
    @Override
    public Pair getBestPair() {
        return positionPointerList.get(0).get(0);
    }

    /**
     * get best reference position.
     * @return 
     */
    @Override
    public int getBestPosition() {
        return positionPointerList.get(0).getRefPosition();
    }

    /**
     * sort positions pointer (through their 1st associated Pair_16_32_bit PP*), then
 sort the Pair_16_32_bit list itself
     */
    @Override
    public void sort() {
        //sort each list of Pair_16_32_bit objects
        this.positionPointerList.stream().forEach(p->Collections.sort(p));
//        for (Iterator<PositionPointer> iterator = positionPointerList.iterator(); iterator.hasNext();) {
//            PositionList ppl = iterator.next();
//            Collections.sort(ppl);
//        }
        //then sort these list (= sort the position by the PP* at 1st 
        //position of the list)
        Collections.sort(positionPointerList);            
    }  

    
    /**
     * from the "full" hash, remove all (nodeId,PP*) pairs associated to positions
     * which are not the position holding the best PP* --> "best position" hash
     */
    public void clearPairsOfWorsePositions() {
        
        for (int i = positionPointerList.size()-1;i>0; i--) {
            positionPointerList.get(i).clear();
            positionPointerList.remove(i);
        }
    }
    
    /**
     * from the "full" or "medium" hash, remove all (nodeId,PP*) 
     * pairs associated to positions which are not the position holding
     * the best PP*, then keeps only X (nodeid,PP*) pairs
     * @param X
     * @return boolean that if true, notifies that this position should be deleted
     * because biased after the pair trimming
     */
    @Deprecated
    public boolean limitToXPairsPerPosition(int X) {
        if (positionPointerList.get(0).size()>X) {
            //do not forget X is shifted by -1 for array coordinates (idx9==10th elt)
            //remove nth to (X+2)th Pair_16_32_bit
            for (int i = positionPointerList.get(0).size()-1; i > X; i--) {
                positionPointerList.get(0).remove(i);
            }
            //now check the state of the [1st to (X+1)] pairs
            //(from right to left) memorize the index of the last PP* change
            // X.PP*==(X+1).PP*? T-> reachIndex0(no change)?  T->remove ALL
            //                                                F->remove [index;X+1]
            //                   F->remove [X+1;X+1]
            if (positionPointerList.get(0).get(X-1).getPPStar()
                !=positionPointerList.get(0).get(X).getPPStar()) {
                //no need to go further, the X Pairs are kept
                //just remove the X+1 pair
                positionPointerList.get(0).remove(X);
                return false;
            } else {
                for (int i = X-2; i > -1; i--) { 
                    //test value change between i and i+1, with i in [0,X-1]
                    if ((positionPointerList.get(0).get(i).getPPStar()-positionPointerList.get(0).get(i+1).getPPStar())!=0.0f) {
                        //rightmost change found, remove [rightmostchangeIndex;X+1]
                        for (int j = X; j > i; j--) {
                            positionPointerList.get(0).remove(j);
                        }
                        return false;
                    }
                }
                //no change found ! removes everything, this kmer is biased
                return true;
            }
        }
        return false;
    }
        
    

    /**
     * internal class just to link a reference position to a LinkedList
     */
    private class PositionList extends ArrayList<Pair> implements Comparable<PositionList>,Serializable {
        
        private static final long serialVersionUID = 7210L;
        
        private int refPosition=-1;

        public PositionList(int refPosition) {
            super(1);
            this.refPosition=refPosition;
        }

        public int getRefPosition() {
            return refPosition;
        }

        /**
         * Comparable retrieving inversed order to set. Work considering 
 that the list of Pair_16_32_bit object is already sorted (best PP* as
 index 0).
         * @param o
         * @return 
         */
        @Override
        public int compareTo(PositionList o) {
            return -Float.compare(this.get(0).getPPStar(), o.get(0).getPPStar());
        }

        @Override
        public boolean add(Pair e) {
            boolean ok= super.add(e);
            //this.trimToSize(); //force arrylist capacity to preserve a max of memory
            return ok;
        }
        
        


    }


}
