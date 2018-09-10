/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import java.io.Serializable;
import java.util.Collections;
import java.util.List;


/**
 * Intermediate table used to select the (nodeId,PP*) pairs through 
 * associated reference position.
 * @author ben
 */
public class UnionPointer implements HashPointer,Serializable {
    
    private static final long serialVersionUID = 7300L;

    ObjectList<Pair> list=null;

    /**
     *
     */
    public UnionPointer() {
        list=new ObjectArrayList<>();
    }

    @Override
    public void registerTuple(int nodeId, int refPosition, float PPStar) {
        //No pairs associatd to the word yet
//        if (list.size()>0) {
//            OptionalInt finding =   IntStream.rangeClosed(0, list.size()-1)
//                                    .filter(i->list.get(i).getNodeId()==nodeId)
//                                    .findFirst();
//            if (finding.isPresent()) {
//                int i=finding.getAsInt();
//                //System.out.println("In node: Pair of same nodeId already existing at index: "+i);
//                //replace previous pair if better PP*
//                if (PPStar > list.get(i).getPPStar()) {
//                    list.get(i).setPPStar(PPStar);
//                }   
//            } else {
//                //System.out.println("In node: New pair created: ");
//                list.add(new Pair_16_32_bit(nodeId, PPStar));
//            }
//        } else {
//            //System.out.println("In node: First pair created: ");
//            //assign new pair if node not already registered
//            list.add(new Pair_16_32_bit(nodeId, PPStar));
//        }        
        if (list.size()>0) {
            boolean found=false;
            for (Pair p:list) {
                //nodeid already registerd
                if (p.getNodeId()==nodeId) {
                    found=true;
                    //replace previous pair if better PP*
                    if (PPStar > p.getPPStar()) {
                        p.setPPStar(PPStar);
                    }
                    break;
                }
            }
            //node not yet registered, do it
            if (!found) {
                list.add(new Pair_16_32_bit(nodeId, PPStar));
            }
        } else {
            //create very first pair
            list.add(new Pair_16_32_bit(nodeId, PPStar));
        }

    }
    
    /**
     * ref alignment positions associated to this node, order by underlying max(PP*)
     * @return 
     */
    @Override
    public int[] getPositions() {
        int[] t={0};
        return t;
    }

    /**
     * get all Pairs associated to a particular reference position.
     * @param refPosition
     * @return 
     */
    @Override
    public List<Pair> getPairList(int refPosition) {
        return list;
    }
    
    /**
     * get number of nodes associated to best position
     * @return 
     */
    @Override
    public int getPairCountInTopPosition() {
        return list.size();
    }

    /**
     * get best (nodeId;PP*).
     * @return 
     */
    @Override
    public Pair getBestPair() {
        return list.get(0);
    }

    /**
     * get best reference position.
     * @return 
     */
    @Override
    public int getBestPosition() {
        return 0;
    }

    @Override
    public void sort() {
        Collections.sort(list);            
    }  
}
