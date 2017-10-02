/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import java.io.Serializable;


/**
 * Pair_32_32_bit representing a (nodeId,PP*) association as (int,float).
 * These are used in the LinkedLists
 * @author ben
 */
public class Pair_32_32_bit implements Pair,Serializable {
    
    private static final long serialVersionUID = 7102L;

    private int nodeId=0; //init to 0
    private float PPStar=Float.NEGATIVE_INFINITY;

    public Pair_32_32_bit(int nodeId, float PPStar) {
        this.nodeId=nodeId;
        this.PPStar=PPStar;
    }

    @Override
    public int getNodeId() {
        return nodeId;
    }

    @Override
    public float getPPStar() {
        return PPStar;
    }

    @Override
    public void setNodeId(int nodeId) {
        this.nodeId = nodeId;
    }

    @Override
    public void setPPStar(float PPStar) {
        this.PPStar = PPStar;
    }

    /**
     * the comparator is inversed to put highest values first
     * @param o
     * @return 
     */
    public int compareTo(Pair_32_32_bit o) {
        return -Float.compare(this.PPStar,o.getPPStar());
    }

    @Override
    public String toString() {
        return "nodeId="+nodeId+" PPStar="+PPStar;
    }
    
    

}
