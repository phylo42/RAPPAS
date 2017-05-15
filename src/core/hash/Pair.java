/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import java.io.Serializable;


/**
 * Pair representing a (nodeId,PP*) association. These are used in the 
 * LinkedLists
 * @author ben
 */
public class Pair implements Comparable<Pair>,Serializable {
    
    private static final long serialVersionUID = 7100L;

    
    int nodeId=-1;
    float PPStar=-1.0f;

    public Pair(int nodeId, float PPStar) {
        this.nodeId=nodeId;
        this.PPStar=PPStar;
    }

    public int getNodeId() {
        return nodeId;
    }

    public float getPPStar() {
        return PPStar;
    }

    public void setNodeId(int nodeId) {
        this.nodeId = nodeId;
    }

    public void setPPStar(float PPStar) {
        this.PPStar = PPStar;
    }

    @Override
    /**
     * the comparator is inversed to put highest values first
     * @param o
     * @return 
     */
    public int compareTo(Pair o) {
        if ((this.PPStar-o.getPPStar())<0.0) {
            return 1;
        } else if ((this.PPStar-o.getPPStar())>0.0){
            return -1;
        } else {
            return 0;
        }
    }

    @Override
    public String toString() {
        return "nodeId="+nodeId+" PPStar="+PPStar;
    }
    
    

}
