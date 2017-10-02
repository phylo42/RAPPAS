/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import java.io.Serializable;


/**
 * Pair_16_32_bit representing a (nodeId,PP*) association as (unsigned int, float).
 * These are used in the LinkedLists
 * @author ben
 */
public class Pair_16_32_bit implements Pair,Serializable {
    
    private static final long serialVersionUID = 7101L;

    //as java as no unsigned int, we simply use char, to store the nodeIds
    //this allows a max of  2^16-1 node ids and is 16bit instead of the 32bits of int
    private char nodeId='\u0000'; //init to 0
    private float PPStar=Float.NEGATIVE_INFINITY;

    public Pair_16_32_bit(int nodeId, float PPStar) {
        if (nodeId>='\uFFFF') {
            System.out.println("Extended tree node number reached unsigned short limit (2^16-1). RAPPAS currently not designed to handle such big trees.");
            System.exit(1);
        }
        this.nodeId=(char)nodeId;
        this.PPStar=PPStar;
    }

    @Override
    public int getNodeId() {
        return (int)nodeId;
    }

    @Override
    public float getPPStar() {
        return PPStar;
    }

    @Override
    public void setNodeId(int nodeId) {
        if (nodeId>='\uFFFF') {
            System.out.println("Extended tree node number reached unsigned short limit (2^16-1). RAPPAS currently not designed to handle such big trees.");
            System.exit(1);
        }
        this.nodeId = (char)nodeId;
    }

    @Override
    public void setPPStar(float PPStar) {
        this.PPStar = PPStar;
    }


    @Override
    public String toString() {
        return "nodeId="+(int)nodeId+" PPStar="+PPStar;
    }
    
    

}
