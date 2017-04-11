/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

/**
 * object to represent the scoring of the nodes
 * @author ben
 */
public class Score implements Comparable<Score> {

    private int nodeId;
    private float PPStar;

    public Score(int nodeId, float PPStar) {
        this.nodeId = nodeId;
        this.PPStar = PPStar;
    }

    public int getNodeId() {
        return nodeId;
    }

    public float getPPStar() {
        return PPStar;
    }

    @Override
    public int compareTo(Score o) {
        if (this.PPStar<o.PPStar) {
            return 1;
        } else if (this.PPStar>o.PPStar) {
            return -1;
        } else {
            return 0;
        }
    }
    
    
}
