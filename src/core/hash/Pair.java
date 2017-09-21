/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

/**
 * A (nodeId,PP*) association
 * @author ben
 */
public interface Pair extends Comparable<Pair> {
    public int getNodeId();
    public float getPPStar();
    public void setNodeId(int nodeId);
    public void setPPStar(float PPStar);

    @Override
    public default int compareTo(Pair o) {
        return -Float.compare(this.getPPStar(),o.getPPStar());
    }
   

}
