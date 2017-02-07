/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

/**
 *
 * @author ben
 */
public class SiteProba implements Comparable<SiteProba> {
    public float proba;
    public byte state;

    @Override
    /**
     * the comparator is inversed to put highest values first
     * @param o
     * @return 
     */
    public int compareTo(SiteProba o) {
        if ((this.proba-o.proba)<0.0) {
            return 1;
        } else if ((this.proba-o.proba)>0.0){
            return -1;
        } else {
            return 0;
        }
    }

    @Override
    public String toString() {
        return super.toString()+" SiteProba: "+state+":"+proba;
    }
    
    
}
