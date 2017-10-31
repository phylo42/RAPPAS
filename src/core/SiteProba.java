/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

/**
 * object used by Wrappers, during parsing of the PP. Simply represent a site
 * of the reference alignment and its associated PP, this can be considered 
 * as the tuple (proba,site), it is comparable and is compared through the 'proba'
 * elements of the tuple.
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
        return " SiteProba: state="+state+", PP="+proba;
    }
    
    
}
