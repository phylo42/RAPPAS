/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.util.Comparator;

/**
 *
 * @author ben
 */
public class EstimatedWordComparator implements Comparator<EstimatedWord> {
    @Override
    public int compare(EstimatedWord o1, EstimatedWord o2) {
        if ((o1.getPpValue()-o2.getPpValue())>0.0) {
            return 1;
        } else if ((o1.getPpValue()-o2.getPpValue())<0.0){
            return -1;
        } else {
            return 0;
        }
    } 
}
