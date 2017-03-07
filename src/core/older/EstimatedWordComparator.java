/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.older;

import core.ProbabilisticWord;
import java.util.Comparator;

/**
 *
 * @author ben
 */
public class EstimatedWordComparator implements Comparator<ProbabilisticWord> {
    @Override
    public int compare(ProbabilisticWord o1, ProbabilisticWord o2) {
        if ((o1.getPpStarValue()-o2.getPpStarValue())>0.0) {
            return 1;
        } else if ((o1.getPpStarValue()-o2.getPpStarValue())<0.0){
            return -1;
        } else {
            return 0;
        }
    } 
}
