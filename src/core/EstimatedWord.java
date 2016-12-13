/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.io.Serializable;
import java.util.AbstractMap;
import java.util.Arrays;
import java.util.Comparator;

/**
 *
 * @author ben
 */
public class EstimatedWord extends AbstractWord implements Comparable<EstimatedWord>,Serializable  {
        
    private static final long serialVersionUID = 4021L;
    
    private double ppValue;
    
    public EstimatedWord(byte[] key, Double value) {
        this.word=key;
        this.ppValue=value;
    }

    /**
     * the comparator is inversed to put highest values first
     * @param o
     * @return 
     */
    public int compareTo(EstimatedWord o) {
        if ((this.getPpValue()-o.getPpValue())<0.0) {
            return 1;
        } else if ((this.getPpValue()-o.getPpValue())>0.0){
            return -1;
        } else {
            return 0;
        }
    }

    public double getPpValue() {
        return ppValue;
    }

    
    @Override
    public String toString() {
        return Arrays.toString(word)+" "+ppValue;
    }

    public String toStringNice(States s) {
        char[] wordChar=new char[word.length];
        for (int i = 0; i < wordChar.length; i++) {
            wordChar[i]=s.byteToState(word[i]);
        }
        
        return Arrays.toString(word)+" "+Arrays.toString(wordChar)+" "+ppValue;
    }


    
}