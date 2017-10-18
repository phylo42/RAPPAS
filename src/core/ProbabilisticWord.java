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
public class ProbabilisticWord extends AbstractWord implements Comparable<ProbabilisticWord>,Serializable  {
        
    private static final long serialVersionUID = 4021L;
    
    private float ppStarValue;
    private int originalPosition;
    
    public ProbabilisticWord(byte[] key, float value) {
        this.word=key;
        this.ppStarValue=value;
    }
    
    public ProbabilisticWord(byte[] key, float value, int originalPosition) {
        this.word=key;
        this.ppStarValue=value;
        this.originalPosition=originalPosition;
    }

    /**
     * the comparator is inversed to put highest values first
     * @param o
     * @return 
     */
    public int compareTo(ProbabilisticWord o) {
        if ((this.getPpStarValue()-o.getPpStarValue())<0.0) {
            return 1;
        } else if ((this.getPpStarValue()-o.getPpStarValue())>0.0){
            return -1;
        } else {
            return 0;
        }
    }

    public float getPpStarValue() {
        return ppStarValue;
    }

    public int getOriginalPosition() {
        return originalPosition;
    }

    public void update(byte[] key, float value, int originalPosition) {
        this.word=key;
        this.ppStarValue=value;
        this.originalPosition=originalPosition;
    }
    
    
    @Override
    public String toString() {
        return Arrays.toString(word)+" "+ppStarValue;
    }

    public String toStringNice(States s) {
        char[] wordChar=new char[word.length];
        for (int i = 0; i < wordChar.length; i++) {
            wordChar[i]=s.byteToState(word[i]);
        }
        
        return Arrays.toString(word)+" "+Arrays.toString(wordChar)+" "+ppStarValue;
    }


    
}