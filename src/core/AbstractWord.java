/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import core.algos.JenkinsHash;
import java.io.Serializable;
import java.util.Arrays;

/**
 *
 * @author ben
 */
public abstract class AbstractWord implements Word,Serializable {
    
    private static final long serialVersionUID = 4020L;
    
    protected byte[] word;

    @Override
    public byte[] getWord() {
        return word;
    }

    /**
     * Note overriding hascode() obliges to override equals()\n
     * If note, many functions of the JDK will behave inconsistently
     * @param obj
     * @return 
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null) return false;
        if (obj == this) return true;
        if (!(obj instanceof Word))return false;
        Word otherObj = (Word)obj;
        return Arrays.equals(this.word,otherObj.getWord());
    }

    /**
     * Note overriding hascode() obliges to override equals()\n
     * If note, many functions of the JDK will behave inconsistently
     * @return 
     */
    @Override
    public int hashCode() {
        return JenkinsHash.hash32(this.word, 1);
    }


    
}
