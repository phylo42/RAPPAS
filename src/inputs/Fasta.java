/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import core.Word;
import core.algos.JenkinsHash;
import java.util.Arrays;

/**
 * Fasta representation, not that the '' is not in the header
 * @author ben
 */
public class Fasta {
    
    String sequence=null;
    String header=null;

    public Fasta(String header, String sequence) {
        this.sequence=sequence;
        this.header=header;
    }

    public String getHeader() {
        return header;
    }

    /**
     * the sequence itself, without trailing characters
     * @param removeGaps 
     */
    public String getSequence(boolean removeGaps) {
        if (!removeGaps)
            return sequence;
        else 
            return sequence.replaceAll("-", "");
    }

    public void setHeader(String header) {
        this.header = header;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
    
    /**
     * returns >header\nsequence
     * @return 
     */
    public String getFormatedFasta() {
        return ">"+header+"\n"+sequence;
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
        if (!(obj instanceof Fasta))return false;
        Fasta otherObj = (Fasta)obj;
        return this.sequence.equals(otherObj.getSequence(false));
    }

    /**
     * Note overriding hascode() obliges to override equals()\n
     * If note, many functions of the JDK will behave inconsistently
     * @return 
     */
    @Override
    public int hashCode() {
        return JenkinsHash.hash32(this.sequence.getBytes(), 1);
    }

    @Override
    public String toString() {
        return header+"(l="+sequence.length()+")";
    }
    
    
    
}
