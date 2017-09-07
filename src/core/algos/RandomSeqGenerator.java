/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import core.DNAStates;
import core.States;
import inputs.Fasta;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * generate random sequences
 * @author ben
 */
public class RandomSeqGenerator {
    
    States s=null;
    //pseudorandom generator
    Long seed=new Long(1);
    Random rand = new Random(seed);
    int sequenceLength=0;
    //standard dev arounf the mean read length
    double Rsd=5;
    Double v_length = 0.0;
    
    //variables used for Fasta creation
    StringBuffer sb=null;
    StringBuffer header=null;
    Fasta fasta=null;
    
    //int counter
    int counter=-1;
    
    public RandomSeqGenerator(States s,int sequenceLength) {
        this.s=s;
        this.sb= new StringBuffer(sequenceLength+new Double(4*Rsd).intValue());
        this.header=new StringBuffer(4+String.valueOf(sequenceLength).length());
        this.sequenceLength=sequenceLength;
        this.Rsd=sequenceLength*0.1;
    }
    
    public Fasta generateSequence() {
        sb.delete(0, sb.length());
        //select normally distibuted read length v_length
        //centered around R and with standard dev of value Rsd
        v_length = 0.0+sequenceLength+rand.nextGaussian()*Rsd;
        for (int i = 0; i < v_length; i++) {
            sb.appendCodePoint(s.byteToState((byte)(rand.nextInt(s.getNonAmbiguousStatesCount()))));            
        }
        counter++;
        return new Fasta("rand"+counter, sb.toString());         
    }
    
    
    public static void main(String[] args) {
        RandomSeqGenerator rs=new RandomSeqGenerator(new DNAStates(),20);
        for (int i = 0; i < 10; i++) {
            Fasta f=rs.generateSequence();
            System.out.println(f.getFormatedFasta());
        }
    }
    
}
