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
    Random rand = new Random();
    //standard dev arounf the mean read length
    int Rsd=5;
    
    public RandomSeqGenerator(States s) {
        this.s=s;
    }
    
    public List<Fasta> generateSequence(int amount, int sequenceLength) {
        ArrayList<Fasta> list=new ArrayList<>(amount);
        for (int i = 0; i < amount; i++) {
            //select normally distibuted read length v_length
            //centered around R and with standard dev of value Rsd
            double Rsd=sequenceLength*0.1;
            int v_length = new Double(0.0+sequenceLength+rand.nextGaussian()*Rsd).intValue();
            StringBuilder sb = new StringBuilder();
            rand.ints(v_length, 0, s.getNonAmbiguousStatesCount()).boxed().map(intVal->s.byteToState(intVal.byteValue())).forEachOrdered(sb::appendCodePoint);
            list.add(new Fasta("rand"+i, sb.toString()));
            //System.out.println(f.getFormatedFasta());
        }
        return list;
    }
    
    
    public static void main(String[] args) {
        RandomSeqGenerator rs=new RandomSeqGenerator(new DNAStates());
        rs.generateSequence(10, 100);
    }
    
}
