/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import alignement.Alignment;
import core.DNAStates;
import core.PProbasSorted;
import core.ProbabilisticWord;
import core.States;
import core.hash.CustomHash_v2;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.PAMLWrapper;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class WordExplorer {

    //set the array to the number of words that will be tested
    //initial capacity of the arraylist is arbitrary, but tests showed that
    //only ~ 10000 words survived the 2 tresholds if both set at 1e-6
    //table realloc will occur if we get behind this size,
    //I was initally using testedWordCount as the maximum size, but this is
    //a memory killer for nothing...
    int initalCapacity=1000;
    ArrayList<ProbabilisticWord> words=new ArrayList(initalCapacity);
     
    //external parameters
    int k=-1;
    int refPosition=-1;
    int nodeId=-1;
    float wordThresholdAsLog=0.0f;
    
    //table to register proba proportion that is left on 
    float[] probaLeftPerCol=null;
    
    //variables for the recursive algorithm
    float currentLogSum=0.0f;
    byte[] word=null;
    PProbasSorted ppSet=null;
    boolean boundReached=false;
    boolean wordCompression=false;
    States s=null;


    /**
     *
     * @param k
     * @param refPosition
     * @param nodeId
     * @param ppSet
     * @param wordThresholdAsLog
     */
    public WordExplorer(int k, int refPosition, int nodeId, PProbasSorted ppSet, float wordThresholdAsLog, boolean wordCompression, States s) {
        this.refPosition=refPosition;
        this.wordThresholdAsLog=wordThresholdAsLog;
        this.word=new byte[k];
        this.k=k;
        this.nodeId=nodeId;
        this.ppSet=ppSet;
        this.wordCompression=wordCompression;
        this.s=s;
    }
    
    public ArrayList<ProbabilisticWord> getRetainedWords() {
        return words;
    }
    
    /**
     * 
     * @param i
     * @param j
     */
    public void exploreWords(int i, int j) {
        //for the algorithm below, consider ppSet as a (i,j) matrix
        //- with i the position on the ref alignment, refPosition being the i=0
        //  and refPosition+k being the maximumi, with i=k-1.
        //- with j the current state, descending ordered by their PP, from j=0
        //  to a maxumim j=(#states)
        //System.out.println("IN: "+i+" "+j);
        word[i-refPosition]=ppSet.getState(nodeId, i, j);
        //System.out.println("sumCurrentWord="+currentLogSum+"+"+ppSet.getPP(nodeId, i, j));
        currentLogSum+=ppSet.getPP(nodeId, i, j);
        boundReached = currentLogSum<wordThresholdAsLog;
        
        //register word if k-th position
        if (i==(refPosition+k-1)) {
            //register word
            if (!boundReached) {
                if (wordCompression) {
                    words.add(new ProbabilisticWord(s.compressMer(Arrays.copyOf(word, word.length)), currentLogSum, refPosition ));
                } else {
                    words.add(new ProbabilisticWord(Arrays.copyOf(word, word.length), currentLogSum, refPosition ));
                }
                //System.out.println("REGISTER: "+Arrays.toString(word)+" log10(PP*)="+currentLogSum);
            }
            //decrease before return
            //System.out.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
            currentLogSum-=ppSet.getPP(nodeId, i, j);
            //force return to parent
            return;
        } else {
            //go down recursively
            for (int j2 = 0; j2 < ppSet.getStateCount(); j2++) {
                if (boundReached) {break;}
                exploreWords(i+1, j2);
                //System.out.println("OUT: "+(i+1)+" "+j2);
            }
        }
        //decrease before natural return
        //System.out.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
        currentLogSum-=ppSet.getPP(nodeId, i, j);
        
    }
    
    
    
    
}
