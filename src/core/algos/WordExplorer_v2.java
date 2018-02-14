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
import inputs.PHYMLWrapper;
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
 * second version of word explorer, which takes into account empty columns to generate k-mers
 * the option 'onlyOneJump', restrain the kmer generation to a single jump
 * when several jumps are possible inside the kmer
 * @author ben
 */
public class WordExplorer_v2 {

    //set the array to the number of words that will be tested
    //initial capacity of the arraylist is arbitrary, but tests showed that
    //only ~ 10000 words survived the 2 tresholds if both set at 1e-6
    //table realloc will occur if we get behind this size,
    //I was initally using testedWordCount as the maximum size, but this is
    //a memory killer for nothing...
    int initalCapacity=1000;
    ArrayList<ProbabilisticWord> words=new ArrayList(initalCapacity);
     
    //external variables used in the recursion
    int k=-1;
    int refPosition=-1;
    int current_k=-1;
    int nodeId=-1;
    float wordThresholdAsLog=0.0f;
    int idxOfFirstJump=-1;
    
    //table to register proba proportion that is left on 
    float[] probaLeftPerCol=null;
    
    //variables for the recursive algorithm
    float currentLogSum=0.0f;
    byte[] word=null;
    PProbasSorted ppSet=null;
    boolean boundReached=false;
    boolean wordCompression=false;
    States s=null;
    //represents  gap intervals
    boolean doGapJumps=false;
    ArrayList<Integer>[] gapIntervals =null;
    boolean limitTo1Jump=false;

    /**
     *
     * @param k
     * @param refPosition
     * @param nodeId
     * @param ppSet
     * @param align
     * @param wordThresholdAsLog
     * @param wordCompression
     * @param s
     * @param doGapJumps
     * @param limitTo1Jump
     */
    public WordExplorer_v2( int k,
                            int refPosition,
                            int nodeId,
                            PProbasSorted ppSet,
                            Alignment align,
                            float wordThresholdAsLog,
                            boolean wordCompression,
                            States s,
                            boolean doGapJumps,
                            boolean limitTo1Jump) {
        this.refPosition=refPosition;
        this.wordThresholdAsLog=wordThresholdAsLog;
        this.word=new byte[k];
        this.k=k;
        this.nodeId=nodeId;
        this.ppSet=ppSet;
        this.wordCompression=wordCompression;
        this.s=s;
        this.current_k=0;
        this.gapIntervals = align.getGapIntervals();
        this.doGapJumps=doGapJumps;
        this.limitTo1Jump=limitTo1Jump;
        
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
        //- with i the position on the ref alignment, with refPosition=i at start
        //  and refPosition+k being the maximum where i=k-1.
        //- with j the current state, in descending ordered by their PP, from j=0
        //  to a maxumim j=(#states)
        

        //Infos.println("ENTERING : "+i+" "+j+" , current_k:"+current_k);
        //if after the alignment limit (can happen because of gap jumps)
        //then return, assigning a value to this kmer is impossible.
        if (i>ppSet.getSiteCount()-1) {
            return;
        }
        //reset the jump counter if we come back to 1st mer position
        if (current_k==0) {
            idxOfFirstJump=-1;
        }
        
        word[current_k]=ppSet.getState(nodeId, i, j);
        //Infos.println("sumCurrentWord="+currentLogSum+"+"+ppSet.getPP(nodeId, i, j));
        currentLogSum+=ppSet.getPP(nodeId, i, j);
        boundReached = currentLogSum<wordThresholdAsLog;
        
        //register word if k-th position
        if (current_k==k-1) {
            //register word
            if (!boundReached) {
                if (wordCompression) {
                    words.add(new ProbabilisticWord(s.compressMer(Arrays.copyOf(word, word.length)), currentLogSum, refPosition ));
                } else {
                    words.add(new ProbabilisticWord(Arrays.copyOf(word, word.length), currentLogSum, refPosition ));
                }
                //Infos.println("REGISTER: "+Arrays.toString(word)+" log10(PP*)="+currentLogSum);
            }
            //decrease before return
            //Infos.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
            currentLogSum-=ppSet.getPP(nodeId, i, j);
            //force return to parent
            return;
        } else {
            
            //go down recursively, testing each state at position i+1
            for (int j2 = 0; j2 < ppSet.getStateCount(); j2++) {
                if (boundReached) {break;}
                
                //do this exploration, do not consider gap interval
                //Infos.println("EXPORATION NORMAL");
                //we move to next kmer position
                current_k++;
                exploreWords(i+1, j2);
                current_k--;
                //Infos.println("OUT from : "+(i+1)+" "+j2+" , current_k:"+current_k+" , currently HERE: "+i+" "+j);
                //redo exploration, by jumping over the gap intervals if any
                //do condition on i because i+1 cannot be outside the alignment
                if (doGapJumps && i<ppSet.getSiteCount()-1) {
                    //verify if next i is associated to a gap interval
                    if (gapIntervals[i+1]!=null) { 
                        //do all jumps combinations 
                        if (!limitTo1Jump) {
                            for (int i_interval = 0; i_interval < gapIntervals[i+1].size(); i_interval++) {
                                //we move to next kmer position and explore
                                current_k++;
                                exploreWords((i+1)+gapIntervals[i+1].get(i_interval), j2); //jump from i to i+1+gapInterval_length
                                current_k--;
                            }
                        //or allow only 1 jump
                        } else {
                            if ((idxOfFirstJump==-1) ) { //do a jump only if not previously done
                                idxOfFirstJump=i;
                                for (int i_interval = 0; i_interval < gapIntervals[i+1].size(); i_interval++) {
                                    //Infos.println("EXPLORATION WITH GAP JUMP: gap_i+1="+(i+1)+" length="+gapIntervals[i+1].get(i_interval));
                                    //Infos.println("NEXT i jump to --> "+((i+1)+gapIntervals[i+1].get(i_interval)));
                                    //we move to next kmer position
                                    current_k++;
                                    exploreWords((i+1)+gapIntervals[i+1].get(i_interval), j2); //jump from i to i+1+gapInterval_length
                                    current_k--;
                                    //Infos.println("OUT from : "+(i+gapIntervals[i+1].get(i_interval))+" "+j2+" , current_k:"+current_k+" , currently HERE: "+i+" "+j);
                                }
                            } else {
                                //Infos.println("DO NOT JUMP, already jumped at idxOfFirstJump="+idxOfFirstJump);
                            }
                        }
                    }
                }
            }
            
            
            
        }
        //decrease before natural return
        //Infos.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
        currentLogSum-=ppSet.getPP(nodeId, i, j);
    }
    
    
    
}
