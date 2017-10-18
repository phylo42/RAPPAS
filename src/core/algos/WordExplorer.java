/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import alignement.Alignment;
import core.DNAStates;
import core.older.PProbas;
import core.PProbasSorted;
import core.ProbabilisticWord;
import core.hash.CustomHash;
import core.States;
import core.Word;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.PAMLWrapper;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
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
    //only ~ 100-1000 words survived the 2 tresholds if both set at 1e-6
    //table realloc will occur if we get behind this size,
    //I was initally using testedWordCount as the maximum size, but this is
    //a memory killer ...
    int initalCapacity=100;
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
    
    //variables to manage words list 
    //i.e. to avoid too many memory footprint, this list need to be created at 
    //instanciation of WordExplorer, then the ProbabilisticWord created in the list
    //should not be replaced with new instanticiations but the primitive content
    //is updated instaed.
    //as the ProbabilisticWord list return by WordExplorer will have a different
    //size for different i,j positions (see exploreWords() ), we need variables
    //to determine which fraction of the list [0,n] should be returned
    int lastIndexOfCurrentRecursion=-1; //increases while exploreWords() finds new words
    boolean newCall=true; //to know when we are at first recursive call

    /**
     *
     * @param k
     * @param refPosition
     * @param nodeId
     * @param ppSet
     * @param wordThresholdAsLog
     */
    public WordExplorer(int k, int refPosition, int nodeId, PProbasSorted ppSet, float wordThresholdAsLog) {
        this.refPosition=refPosition;
        this.wordThresholdAsLog=wordThresholdAsLog;
        this.word=new byte[k];
        this.k=k;
        this.nodeId=nodeId;
        this.ppSet=ppSet;
    }
    
    public List<ProbabilisticWord> getRetainedWords() {
        if (lastIndexOfCurrentRecursion>-1) {
            return words.subList(0, lastIndexOfCurrentRecursion);
        } else {
            return new ArrayList<>(1); //empty list
        }
    }
    
    /**
     * return list of words retained for the k-mer determine by the i,j reference
     * alignment coordinates
     * @param i
     * @param j
     */
    public void exploreWords(int i, int j) {
        //for the algorithm below, consider ppSet as a (i,j) matrix
        //- with i the position on the ref alignment, refPosition being the i=0
        //  and refPosition+k being the maximumi, with i=k-1.
        //- with j the current state, descending ordered by their PP, from j=0
        //  to a maxumim j=(#states)
        
        //init index if this is a new call for these i,j interval
        if (newCall) {
            newCall=false;
            lastIndexOfCurrentRecursion=-1;
        }
        
        //System.out.println("IN: "+i+" "+j);
        word[i-refPosition]=ppSet.getState(nodeId, i, j);
        //System.out.println("sumCurrentWord="+currentLogSum+"+"+ppSet.getPP(nodeId, i, j));
        currentLogSum+=ppSet.getPP(nodeId, i, j);
        if(currentLogSum<wordThresholdAsLog)
            boundReached=true;
        else
            boundReached=false;
        
        //register word if k-th position
        if (i==(refPosition+k-1)) {
            //register word
            if (!boundReached) {
                //words list is too small to fit new word
                if (words.size()-1<lastIndexOfCurrentRecursion) {
                    lastIndexOfCurrentRecursion++;
                    words.add(new ProbabilisticWord(Arrays.copyOf(word, word.length), currentLogSum, refPosition ));
                //words list big enough, just update value in the word object    
                } else {
                    words.get(lastIndexOfCurrentRecursion).update(Arrays.copyOf(word, word.length), currentLogSum, refPosition);
                    lastIndexOfCurrentRecursion++;
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
    
    
    
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose", "1");
        
        try {
            
            
            String a="/media/ben/STOCK/SOURCES/NetBeansProjects/ViromePlacer/WD/relaxed_trees/relaxed_align_BrB_minbl0.001_1peredge.fasta";
            String rst="/media/ben/STOCK/SOURCES/NetBeansProjects/ViromePlacer/WD/AR/finished/rst";
            
            
            int k=8;
            float sitePPThreshold=Float.MIN_VALUE;
            int thresholdFactor=10;
            
            Infos.println("k="+k);
            Infos.println("factor="+thresholdFactor);
            float wordPPStarThreshold=(float)(thresholdFactor*Math.pow(0.25,k));
            Infos.println("wordPPStarThreshold="+wordPPStarThreshold);
            float thresholdAsLog=(float)Math.log10(wordPPStarThreshold);
            Infos.println("wordPPStarThreshold(log10)="+thresholdAsLog);
            
            //////////////////////
            //States: DNA or AA
            States s=new DNAStates();
            //////////////////////
            //LOAD ORIGINAL ALIGNMENT
            Infos.println("Loading Alignment...");
            FASTAPointer fp=new FASTAPointer(new File(a), false);
            Fasta fasta=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(fasta);
    }
            Alignment align=new Alignment(fastas);
            Infos.println(align.describeAlignment(false));
            fp.closePointer();
            //////////////////////////////////////////////
            //LOAD THE POSTERIOR PROBAS AND PAML TREE IDS
            Infos.println("Loading PAML tree ids and Posterior Probas...");
            double startParsingTime=System.currentTimeMillis();
            PAMLWrapper pw=new PAMLWrapper(align, s); //align extended or not by the relaxed bloc
            FileInputStream input = null;
            //input = new FileInputStream(new File(pp));
            input = new FileInputStream(new File(rst));
            //input = new FileInputStream(new File(ARPath+"rst"));
            PhyloTree tree= pw.parseTree(input);
            Infos.println("Parsing posterior probas..");
            input = new FileInputStream(new File(rst));
            PProbasSorted pprobas = pw.parseSortedProbas(input,sitePPThreshold,true,Integer.MAX_VALUE); //parse less nodes for debug
            input.close();
            double endParsingTime=System.currentTimeMillis();
            Infos.println("Word search took "+(endParsingTime-startParsingTime)+" ms");
    
    
            //positions for which word are checked
            SequenceKnife knife=new SequenceKnife(new String(align.getCharMatrix()[0]), k, k, s, SequenceKnife.SAMPLING_LINEAR);
            int[] refPositions=knife.getMerOrder();     
            
            
            //to raidly check that sorted probas are OK
            Infos.println("NodeId=0, 5 first PP:"+Arrays.deepToString(pprobas.getPPSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateSet(0, 0, 5)));
           
            
            //prepare simplified hash
            CustomHash hash=new CustomHash();

            Infos.println("Word generator threshold will be:"+thresholdAsLog);
            //Word Explorer
            int totalWordsInHash=0;
            double startHashBuildTime=System.currentTimeMillis();
            for (int nodeId:tree.getInternalNodesByDFS()) {
                Infos.println("NodeId: "+nodeId+" "+tree.getById(nodeId).toString() );
                
                
                double startMerScanTime=System.currentTimeMillis();
                WordExplorer wd =null;
                int totalWordsInNode=0;
                for (int pos:knife.getMerOrder()) {

                    if(pos+k-1>align.getLength()-1)
                        continue;
                    //DEBUG
//                    if(pos<30 || pos>36)
//                        continue;
                    //DEBUG
                    
                    //System.out.println("Current align pos: "+pos +" to "+(pos+(k-1)));
                    //double startScanTime=System.currentTimeMillis();
                    wd =new WordExplorer(   k,
                                            pos,
                                            nodeId,
                                            pprobas,
                                            thresholdAsLog
                                        );
                    
                    for (int j = 0; j < pprobas.getStateCount(); j++) {
                        wd.exploreWords(pos, j);
}
                    
                    //register the words in the hash
                    wd.words.stream().forEach((w)-> {hash.addTuple(w, w.getPpStarValue(), nodeId, w.getOriginalPosition());});
                    totalWordsInNode+=wd.words.size();
                    
                    //Infos.println("Words in this position:"+wd.words.size());
                    
                    //double endScanTime=System.currentTimeMillis();
                    //Infos.println("Word search took "+(endScanTime-startScanTime)+" ms");
                    
                    wd=null;
                }
                totalWordsInHash+=totalWordsInNode;
                
                
                //register all words in the hash
                Infos.println("Words in this node:"+totalWordsInNode);
                double endMerScanTime=System.currentTimeMillis();
                Infos.println("Word search took "+(endMerScanTime-startMerScanTime)+" ms");
                Environement.printMemoryUsageDescription();
                
            }

            
            Infos.println("Resorting hash components...");
            hash.sortTuples();
            
            double endHashBuildTime=System.currentTimeMillis();
            Infos.println("Overall, hash built took: "+(endHashBuildTime-startHashBuildTime)+" ms");
            
            Infos.println("Words in the hash:"+totalWordsInHash);
                
            Infos.println("FINISHED.");
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(WordExplorer.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(WordExplorer.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
    }
    
    
}
