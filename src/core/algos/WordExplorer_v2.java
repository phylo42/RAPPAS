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
import core.hash.CustomHash;
import core.States;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.PAMLWrapper;
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
    ArrayList<Integer>[] gapIntervals =null;

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
     */
    public WordExplorer_v2(int k, int refPosition, int nodeId, PProbasSorted ppSet, Alignment align,float wordThresholdAsLog, boolean wordCompression, States s) {
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
        

        Infos.println("IN: "+i+" "+j);
        Infos.println("current_k:"+current_k);
        //if after the alignment limit (can happen because of jumps
        //over gapIntervals colulmns), then return, this word is not possible.
        if (i>ppSet.getSiteCount()-1) {
            return;
        }
        //if we start the exploration (current_k==0)at a i which is in a 
        //gap interval, we avoid this exporation.
        if ( (current_k==0) && gapIntervals[i]!=null ) {
            Infos.println("starting in interval, return");
            return;
        }
        
        
        word[current_k]=ppSet.getState(nodeId, i, j);
        Infos.println("sumCurrentWord="+currentLogSum+"+"+ppSet.getPP(nodeId, i, j));
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
                Infos.println("REGISTER: "+Arrays.toString(word)+" log10(PP*)="+currentLogSum);
            }
            //decrease before return
            Infos.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
            currentLogSum-=ppSet.getPP(nodeId, i, j);
            //force return to parent
            return;
        } else {
            //go down recursively
            for (int j2 = 0; j2 < ppSet.getStateCount(); j2++) {
                if (boundReached) {break;}
                
                //remember with move to next kmer position
                current_k++;
                
                //do this exploration, do not consider gap interval
                Infos.println("EXPORATION NORMAL");
                System.setProperty("debug.verbose", "1");
                exploreWords(i+1, j2);
                System.setProperty("debug.verbose", "1");
                //redo exploration, by jumping over the gap intervals 
                if (gapIntervals[i+1]!=null) {
                    for (int i_interval = 0; i_interval < gapIntervals[i+1].size(); i_interval++) {
                        Infos.println("EXPLORATION WITH GAP JUMP: gap_i+1="+(i+1)+" length="+gapIntervals[i+1].get(i_interval));
                        Infos.println("NEXT i jump to --> "+(i+gapIntervals[i+1].get(i_interval)));
                        exploreWords(i+gapIntervals[i+1].get(i_interval), j2);
                    }
                }
                System.setProperty("debug.verbose", "1");
                
                //remember with move back to previous kmer position
                current_k--;
                
                
                Infos.println("OUT: "+(i+1)+" "+j2);
                Infos.println("HERE: "+i+" "+j2);
                Infos.println("current_k:"+current_k);
            }
        }
        //decrease before natural return
        Infos.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
        currentLogSum-=ppSet.getPP(nodeId, i, j);
        
    }
    
    
    
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose", "1");
        
        try {
            
            
            String a="/home/ben/Dropbox/viromeplacer/test_datasets/mod_matK.fasta.linsi.aln.reduced.fasta";
            String ARTree="/home/ben/Dropbox/viromeplacer/test_datasets/mod_matK.fasta.linsi.aln.reduced_phyml_ancestral_tree.txt";
            String ARStats="/home/ben/Dropbox/viromeplacer/test_datasets/mod_matK.fasta.linsi.aln.reduced_phyml_ancestral_seq.txt";
            
            int k=6;
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
            PHYMLWrapper pw=new PHYMLWrapper(align, s);
            //PAMLWrapper pw=new PAMLWrapper(align, s); //align extended or not by the relaxed bloc
            FileInputStream input = null;
            //input = new FileInputStream(new File(pp));
            input = new FileInputStream(new File(ARTree));
            //input = new FileInputStream(new File(ARPath+"tree"));
            PhyloTree tree= pw.parseTree(input, false);
            Infos.println("Parsing posterior probas..");
            input = new FileInputStream(new File(ARStats));
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
            //to rapidly check gap intervals
            ArrayList<Integer>[] gapIntervals1 = align.getGapIntervals();
            for (int i = 0; i < 107; i++) {
                ArrayList<Integer> arrayList = gapIntervals1[i];
                Infos.println("Gap interval i="+i+" (pos="+i+1+") : "+gapIntervals1[i]);
            }
           
            
            //prepare simplified hash
            CustomHash hash=new CustomHash();

            Infos.println("Word generator threshold will be:"+thresholdAsLog);
            //Word Explorer
            int totalWordsInHash=0;
            double startHashBuildTime=System.currentTimeMillis();
            for (int nodeId:tree.getInternalNodesByDFS()) {
                Infos.println("NodeId: "+nodeId+" "+tree.getById(nodeId).toString() );
                
                
                double startMerScanTime=System.currentTimeMillis();
                WordExplorer_v2 wd =null;
                int totalWordsInNode=0;
                for (int pos:knife.getMerOrder()) {

                    if(pos+k-1>align.getLength()-1)
                        continue;
                    //DEBUG
                    if(pos<79 || pos>107)
                        continue;
                    //DEBUG
                    System.setProperty("debug.verbose", "1");
                    
                    //double startScanTime=System.currentTimeMillis();
                    wd =new WordExplorer_v2(   k,
                                            pos,
                                            nodeId,
                                            pprobas,
                                            align,
                                            thresholdAsLog,
                                            false,
                                            s
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
            Logger.getLogger(WordExplorer_v2.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(WordExplorer_v2.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
    }
    
    
}
