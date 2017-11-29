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
    int[] gapIntervals=null;

    /**
     *
     * @param k
     * @param refPosition
     * @param nodeId
     * @param ppSet
     * @param wordThresholdAsLog
     */
    public WordExplorer_v2(int k, int refPosition, int nodeId, PProbasSorted ppSet, float wordThresholdAsLog, boolean wordCompression, States s) {
        this.refPosition=refPosition;
        this.wordThresholdAsLog=wordThresholdAsLog;
        this.word=new byte[k];
        this.k=k;
        this.nodeId=nodeId;
        this.ppSet=ppSet;
        this.wordCompression=wordCompression;
        this.s=s;
        this.current_k=0;
        
        //example of regions where 1 gapIntervals interval of 1 site
        //and 1 gapIntervals interval of 2 sites
        int[] intervalStarts={2,4};
        int[] gapIntervals=new int[ppSet.getSiteCount()];
        gapIntervals[2]=1;
        gapIntervals[4]=2;
        //{0,0,1,0,2,0,0,0}; //this table should be built in constructor
        
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
        

        System.out.println("IN: "+i+" "+j);
        System.out.println("current_k:"+current_k);
        //if after the alignment limit (can happen because of jumps
        //over gapIntervals colulmns), then return, this word is not possible.
        if (i>ppSet.getSiteCount()-1) {
            return;
        }
        //if we start the exploration (current_k==0)at a i which is in a 
        //gap interval, we avoid this exporation.
        if ( (current_k==0) && gapIntervals[i]>0 ) {
            System.out.println("statring in interval, return");
            return;
        }
        
        
        word[current_k]=ppSet.getState(nodeId, i, j);
        System.out.println("sumCurrentWord="+currentLogSum+"+"+ppSet.getPP(nodeId, i, j));
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
                System.out.println("REGISTER: "+Arrays.toString(word)+" log10(PP*)="+currentLogSum);
            }
            //decrease before return
            System.out.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
            currentLogSum-=ppSet.getPP(nodeId, i, j);
            //force return to parent
            return;
        } else {
            //go down recursively
            for (int j2 = 0; j2 < ppSet.getStateCount(); j2++) {
                if (boundReached) {break;}
                
                //remember with move to next kmer position
                current_k++;
                
                //continue exploration by jumping over the next column
                //which is flagged as a gapIntervals interval
                if (gapIntervals[i+1]!=0) {
                    exploreWords((i+1)+gapIntervals[i+1], j2);
                }
                //redo this exploration, but now considering this gapIntervals zone
                exploreWords(i+1, j2);
                
                //remember with move back to previous kmer position
                current_k--;
                
                
                System.out.println("OUT: "+(i+1)+" "+j2);
                System.out.println("current_k:"+current_k);
            }
        }
        //decrease before natural return
        System.out.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
        currentLogSum-=ppSet.getPP(nodeId, i, j);
        
    }
    
    
    
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose", "1");
        
        try {
            
            
            String a="/media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/mod_matK.fasta.linsi.aln.reduced.fasta";
            String ARTree="/media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/mod_matK.fasta.linsi.aln.reduced_phyml_ancestral_tree.txt";
            String ARStats="/media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/mod_matK.fasta.linsi.aln.reduced_phyml_ancestral_seq.txt";
            
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
//                    if(pos<30 || pos>36)
//                        continue;
                    //DEBUG
                    
                    //System.out.println("Current align pos: "+pos +" to "+(pos+(k-1)));
                    //double startScanTime=System.currentTimeMillis();
                    wd =new WordExplorer_v2(   k,
                                            pos,
                                            nodeId,
                                            pprobas,
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
