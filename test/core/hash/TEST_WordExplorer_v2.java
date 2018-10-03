/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import alignement.Alignment;
import core.DNAStates;
import core.PProbasSorted;
import core.States;
import core.algos.SequenceKnife;
import core.algos.WordExplorer;
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
public class TEST_WordExplorer_v2 {

    
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
            Alignment align=new Alignment(s,fastas);
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
            CustomHash_v2 hash=new CustomHash_v2(k, s, CustomHash_v2.NODES_POSITION);

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
                                            thresholdAsLog,
                                            false,
                                            s
                                        );
                    
                    for (int j = 0; j < pprobas.getStateCount(); j++) {
                        wd.exploreWords(pos, j);
                    }
                    
                    //register the words in the hash
                    wd.getRetainedWords().stream().forEach((w)-> {hash.addTuple(w.getWord(), w.getPpStarValue(), nodeId, w.getOriginalPosition());});
                    totalWordsInNode+=wd.getRetainedWords().size();
                    
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
            hash.sortData();
            
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
