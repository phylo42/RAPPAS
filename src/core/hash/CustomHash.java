/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import core.older.PProbas;
import alignement.Alignment;
import core.DNAStates;
import core.ProbabilisticWord;
import core.States;
import core.Word;
import core.algos.SequenceKnife;
import core.algos.WordGenerator;
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
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import tree.NewickWriter;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class CustomHash implements Serializable{
    
    private static final long serialVersionUID = 7000L;
    
    HashMap<Word,LinkedList<Tuple>> hash=null;

    
    public CustomHash() {
        hash=new HashMap<>();
    }
    
    
    public void addTuple(Word w, float PPStar,int nodeId,int refPos) {
        if (!hash.containsKey(w)) {
            hash.put(w, new LinkedList<>());
        }
        hash.get(w).add(new Tuple(PPStar,nodeId,refPos));
    }
    
    public LinkedList<Tuple> getAllTuples(Word w) {
        return hash.get(w);
    }
    
    /**
     * return the top tuples of PP*> given threshold
     * @param w
     * @param PPStarTresholdAsLog10
     * @return 
     */
    public List<Tuple> getTopTuples(Word w,float PPStarTresholdAsLog10) {
        return hash.get(w).stream().filter(t -> t.PPStar>=PPStarTresholdAsLog10).collect(Collectors.toList());
    }
    
    /**
     * return the top tuples of PP*> given threshold and for the given nodes
     * @param w
     * @param PPStarTresholdAsLog10
     * @param nodeIdsTested  a boolean table, with true in the positions corresponding to the nodeIds which are used for preplacement
     * @return 
     */
    public List<Tuple> getTopTuplesUnderNodeShift(Word w,float PPStarTresholdAsLog10, boolean[] nodeIdsTested) {
        return hash.get(w).stream().filter(t -> ( t.PPStar>=PPStarTresholdAsLog10 && nodeIdsTested[t.nodeId] )).collect(Collectors.toList());
    }
    /**
     * return the tuple of highest PP*, whatever the node and position
     * @param w
     * @return 
     */    
    public Tuple getTopTuple(Word w) {
        LinkedList<Tuple> l;
        return (l=hash.get(w))==null ? null : l.getFirst();
    }    
    
    /**
     * return a tuple selected by reference position and nodeId or null is not
     * in hash
     * @param w
     * @param nodeId
     * @param refPosition
     * @return 
     */
    public Tuple getTuplePerNodeAndRefPosition(Word w,int nodeId, int refPosition) {
        return hash.get(w).stream().filter(t -> ( t.refPos==refPosition && t.nodeId==nodeId )).findFirst().orElse(null);
    }
    
    /**
     * return a list of tuples selected by reference position
     * @param refPosition
     * @return 
     */
    public List<Tuple> getTuplePerRefPosition(Word w,int refPosition) {
        return hash.get(w).stream().filter(t -> ( t.refPos==refPosition)).collect(Collectors.toList());
    }
    
    public void sortTuples() {
        double startTime=System.currentTimeMillis();
        hash.values().stream().forEach( l -> {Collections.sort(l);} );
        double endTime=System.currentTimeMillis();
        Infos.println("Tuples sorting took "+(endTime-startTime)+" ms");
    }
    
    public Set<Word> getRegisteredWords() {
        return hash.keySet();
    }
          
    public List<Tuple> getTuples(Word w) {
        return hash.get(w);
    }

    public Set<Word> getKeys() {
        return hash.keySet();
    }
    
    
    
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose", "1");
        
        try {
            
            
            String a="/media/ben/STOCK/SOURCES/NetBeansProjects/ViromePlacer/WD/relaxed_trees/relaxed_align_BrB_minbl0.001_1peredge.fasta";
            String rst="/media/ben/STOCK/SOURCES/NetBeansProjects/ViromePlacer/WD/AR/finished/rst";
            
            
            int k=7;
            float sitePPThreshold=1e-45f;
            float wordPPStarThreshold=(float)Math.pow(0.25,k);
            int thresholdFactor=10;
            
            Infos.println("k="+k);
            Infos.println("factor="+thresholdFactor);
            Infos.println("wordPPStatThreshold="+wordPPStarThreshold);
            
            
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
            PAMLWrapper pw=new PAMLWrapper(align, s); //align extended or not by the relaxed bloc
            
            FileInputStream input = null;
            //input = new FileInputStream(new File(pp));
            input = new FileInputStream(new File(rst));
            //input = new FileInputStream(new File(ARPath+"rst"));
            PhyloTree tree= pw.parseTree(input);
            Infos.println("Parsing posterior probas..");
            input = new FileInputStream(new File(rst));
            PProbas pprobas = pw.parseProbas(input,sitePPThreshold,false);
            input.close();
            
            
            //positions for which word are checked
            SequenceKnife knife=new SequenceKnife(new String(align.getCharMatrix()[0]), k, k, s, SequenceKnife.SAMPLING_LINEAR);
            int[] refPositions=knife.getMerOrder();            
            //Word generator
            WordGenerator wg=new WordGenerator();
            
            //write log of word count per node/position
            Infos.println("Writing word counts in log...");
            File logWordCount=new File("log_word_count_k"+k+"_fact"+thresholdFactor);
            FileWriter fw=new FileWriter(logWordCount);
            fw.append("nodeId");
            for (int i=0;i<refPositions.length;i++) {
                fw.append("\t"+i);
            }
            fw.append("\n");
            
            //HASH
            CustomHash sh=new CustomHash();
            
            for (int nodeId=0;nodeId<pprobas.accessTable().length;nodeId++) {
                Infos.println("##### NodeId="+nodeId);

                fw.append(String.valueOf(nodeId));
                
                for (int i=0;i<refPositions.length;i++) {
                    int refPosition=refPositions[i];
                    
                    if ((refPosition+k) > (align.getLength()-1))
                        break;
                    
                    ArrayList<ProbabilisticWord> words = wg.generateProbableWords3(nodeId, pprobas.getPPSet(nodeId, refPosition, refPosition+k),1e-323,wordPPStarThreshold);
                    
//                    if (words.size()>5000) {
//                        Infos.println("#words > treshold at "+i+": "+words.size());
//                        words.stream().forEach(w->{System.out.println(Arrays.toString(w.word)+":"+w.getPpStarValue()+" ("+w.getOriginalPosition()+")");});
//                        System.exit(1);
//                    }
                    //Infos.println("#words at "+i+": "+words.size());
                    int wordPassingThreshold=0;
                    for (Iterator<ProbabilisticWord> iterator = words.iterator(); iterator.hasNext();) {
                        ProbabilisticWord w = iterator.next();
                        if (w.getPpStarValue()>wordPPStarThreshold*thresholdFactor) {
                            sh.addTuple(w, w.getPpStarValue(), nodeId, refPosition );
                            wordPassingThreshold++;
                        }
                    }
                    fw.append("\t"+wordPassingThreshold);
  
                }
                
                Environement.printMemoryUsageDescription();
                fw.append("\n");
            }
            //sorting all nodes
            Infos.println("Sorting PP* per node...");
            for (Iterator<LinkedList<Tuple>> it=sh.hash.values().iterator();it.hasNext();) {
                LinkedList l=it.next();
                Collections.sort(l);
            }
            
            
            

            
            
            fw.close();
            
            Infos.println("FINISHED.");
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(CustomHash.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(CustomHash.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        
        
        
    }
    
    
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * simple object defining the caracteristics of particular word
     */
    public class Tuple implements Comparable<Tuple>,Serializable {
        
        private static final long serialVersionUID = 7010L;

        
        protected float PPStar=-1.0f;
        protected int nodeId=-1;
        protected int refPos=-1; 
        
        public Tuple(float PPStar,int nodeId, int refPos) {
            this.PPStar=PPStar;
            this.nodeId=nodeId;
            this.refPos=refPos;
        }

        public int getNodeId() {
            return nodeId;
        }

        public float getPPStar() {
            return PPStar;
        }

        public int getRefPos() {
            return refPos;
        }
        
        /**
        * the comparator is inversed to return highest values first when sorting
        * @param o
        * @return 
        */
        @Override
       public int compareTo(Tuple o) {
           if (this.PPStar-o.PPStar<0.0) {
               return 1;
           } else if (this.PPStar-o.PPStar>0.0){
               return -1;
           } else {
               return 0;
           }
       }

        @Override
        public String toString() {
            return "Tuple: nodeId="+nodeId+" refPos="+refPos+" PPStar="+PPStar;
        }
        public String toStringCSV() {
            return nodeId+","+refPos+","+PPStar;
        }
       
    }
    
    
    
    
}
