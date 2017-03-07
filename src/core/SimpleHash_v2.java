/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import core.older.PProbas;
import alignement.Alignment;
import core.algos.SequenceKnife;
import core.algos.WordExplorer;
import core.algos.WordGenerator;
import core.hashmap.CustomHashMap;
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
import tree.Tree;

/**
 *
 * @author ben
 */
public class SimpleHash_v2 implements Serializable{
    
    private static final long serialVersionUID = 7000L;
    
    CustomHashMap<Word,LinkedList<Tuple>> hash=null;

    public SimpleHash_v2() {
        hash=new CustomHashMap<>();
    }
    
    
    public SimpleHash_v2(int k,States s) {
        hash=new CustomHashMap<>(new Double(Math.pow(s.getNonAmbiguousStatesCount(), k)*0.75).intValue());
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
    
    public List<Tuple> getTopTuples(Word w,float PPStarTresholdAsLog10) {
        return hash.get(w).stream().filter(t -> t.PPStar>=PPStarTresholdAsLog10).collect(Collectors.toList());
    }

    public CustomHashMap<Word, LinkedList<Tuple>> getHash() {
        return hash;
    }
    
    
    
    /**
     * 
     * @param w
     * @param PPStarTresholdAsLog10
     * @param nodeIdsTested  a boolean table, with true in the positions corresponding to the nodeIds which are used for preplacement
     * @return 
     */
    public List<Tuple> getTopTuplesUnderNodeShift(Word w,float PPStarTresholdAsLog10, boolean[] nodeIdsTested) {
        return hash.get(w).stream().filter(t -> ( t.PPStar>=PPStarTresholdAsLog10 && nodeIdsTested[t.nodeId] )).collect(Collectors.toList());
    }    
    public Tuple getTopTuple(Word w) {
        LinkedList<Tuple> l;
        return (l=hash.get(w))==null ? null : l.getFirst();
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
            
            

            String a="/media/ben/STOCK/DATA/viromeplacer/WD/relaxed_trees/relaxed_align_BrB_minbl0.001_1peredge.fasta";
            String t="/media/ben/STOCK/DATA/viromeplacer/WD/relaxed_trees/relaxed_tree_BrB_minbl0.001_1peredge_withBL_withoutInternalLabels.tree";
            String rst="/media/ben/STOCK/DATA/viromeplacer/WD/AR/rst";
            
            
            int k=12;
            float sitePPThreshold=1e-45f;
            float thresholdFactor=2.0f;
            float wordPPStarThreshold=(float)Math.pow(thresholdFactor*0.25,k);
            float wordPPStarThresholdAsLog=(float)Math.log10(wordPPStarThreshold);
            
            Infos.println("k="+k);
            Infos.println("factor="+thresholdFactor);
            Infos.println("wordPPStatThreshold="+wordPPStarThreshold);
            Infos.println("wordPPStatThresholdAsLog="+wordPPStarThresholdAsLog);
            
            
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
            Tree tree= pw.parseTree(input);
            Infos.println("Parsing posterior probas..");
            input = new FileInputStream(new File(rst));
            PProbasSorted pprobas = pw.parseSortedProbas(input, sitePPThreshold, true,Integer.MAX_VALUE);  //DEBUG: LIMITED NODE NUMBER
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
            
          //prepare simplified hash
            SimpleHash_v2 hash=new SimpleHash_v2();

            Infos.println("Word generator threshold will be:"+wordPPStarThresholdAsLog);
            Infos.println("Building all words probas...");
            //Word Explorer
            int totalTuplesInHash=0;
            int nodeCounter=0;
            double[] wordsPerNode=new double[tree.getInternalNodesByDFS().size()];
            double startHashBuildTime=System.currentTimeMillis();
            for (int nodeId:tree.getInternalNodesByDFS()) {
                Infos.println("NodeId: "+nodeId+" "+tree.getById(nodeId).toString() );
                    //DEBUG
                    //if(nodeId!=709)
                    //    continue;
                    //DEBUG                
                
                double startMerScanTime=System.currentTimeMillis();
                WordExplorer wd =null;
                int totaTuplesInNode=0;
                for (int pos:knife.getMerOrder()) {

                    if(pos+k-1>align.getLength()-1)
                        continue;
                    //DEBUG
                    //if(pos<1198 || pos>1202)
                    //    continue;
                    //DEBUG
                    
                    //System.out.println("Current align pos: "+pos +" to "+(pos+(k-1)));
                    double startScanTime=System.currentTimeMillis();
                    wd =new WordExplorer(   k,
                                            pos,
                                            nodeId,
                                            pprobas,
                                            wordPPStarThresholdAsLog
                                        );
                    
                    for (int j = 0; j < pprobas.getStateCount(); j++) {
                        wd.exploreWords(pos, j);
                    }
                    
                    //register the words in the hash
                    wd.getRetainedWords().stream().forEach((w)-> {hash.addTuple(w, w.getPpStarValue(), nodeId, w.getOriginalPosition());});
                    totaTuplesInNode+=wd.getRetainedWords().size();
                    
                    //wd.getRetainedWords().stream().forEach((w)->System.out.println(w));
                    //Infos.println("Words in this position:"+wd.getRetainedWords().size());
                    
                    //double endScanTime=System.currentTimeMillis();
                    //Infos.println("Word search took "+(endScanTime-startScanTime)+" ms");
                    
                    wd=null;
                }
                wordsPerNode[nodeCounter]=totaTuplesInNode;
                totalTuplesInHash+=totaTuplesInNode;
                
                
                //register all words in the hash
                //Infos.println("Tuples in this node:"+totaTuplesInNode);
                double endMerScanTime=System.currentTimeMillis();
                //Infos.println("Word search took "+(endMerScanTime-startMerScanTime)+" ms");
                //Environement.printMemoryUsageDescription();
                nodeCounter++;
                
            }

            
            Infos.println("Resorting hash components...");
            hash.sortTuples();
            
            double endHashBuildTime=System.currentTimeMillis();
            Infos.println("Overall, hash built took: "+(endHashBuildTime-startHashBuildTime)+" ms");
            Infos.println("Words in the hash: "+hash.getKeys().size());
            Infos.println("Tuples in the hash:"+totalTuplesInHash);

            
            //HERE Search how many Nodes per hash bucket.
            CustomHashMap.Node<Word, LinkedList<Tuple>>[] accessToHash = hash.getHash().getAccessToHash();
            for (int i = 0; i < accessToHash.length; i++) {
                Object node = accessToHash[i];
                if (node instanceof core.hashmap.CustomHashMap.TreeNode) {
                    CustomHashMap.TreeNode n=(CustomHashMap.TreeNode)node;
                    System.out.println("i"+i+"=tree: "+n.getValue()+" left:"+n.getLeft()+" right:"+n.getRight());
                    System.out.println("Size:"+hash.bucketDFS(n));
                } else {
                    if (node!=null)
                        System.out.println("i"+i+"=other: "+node.getClass().getName());
                    else
                        System.out.println(node);
                }
            }

            
            
            fw.close();
            
            Infos.println("FINISHED.");
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(SimpleHash_v2.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(SimpleHash_v2.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        
        
        
    }
    
    
    
    
    /**
     * Basic DFS to explore the RED/BLACK tree associated to buckets
     * @param root 
     */
    int DFSCount=0;
    private int bucketDFS(CustomHashMap.TreeNode root) {
        if (root.getLeft()!=null)
            bucketDFS(root.getLeft());
        if (root.getRight()!=null)
            bucketDFS(root.getRight());
        return DFSCount++;
    }
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    
    
    public class Bucket implements Serializable {
    
        private static final long serialVersionUID = 7020L;
    
    
    
    }
    
    
    
    
    
    
    
    
    
    
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
