/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import alignement.Alignment;
import core.DNAStates;
import core.PProbasSorted;
import core.SimpleWord;
import core.States;
import core.Word;
import core.algos.SequenceKnife;
import core.algos.WordExplorer;
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
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import tree.Tree;

/**
 * In this second hash version, the Nodes composing the buckets contain a table
 * of reference positions, which point to a list of Pairs representing 
 * (nodeId, PP*)
 * @author ben
 */
public class SimpleHash_v2 implements Serializable{
    
    private static final long serialVersionUID = 7000L;
    
    CustomHashMap<Word,CustomNode> hash=null;

    public SimpleHash_v2() {
        hash=new CustomHashMap<>();
    }
    
    
    public SimpleHash_v2(int k,States s) {
        hash=new CustomHashMap<>(new Double(Math.pow(s.getNonAmbiguousStatesCount()+1, k)).intValue(),1.0f);
        System.out.println("INIT HASH("+(new Double(Math.pow(s.getNonAmbiguousStatesCount(), k)).intValue()+1)+",1.0)");
    }
    
    public void addTuple(Word w, float PPStar,int nodeId,int refPos) {
        if (!hash.containsKey(w)) {
            hash.put(w, new CustomNode());
        }
        hash.get(w).registerTuple(nodeId, refPos, PPStar);
    }

    public CustomHashMap<Word, CustomNode> getHash() {
        return hash;
    }

    public Pair getTopPair(Word w) {
        return hash.get(w).getBestPair();
    }    
    
    public void sortData() {
        double startTime=System.currentTimeMillis();
        hash.values().stream().forEach( l -> {l.sort();} );
        double endTime=System.currentTimeMillis();
        Infos.println("Tuples sorting took "+(endTime-startTime)+" ms");
    }
    
    
    /**
     * pairs whatever the associated position
     * @param w
     * @return 
     */
    public List<Pair> getPairs(Word w) {
        ArrayList<Pair> l =new ArrayList();
        for (int p: hash.get(w).getPositions())
            l.addAll(hash.get(w).getPairList(p));
        return l;
    }
    
    /**
     * pairs selected by associated position
     * @param w
     * @param position
     * @return 
     */
    public List<Pair> getTuples(Word w,int position) {
        if (hash.containsKey(w)) {
            return hash.get(w).getPairList(position);
        } else {
            return null;
        }
    }

    public Set<Word> getKeys() {
        return hash.keySet();
    }
    
    
    
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose", "1");
        
        try {
            
            String HOME = System.getenv("HOME");

            //DATASET BASIC RAPID TESTS:
            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD2";
            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/";
            String a=inputsPath+"mod_mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta";
            String t=inputsPath+"RAxML_bipartitionsBranchLabels.result_alpha_RNApol_REROOTED.tree";


            String rst=workDir+"/AR/rst";
            
            
            int k=8;
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
            
            //write log of word count per node/position
            Infos.println("Writing word counts in log...");
            File logWordCount=new File("log_word_count_k"+k+"_fact"+thresholdFactor);
            FileWriter fw=new FileWriter(logWordCount);
            fw.append("nodeId");
            for (int i=0;i<refPositions.length;i++) {
                fw.append("\t"+i);
            }
            fw.append("\n");
            
            
            
            ///////////////////////////////////////////////////::
            ///// TEST ON HASH V1
            boolean testBuckets=false;
            
            //prepare simplified hash
            SimpleHash hash=new SimpleHash();

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
                int totalTuplesPassingThreshold=0;
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
                    totalTuplesPassingThreshold+=wd.getRetainedWords().size();
                    
                    //wd.getRetainedWords().stream().forEach((w)->System.out.println(w));
                    //Infos.println("Words in this position:"+wd.getRetainedWords().size());
                    
                    //double endScanTime=System.currentTimeMillis();
                    //Infos.println("Word search took "+(endScanTime-startScanTime)+" ms");
                    
                    wd=null;
                }
                wordsPerNode[nodeCounter]=totalTuplesPassingThreshold;
                totalTuplesInHash+=totalTuplesPassingThreshold;
                
                
                //register all words in the hash
                //Infos.println("Tuples in this node:"+totalTuplesPassingThreshold);
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

            fw.close();
            
            
            
            
            
            
            ///////////////////////////////////////////////////::
            ///// TEST ON HASH V2
            testBuckets=false;
            
            //prepare simplified hash
            SimpleHash_v2 hash2=new SimpleHash_v2();

            Infos.println("Word generator threshold will be:"+wordPPStarThresholdAsLog);
            Infos.println("Building all words probas...");
            //Word Explorer
            totalTuplesInHash=0;
            nodeCounter=0;
            wordsPerNode=new double[tree.getInternalNodesByDFS().size()];
            startHashBuildTime=System.currentTimeMillis();
            for (int nodeId:tree.getInternalNodesByDFS()) {
                Infos.println("NodeId: "+nodeId+" "+tree.getById(nodeId).toString() );
                    //DEBUG
                    //if(nodeId!=709)
                    //    continue;
                    //DEBUG                
                
                double startMerScanTime=System.currentTimeMillis();
                WordExplorer wd =null;
                int totalTuplesPassingThreshold=0;
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
                    wd.getRetainedWords().stream().forEach((w)-> {hash2.addTuple(w, w.getPpStarValue(), nodeId, w.getOriginalPosition());});
                    totalTuplesPassingThreshold+=wd.getRetainedWords().size();
                    
                    //wd.getRetainedWords().stream().forEach((w)->System.out.println(w));
                    //Infos.println("Words in this position:"+wd.getRetainedWords().size());
                    
                    //double endScanTime=System.currentTimeMillis();
                    //Infos.println("Word search took "+(endScanTime-startScanTime)+" ms");
                    
                    wd=null;
                }
                wordsPerNode[nodeCounter]=totalTuplesPassingThreshold;
                totalTuplesInHash+=totalTuplesPassingThreshold;
                
                
                //register all words in the hash
                //Infos.println("Tuples in this node:"+totalTuplesPassingThreshold);
                double endMerScanTime=System.currentTimeMillis();
                //Infos.println("Word search took "+(endMerScanTime-startMerScanTime)+" ms");
                //Environement.printMemoryUsageDescription();
                nodeCounter++;
                
            }

            
            Infos.println("Resorting hash_v2 components...");
            hash2.sortData();
            
            endHashBuildTime=System.currentTimeMillis();
            Infos.println("Overall, hash_v2 built took: "+(endHashBuildTime-startHashBuildTime)+" ms");
            Infos.println("Words in the hash_v2: "+hash.getKeys().size());
            Infos.println("Tuples in the hash_v2:"+totalTuplesInHash);

            
            //HERE Search how many Nodes per hash bucket.
            if (testBuckets) {
                CustomHashMap.Node<Word, CustomNode>[] accessToHash = hash2.getHash().getAccessToHash();
                for (int i = 0; i < accessToHash.length; i++) {
                    Object node = accessToHash[i];
                    if (node instanceof core.hash.CustomHashMap.TreeNode) {
                        CustomHashMap.TreeNode n=(CustomHashMap.TreeNode)node;
                        System.out.println("i"+i+"=TreeNode: "+n.getValue()+" left:"+n.getLeft()+" right:"+n.getRight());
                        System.out.println("Size:"+hash2.bucketDFS(n));
                    } else if (node instanceof core.hash.CustomHashMap.Node) {
                        if (node!=null)
                            System.out.println("i"+i+"=Node: "+node.getClass().getName());
                        else
                            System.out.println(node);
                    } else {
                        System.out.println(node);
                    }
                }
            }
            fw.close();
            
            
            ///////////////////////////////////////////////////::
            ///// CHECK IF THEY CONTAIN SAME VALUES
            System.out.println("##########################################");
            
            byte[] word={1, 0, 0, 0, 3, 1, 2, 0}; //matches positions 1717,983,2149
            Word w=new SimpleWord(word);
            
            //all words in hash
            //hash.getKeys().stream().forEach(key->{System.out.println(key);});
            
            //number of positions per word
            //hash2.getHash().entrySet().stream().forEach(e->{System.out.println(e.getKey()+" --> "+e.getValue().getPositions().length);});
            
            System.out.println("-----");
            System.out.println(hash.getTuples(w).toString().replaceAll(",", "\n"));
            System.out.println("-----");
            System.out.println(hash2.getTuples(w, 1717).toString().replaceAll(",", "\n"));
            System.out.println(hash2.getTuples(w, 983).toString().replaceAll(",", "\n"));
            System.out.println(hash2.getTuples(w, 2149).toString().replaceAll(",", "\n"));
            
            
            
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
 
    
    

    

    
    
    
    
}
