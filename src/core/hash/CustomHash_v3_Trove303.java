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
import etc.Infos;
import gnu.trove.map.hash.TCustomHashMap;
import gnu.trove.strategy.HashingStrategy;
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
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import tree.PhyloTree;

/**
 * In this second hash version, the Nodes composing the buckets contain a table
 * of positions pointing to a list of Pairs representing a
 * (nodeId, PP*)
 * @author ben
 */
public class CustomHash_v3_Trove303 implements Serializable{
    
    private static final long serialVersionUID = 7000L;
    
    public static final int NODES_POSITION=1;
    public static final int NODES_UNION=2;
    
    int nodeType=NODES_UNION;
    
    TCustomHashMap<byte[],HashPointer> hash;

    
    ArrayList<Pair> pairsBuffer =new ArrayList(1024);
    HashPointer tempPointer=null;
    
    int maxPointerSize=-1;
    long addTupleCalls=0l;
    long totalContainsTime=0l;
    long totalRegisterTime=0l;
    /**
     *
     * @param k
     * @param s
     * @param nodeType one of NODES_UNION or POSITION_UNION
     */
    public CustomHash_v3_Trove303(int k, States s, int nodeType) {        
        this.maxPointerSize=new Double(Math.pow(s.getNonAmbiguousStatesCount(), k)).intValue();
        this.nodeType=nodeType;
        //internal tests showed that with k<14 we generally get at least 75% of the possible k-mers
        this.hash=new TCustomHashMap<>(
                new HashingStrategy<byte[]>() {
                    @Override
                    public int computeHashCode(byte[] object) {
                        return Arrays.hashCode(object);
                    }

                    @Override
                    public boolean equals(byte[] o1, byte[] o2) {
                        return Arrays.equals(o1, o2);
                    }
                },
                new Double(maxPointerSize*0.75).intValue(),  //intial capacity
                0.8f //inital load factor
        );
        //init first temp pointer
        switch (nodeType) {
            case NODES_UNION:
                tempPointer=new UnionPointer();
                break;
            case NODES_POSITION:
                tempPointer=new PositionPointer();
                break;
        }
        
    }
    
    /**
     * only entry point to fill the hash (used at DB generation)
     * @param word
     * @param PPStar
     * @param nodeId
     * @param refPos 
     */
    public void addTuple(byte[] word, float PPStar,int nodeId,int refPos) {

        long containsStart=System.nanoTime();
        
        if (!hash.containsKey(word)) {
            long containsEnd=System.nanoTime();
            totalContainsTime+=(containsEnd-containsStart);
            
            HashPointer hp=null;
            //we build a new pointer 
            switch (nodeType) {
                case NODES_UNION:
                    hp=new UnionPointer();
                    break;
                case NODES_POSITION:
                    hp=new PositionPointer();
                    break;
                default:
                    hp=new UnionPointer();
                    break;
            }     
            long registerStart=System.nanoTime();
            hp.registerTuple(nodeId, refPos, PPStar);
            hash.put(word, hp);
            long registerEnd=System.nanoTime();
            totalRegisterTime+=(registerEnd-registerStart);            
            
        } else {
            long containsEnd=System.nanoTime();
            totalContainsTime+=(containsEnd-containsStart);
            
            //update list of pairs
            long registerStart=System.nanoTime();
            hash.get(word).registerTuple(nodeId, refPos, PPStar);
            long registerEnd=System.nanoTime();
            totalRegisterTime+=(registerEnd-registerStart);
        }
        

        
        addTupleCalls++;
    }

    /**
     * for debug purposes only (access to CustomHashMap hashtable)
     * @return 
     */
    public TCustomHashMap<byte[], HashPointer> getHash() {
        return hash;
    }
    
    /**
     * to know which type of nodes is backing this hash, 
     * one of CustomHash_v2.NODES_UNION or CustomHash_v2.NODES_POSITION
     * @return 
     */
    public int getHashType() {
        return this.nodeType;
    }

    /**
     * best (nodeId,PP*) associated to this word
     * @param w
     * @return 
     */
    public Pair getTopPair(byte[] w) {
        HashPointer cn=null;
        if ((cn=hash.get(w))!=null) {
            return cn.getBestPair();
        } else {
            return null;
        }
    }    
    
    /**
     * retrieves only (nodeId,PP*) pairs stored under the position
     * associated to the best PP*
     * @param w
     * @return null if word not in present in hash
     */
    public List<Pair> getPairsOfTopPosition(byte[] w) {
        HashPointer cn=null;
        if ((cn=hash.get(w))!=null) {
            return cn.getPairList(cn.getBestPosition());
        } else {
            return null;
        }
    }  
    
    /**
     * reference alignment positions associated to a word
     * @param w
     * @return 
     */
    public int[] getPositions(byte[] w) {
        HashPointer cn=null;
        if ((cn=hash.get(w))!=null) {
            return cn.getPositions();
        } else {
            return null;
        }
    }
    
    /**
     * reference alignment positions associated to a word
     * @param w
     * @return -1 if word not in hash
     */
    public int getTopPosition(byte[] w) {
        HashPointer cn=null;
        if ((cn=hash.get(w))!=null) {
            return cn.getBestPosition();
        } else {
            return -1;
        }
    }
    
    
    
    public void sortData() {
        double startTime=System.currentTimeMillis();
        AtomicInteger pairCount=new AtomicInteger(0);
        hash.values().stream().forEach( l -> {l.sort();pairCount.addAndGet(l.getPairCountInTopPosition());} );
        double endTime=System.currentTimeMillis();
        Infos.println("# Pairs sorted in top positions: "+pairCount.get());
        Infos.println("Pair sorting took "+(endTime-startTime)+" ms");
    }
    
    
    /**
     * pairs whatever the associated position
     * @param w
     * @return 
     */
    public List<Pair> getPairs(byte[] w) {
        pairsBuffer.clear();
        for (int p: hash.get(w).getPositions()) {
            pairsBuffer.addAll(hash.get(w).getPairList(p));
        }
        return pairsBuffer;
    }
    
    /**
     * pairs selected by associated position
     * @param w
     * @param position
     * @return 
     */
    public List<Pair> getPairs(byte[] w,int position) {
        if (hash.containsKey(w)) {
            return hash.get(w).getPairList(position);
        } else {
            return null;
        }
    }

    public Set<byte[]> keySet() {
        return hash.keySet();
    }
    
    
    
    
    /**
     * empty all the positions which where not associated to the best PP*
     */
    public void reduceToMediumHash() {
        
        hash.keySet().stream().forEach((next) -> {
            ((PositionPointer)hash.get(next)).clearPairsOfWorsePositions();
        });
        
    }

    /**
     * in top position, empties all the positions which where not associated to the best PP*,
     * and retains only X pairs at the best position
     * @param X
     */
    @Deprecated
    public void reducetoSmallHash(int X) {
        List<byte[]> collect = hash.keySet()  .stream() 
                                            //.peek((w)->System.out.println("REDUCING:"+w))
                                            .filter((w) -> ((PositionPointer)hash.get(w)).limitToXPairsPerPosition(X))
                                            //.peek((w)->System.out.println("TRASHED!:"+w))
                                            .collect(Collectors.toList());
        collect.stream().forEach((w)-> {hash.remove(w);});
        collect=null;
    }
    
    /**
     * discards all words where top positions is associated to more than X nodes
     * @param X
     */
    public void reducetoSmallHash_v2(int X) {
        assert X>0;
        List<byte[]> collect = hash.keySet()  .stream() 
                                            //.peek((w)->System.out.println("REDUCING:"+w))
                                            .filter((w) -> hash.get(w).getPairCountInTopPosition()>X)
                                            //.peek((w)->System.out.println("TRASHED!:"+w))
                                            .collect(Collectors.toList());
        collect.stream().forEach((w)-> {hash.remove(w);});
        collect=null;
    }    
    
    
    public long getTotalContainsTime() {
        return this.totalContainsTime;
    }

    public long getTotalRegisterTime() {
        return this.totalRegisterTime;
    }
    
    public long getTotalAddTupleCalls() {
        return this.addTupleCalls;
    }
    
    
    ////////////////////////////////////////////////////////////////////////////
 
    
    
    //to test the class
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose", "1");
        
        try {
            
            String HOME = System.getenv("HOME");

            //DATASET BASIC RAPID TESTS:
//            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD2";
//            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/";
//            String a=inputsPath+"mod_mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta";
//            String t=inputsPath+"RAxML_bipartitionsBranchLabels.result_alpha_RNApol_REROOTED.tree";
            
            //DATASET LARGE TESTS:
            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD";
            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL";
            String a=inputsPath+"bv_refs_aln_stripped_99.5.fasta";
            String t=inputsPath+"RAxML_result.bv_refs_aln";
            
            
            
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
            PhyloTree tree= pw.parseTree(input);
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
            boolean doHashV1=true;
            final CustomHash hash=new CustomHash();
            
            if (doHashV1) {
                
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
                                                wordPPStarThresholdAsLog,
                                                false,
                                                s
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
            
            }
            
            
            
            
            ///////////////////////////////////////////////////::
            ///// TEST ON HASH V2
            testBuckets=false;
            
            //prepare simplified hash
            CustomHash_v3_Trove303 hash2=new CustomHash_v3_Trove303(k,s,NODES_POSITION);

            Infos.println("Word generator threshold will be:"+wordPPStarThresholdAsLog);
            Infos.println("Building all words probas...");
            //Word Explorer
            int totalTuplesInHash=0;
            int nodeCounter=0;
            double[] wordsPerNode=new double[tree.getInternalNodesByDFS().size()];
            long startHashBuildTime=System.currentTimeMillis();
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
                                            wordPPStarThresholdAsLog,
                                            false,
                                            s
                                        );
                    
                    for (int j = 0; j < pprobas.getStateCount(); j++) {
                        wd.exploreWords(pos, j);
                    }
                    
                    //register the words in the hash
                    wd.getRetainedWords().stream().forEach((w)-> {hash2.addTuple(w.getWord(), w.getPpStarValue(), nodeId, w.getOriginalPosition());});
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
            
            long endHashBuildTime=System.currentTimeMillis();
            Infos.println("Overall, hash_v2 built took: "+(endHashBuildTime-startHashBuildTime)+" ms");
            Infos.println("Words in the hash_v2: "+hash2.keySet().size());
            Infos.println("Tuples in the hash_v2:"+totalTuplesInHash);

            
            //HERE Search how many Nodes per hash bucket.
//            if (testBuckets) {
//                CustomHashMap.Node<byte[], HashPointer>[] accessToHash = hash2.getHash().getAccessToHash();
//                for (int i = 0; i < accessToHash.length; i++) {
//                    Object node = accessToHash[i];
//                    if (node instanceof core.hash.CustomHashMap.TreeNode) {
//                        CustomHashMap.TreeNode n=(CustomHashMap.TreeNode)node;
//                        System.out.println("i"+i+"=TreeNode: "+n.getValue()+" left:"+n.getLeft()+" right:"+n.getRight());
//                        System.out.println("Size:"+hash2.bucketDFS(n));
//                    } else if (node instanceof core.hash.CustomHashMap.Node) {
//                        if (node!=null)
//                            System.out.println("i"+i+"=Node: "+node.getClass().getName());
//                        else
//                            System.out.println(node);
//                    } else {
//                        System.out.println(node);
//                    }
//                }
//            }
//            fw.close();
            
            
            ///////////////////////////////////////////////////::
            ///// CHECK IF THEY CONTAIN SAME VALUES
            System.out.println("##########################################");
            
            byte[] word={1, 0, 0, 0, 3, 1, 2, 0}; //matches positions 1717,983,2149 for 26,4,1 nodes respectively
            //byte[] word={1, 1, 3, 0, 1, 3, 0, 1}; // matches 11 positions ! 388,895,898,1276,2200,2215,2253,2266,2452,2521,2581

            //byte[] word={1, 0, 0, 0, 3, 1, 2, 0};
            //byte[] word={1, 0, 0, 0, 3, 1, 2, 0};
                    
            
            //all words in hash
            //hash.keySet().stream().forEach(key->{System.out.println(key);});
            
            //number of positions per word
            //hash2.getHash().entrySet().stream().forEach(e->{System.out.println(e.getKey()+" --> "+e.getValue().getPositions().length);});
            
            System.out.println("-----");
            if (doHashV1)
                System.out.println(hash.getTuples(word).toString().replaceAll(",", "\n"));
            System.out.println("-----");
            int[] positions=hash2.getPositions(word);
            for (int p:positions) {
                System.out.println(p+":\n"+hash2.getPairs(word, p).toString().replaceAll(",", "\n"));
            }
            
            int i=0;
            for (byte[] wo:hash2.keySet()) {
                System.out.println(i+":"+hash2.getPairsOfTopPosition(wo).size());
                i++;
            }
            
            

            Infos.println("FINISHED.");
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(CustomHash_v3_Trove303.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(CustomHash_v3_Trove303.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        
        
        
    }
    
    
    
    

    

    

    
    
    
    
}
