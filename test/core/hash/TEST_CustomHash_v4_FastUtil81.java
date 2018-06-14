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
import core.hash.CustomHash_v4_FastUtil81;
import static core.hash.CustomHash_v4_FastUtil81.NODES_POSITION;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.PAMLWrapper;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
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
public class TEST_CustomHash_v4_FastUtil81 {
 
    
    
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
//            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD";
//            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL";
//            String a=inputsPath+"bv_refs_aln_stripped_99.5.fasta";
//            String t=inputsPath+"RAxML_result.bv_refs_aln";
            String workDir=HOME+"/Bureau/Data/test_CustomHash_v4_23052018/";
            String inputsPath=HOME+"/Bureau/Data/test_CustomHash_v4_23052018/";
            String a=inputsPath+"aln.fasta";
            String t=inputsPath+"ref_rooted2.tre";
            
            
            
            String rst=workDir+"AR/rst";
            
            
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
            Alignment align=new Alignment(s,fastas);
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
            ///// TEST ON HASH V2
//            boolean testBuckets=false;
//            
//            //prepare simplified hash
//            CustomHash_v4_FastUtil81 hash2=new CustomHash_v4_FastUtil81(k,s,NODES_POSITION);
//
//            Infos.println("Word generator threshold will be:"+wordPPStarThresholdAsLog);
//            Infos.println("Building all words probas...");
//            //Word Explorer
//            int totalTuplesInHash=0;
//            int nodeCounter=0;
//            double[] wordsPerNode=new double[tree.getInternalNodesByDFS().size()];
//            long startHashBuildTime=System.currentTimeMillis();
//            for (int nodeId:tree.getInternalNodesByDFS()) {
//                Infos.println("NodeId: "+nodeId+" "+tree.getById(nodeId).toString() );
//                    //DEBUG
//                    //if(nodeId!=709)
//                    //    continue;
//                    //DEBUG                
//                
//                double startMerScanTime=System.currentTimeMillis();
//                WordExplorer wd =null;
//                int totalTuplesPassingThreshold=0;
//                for (int pos:knife.getMerOrder()) {
//
//                    if(pos+k-1>align.getLength()-1)
//                        continue;
//                    //DEBUG
//                    //if(pos<1198 || pos>1202)
//                    //    continue;
//                    //DEBUG
//                    
//                    //System.out.println("Current align pos: "+pos +" to "+(pos+(k-1)));
//                    double startScanTime=System.currentTimeMillis();
//                    wd =new WordExplorer(   k,
//                                            pos,
//                                            nodeId,
//                                            pprobas,
//                                            wordPPStarThresholdAsLog,
//                                            false,
//                                            s
//                                        );
//                    
//                    for (int j = 0; j < pprobas.getStateCount(); j++) {
//                        wd.exploreWords(pos, j);
//                    }
//                    
//                    //register the words in the hash
//                    wd.getRetainedWords().stream().forEach((w)-> {hash2.addTuple(w.getWord(), w.getPpStarValue(), nodeId, w.getOriginalPosition());});
//                    totalTuplesPassingThreshold+=wd.getRetainedWords().size();
//                    
//                    //wd.getRetainedWords().stream().forEach((w)->System.out.println(w));
//                    //Infos.println("Words in this position:"+wd.getRetainedWords().size());
//                    
//                    //double endScanTime=System.currentTimeMillis();
//                    //Infos.println("Word search took "+(endScanTime-startScanTime)+" ms");
//                    
//                    wd=null;
//                }
//                wordsPerNode[nodeCounter]=totalTuplesPassingThreshold;
//                totalTuplesInHash+=totalTuplesPassingThreshold;
//                
//                
//                //register all words in the hash
//                //Infos.println("Tuples in this node:"+totalTuplesPassingThreshold);
//                double endMerScanTime=System.currentTimeMillis();
//                //Infos.println("Word search took "+(endMerScanTime-startMerScanTime)+" ms");
//                //Environement.printMemoryUsageDescription();
//                nodeCounter++;
//                
//            }

            
//            Infos.println("Resorting hash_v2 components...");
//            hash2.sortData();
//            
//            long endHashBuildTime=System.currentTimeMillis();
//            Infos.println("Overall, hash_v2 built took: "+(endHashBuildTime-startHashBuildTime)+" ms");
//            Infos.println("Words in the hash_v2: "+hash2.keySet().size());
//            Infos.println("Tuples in the hash_v2:"+totalTuplesInHash);

            
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
            
            //byte[] word={1, 0, 0, 0, 3, 1, 2, 0}; //matches positions 1717,983,2149 for 26,4,1 nodes respectively
            //byte[] word={1, 1, 3, 0, 1, 3, 0, 1}; // matches 11 positions ! 388,895,898,1276,2200,2215,2253,2266,2452,2521,2581

            //byte[] word={1, 3, 3, 0, 2, 0, 0, 0};
            //byte[] word={1, 0, 0, 0, 3, 1, 2, 0};
            byte[] word={0, 1, 0, 0, 2, 2, 3, 3};
                    
            
            //all words in hash
            //hash2.keySet().stream().forEach(key->{System.out.println(Arrays.toString(key));});
            //hash2.getHash().entrySet().stream().forEach(e->{System.out.println(Arrays.toString(e.getKey())+" --> "+e.getValue().getPairCountInTopPosition()+" --> "+e.getValue().map);});
            //hash2.getHash().entrySet().stream().forEach(e->{System.out.println(Arrays.toString(e.getKey())+" --> "+e.getValue().getNode().length+" --> "+e.getValue().list);});
            //hash2.getHash().entrySet().stream().forEach(e->{System.out.println(Arrays.toString(e.getKey())+" --> "+e.getValue().+" --> "+e.getValue().map);});
            
            //number of positions per word
            //hash2.getHash().entrySet().stream().forEach(e->{System.out.println(Arrays.toString(e.getKey())+" --> "+e.getValue().getPositions().length);});
            
//            System.out.println("-----");
//            int[] positions=hash2.getPositions(word);
//            System.out.println("positions: "+positions);
//            for (int p:positions) {
//                System.out.println(p+":\n"+hash2.getPairs(word, p).toString().replaceAll(",", "\n"));
//            }
            
//            int i=0;
//            for (byte[] wo:hash2.keySet()) {
//                System.out.println(i+":"+hash2.getPairsOfTopPosition2(wo).size());
//                i++;
//            }
            
            CustomHash_v4_FastUtil81 hash3=new CustomHash_v4_FastUtil81(k,s,NODES_POSITION);
            UnionPointerWithMap hp=new UnionPointerWithMap();
            //hp.registerTuple(155, 18, 0.00009f);
            //hash3.hash.put(word, hp);
            byte[] word2={1, 0, 0, 0, 3, 1, 2, 0};
            //hash3.addTuple(word2, 0.001f, 2, 2);
            
//            System.out.println("hash3: "+hash3.getPairs(word2));
            
            byte[] word3={1, 0, 0, 0, 3, 1, 2, 1};
            byte[] word4={1, 1, 3, 0, 1, 3, 0, 1};
            byte[] word5={1, 0, 0, 0, 3, 1, 2, 0};
            
            hash3.addTuple(word, 0.00004f, 1, 16);
            hash3.addTuple(word2, 0.00018f, 2, 25);
            hash3.addTuple(word3, 0.000000007f, 2, 5);
            hash3.addTuple(word4, 0.8f, 100, 250);
            hash3.addTuple(word5, 0.0005f, 3, 20);
            hash3.addTuple(word2, 0.001817f, 3, 25);
            hash3.addTuple(word2, 0.00015f, 3, 200);
            hash3.addTuple(word2, 0.00000798f, 3, 25);
            hash3.addTuple(word2, 0.0008f, 3, 5);
            hash3.addTuple(word2, 0.0000005f, 7, 25);
            hash3.addTuple(word4, 0.00000798f, 3, 39);
            System.out.println("-----");
//            int[] positions=hash3.getPositions(word);
//            for (int p:positions) {
//                System.out.println(p+":\n"+hash3.getPairs(word, p).toString().replaceAll(",", "\n"));
//            }
            System.out.println("word2: "+hash3.getPairs(word2));
            //System.out.println("hash:"+hash3.hash.object2ObjectEntrySet().toString());
            
            //hash3.getHash().entrySet().stream().forEach(e->{System.out.println(Arrays.toString(e.getKey())+" --> "+e.getValue().getPairCountInTopPosition()+" --> "+e.getValue().map);});
            System.out.println("word2: "+hash3.getPairs(word2)+" --> "+hash3.hash.get(word2).getPairCountInTopPosition());
            System.out.println("word2: "+hash3.getPositions(word2).length);

            Infos.println("FINISHED.");
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(CustomHash_v4_FastUtil81.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(CustomHash_v4_FastUtil81.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        
        
        
    }
    
    
    
}
