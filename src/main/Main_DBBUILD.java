/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import alignement.Alignment;
import core.AAStates;
import core.DNAStates;
import core.PProbasSorted;
import core.SimpleHash;
import core.States;
import core.algos.SequenceKnife;
import core.algos.WordExplorer;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.InputManagerNext;
import inputs.PAMLWrapper;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static main.Main_PLACEMENT.MEMORY_LARGE;
import static main.Main_PLACEMENT.TYPE_DNA;
import static main.Main_PLACEMENT.inputStreamToOutputStream;
import tree.NewickReader;
import tree.NewickWriter;
import tree.PhyloNode;
import tree.PhyloTree;
import tree.ExtendedTree;

/**
 *
 * @author ben
 */
public class Main_DBBUILD {
    
    public static final int TYPE_DNA=1;
    public static final int TYPE_PROT=1;
    
    
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose", "1");
        
        try {
            
            // INPUT FILES//////////////////////////////////////////////////////
            //here,pplacer benchmark
            String inputsPath="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/";
//            String a=wd+"bv_refs_aln.fasta";
            String a=inputsPath+"bv_refs_aln_stripped_99.5.fasta";
            String t=inputsPath+"RAxML_result.bv_refs_aln";

//            a=inputsPath+"bv_refs_aln_stripped_99.5_SMALL_SUBSET.fasta";
//            t=inputsPath+"RAxML_result.bv_refs_aln_SMALL_SUBSET.tree";

            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //PARAMETERS 
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            
            //type of Analysis//////////////////////////////////////////////////
            States s=null; 
            int analysisType=TYPE_DNA;
            
            //States: DNA or AA
            if (analysisType==TYPE_DNA)
                s=new DNAStates();   
            else if (analysisType==TYPE_DNA)
                s=new AAStates();
            
            //base path for outputs/////////////////////////////////////////////
            String path="/media/ben/STOCK/SOURCES/NetBeansProjects/ViromePlacer/WD/";
            //logs
            String logPath=path+"logs/";
            //trees
            String extendedTreePath=path+"extended_trees/";
            //ancestral reconstruciton
            String ARPath=path+"AR/";
            
            
            //build of extended tree/////////////////////////////////////////////
            float minBranchLength=0.001f;
            int numberOfFakeBranchesPerEdge=1;
            String baseMLBinaries="/media/ben/STOCK/SOFTWARE/paml4.9b_hacked/bin/baseml";
            String codeMLBinaries="/media/ben/STOCK/SOFTWARE/paml4.9b_hacked/bin/codeml";
             
            
            //build of AR///////////////////////////////////////////////////////
            boolean verboseAR=true;

            
            //build of Hash/////////////////////////////////////////////////////
            int knifeMode=SequenceKnife.SAMPLING_LINEAR;
            //mers size and word proba thresholds
            int k=5;
            int min_k=5;
            float thresholdFactor=10.0f;
            float sitePPThreshold=Float.MIN_VALUE;
            float wordPPStarThreshold=(float)(thresholdFactor*Math.pow(0.25,k));
            float thresholdAsLog=(float)Math.log10(wordPPStarThreshold);
            float wordAbsent=(float)Math.pow(0.25,k); // not used in hash, but for scoring queries
            Infos.println("k="+k);
            Infos.println("factor="+thresholdFactor);
            Infos.println("wordPPStarThreshold="+wordPPStarThreshold);
            Infos.println("wordPPStarThreshold(log10)="+thresholdAsLog);
            //site and word posterior probas thresholds
 

            
            //debug/////////////////////////////////////////////////////////////
            //skip extended tree reconstruction
            boolean buildRelaxedTree=true;
            //skip paml marginal ancestral reconstruction (made on extended tree)
            boolean launchAR=false;
            
            
            
            
            
            
            
            
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            //////////////////////
            //PREPARE DIRECTORIES
            if (!new File(path).exists()) {new File(path).mkdir();}
            if (!new File(logPath).exists()) {new File(logPath).mkdir();}
            if (!new File(extendedTreePath).exists()) {new File(extendedTreePath).mkdir();}
            if (!new File(ARPath).exists()) {new File(ARPath).mkdir();}
            
            //////////////////////
            //LOAD ORIGINAL ALIGNMENT
            FASTAPointer fp=new FASTAPointer(new File(a), false);
            Fasta fasta=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(fasta);
            }
            Alignment align=new Alignment(fastas);
            Infos.println(align.describeAlignment(false));
            fp.closePointer();
            
            /////////////////////
            //PARSE ORIGINAL TREE
            NewickReader np=new NewickReader();
            PhyloTree tree = np.parseNewickTree(new File(t));
            
            /////////////////////
            //BUILD RELAXED TREE
            File fileRelaxedAlignmentFasta=new File(extendedTreePath+"extended_align_BrB_minbl"+minBranchLength+"_"+numberOfFakeBranchesPerEdge+"peredge.fasta");
            File fileRelaxedAlignmentPhylip=new File(extendedTreePath+"extended_align_BrB_minbl"+minBranchLength+"_"+numberOfFakeBranchesPerEdge+"peredge.phylip");
            File fileRelaxedTreewithBL=new File(extendedTreePath+"extended_tree_BrB_minbl"+minBranchLength+"_"+numberOfFakeBranchesPerEdge+"peredge_withBL.tree");
            File fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreePath+"extended_tree_BrB_minbl"+minBranchLength+"_"+numberOfFakeBranchesPerEdge+"peredge_withBL_withoutInternalLabels.tree");;
            //String extendedTreeForJplace=null;
            ExtendedTree extendedTreeOnBranches=null;
            if (buildRelaxedTree) {
                try {
                    extendedTreeOnBranches=new ExtendedTree(tree,minBranchLength,numberOfFakeBranchesPerEdge);                    
                    extendedTreeOnBranches.initIndexes(); // don't forget to reinit indexes !!!
                    ArrayList<PhyloNode> listOfNewFakeLeaves = extendedTreeOnBranches.getListOfNewFakeLeaves();
                    Infos.println("RelaxedTree contains "+extendedTreeOnBranches.getNodeCount()+ " nodes");
                    Infos.println("RelaxedTree contains "+extendedTreeOnBranches.getLeavesCount()+ " leaves");
                    Infos.println("RelaxedTree contains "+extendedTreeOnBranches.getListOfNewFakeLeaves().size()+ " FAKE_X new leaves");
                    //add new leaves to alignment
                    for (int i = 0; i < listOfNewFakeLeaves.size(); i++) {
                        PhyloNode node = listOfNewFakeLeaves.get(i);
                        char[] gapSeq=new char[align.getLength()];
                        Arrays.fill(gapSeq, '-');
                        align.addSequence(node.getLabel(), gapSeq);
                    }
                    //write alignment and tree for BrB
                    Infos.println("Write extended alignment: "+fileRelaxedAlignmentFasta.getAbsolutePath());
                    align.writeAlignmentAsFasta(fileRelaxedAlignmentFasta);
                    Infos.println("Write extended alignment: "+fileRelaxedAlignmentPhylip.getAbsolutePath());
                    align.writeAlignmentAsPhylip(fileRelaxedAlignmentPhylip);
                    //write extended trees
                    Infos.println("Write extended newick tree: "+fileRelaxedTreewithBL.getAbsolutePath());
                    NewickWriter nw=new NewickWriter(fileRelaxedTreewithBL);
                    nw.writeNewickTree(extendedTreeOnBranches, true, true, false);
                    nw.close();
                    //write version without internal nodes labels
                    Infos.println("Write extended newick tree: "+fileRelaxedTreewithBLNoInternalNodeLabels.getAbsolutePath());
                    nw=new NewickWriter(fileRelaxedTreewithBLNoInternalNodeLabels);
                    nw.writeNewickTree(extendedTreeOnBranches, true, false, false);
                    //extendedTreeForJplace=nw.getNewickTree(extendedTreeOnBranches, true, true, true);
                    nw.close();
                } catch (IOException ex) {
                    ex.printStackTrace();
                    System.out.println("Error raised from extended tree reconstruciton!");
                }
            }
            
            //////////////////////////////////////
            //HERE LAUNCH BASEML ON RELAXED TREE
            //todo
            //can it be done whithout modifying ctl file but through cpmmand parameters ?
            File statsFromRelaxedTree=new File(ARPath+"rst");;
            if (launchAR) {
            
                StringBuilder sb=new StringBuilder();
                
                if (buildRelaxedTree) {
                    sb.append("seqfile = "+fileRelaxedAlignmentPhylip.getAbsolutePath()+"\n");
                    sb.append("treefile = "+fileRelaxedTreewithBLNoInternalNodeLabels.getAbsolutePath()+"\n");
                    
                } else {
                    sb.append("seqfile = "+a+"\n");
                    sb.append("treefile = "+t+"\n");
                }
                sb.append("outfile = "+ARPath+"paml_output"+"\n");
                sb.append("noisy = 2   * 0,1,2,3: how much rubbish on the screen\n");
                sb.append("verbose = 2   * set to 2 to output posterior proba distribution\n");
                sb.append("runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic 3: StepwiseAddition; (4,5):PerturbationNNI\n");
                sb.append("model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu\n");
                sb.append("Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n");
                sb.append("* ndata = 100\n");
                sb.append("clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n");
                sb.append("fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below\n");
                sb.append("kappa = 5  * initial or fixed kappa\n");
                sb.append("fix_alpha = 1   * 0: estimate alpha; 1: fix alpha at value below\n");
                sb.append("alpha = 0.433838   * initial or fixed alpha, 0:infinity (constant rate)\n");
                sb.append("Malpha = 0   * 1: different alpha's for genes, 0: one alpha\n");
                sb.append("ncatG = 25   * # of categories in the dG, AdG, or nparK models of rates\n");
                sb.append("nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK\n");
                sb.append("nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2\n");
                sb.append("getSE = 0   * 0: don't want them, 1: want S.E.s of estimates\n");
                sb.append("RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states\n");
                sb.append("Small_Diff = 7e-6\n");
                sb.append("cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?\n");
                sb.append("* icode = 0  * (with RateAncestor=1. try \"GC\" in data,model=4,Mgene=4)\n");
                sb.append("fix_blength = 2  * 0: ignore, -1: random, 1: initial, 2: fixed\n");
                sb.append("method = 1  * Optimization method 0: simultaneous; 1: one branch a time\n");
                       
                FileWriter fw=new FileWriter(new File(ARPath+"baseml.ctl"));
                Infos.println("Ancestral reconstruciton parameters written in: "+ARPath+"baseml.ctl");
                fw.append(sb);
                fw.close();
                
                //launch paml externally to build the posterior probas on the extended tree
                List<String> com=new ArrayList<>();
                com.add(baseMLBinaries);
                com.add(ARPath+"baseml.ctl");
                Infos.println("Ancestral reconstruct command: "+com);
                
                ProcessBuilder pb = new ProcessBuilder(com);
                //pb.environment().entrySet().stream().forEach((e) ->{ System.out.println(e.getKey()+"="+e.getValue()); });
                //env.put("VAR1", "myValue"); env.remove("OTHERVAR");
                pb.directory(new File(ARPath));                
                pb.redirectErrorStream(false);
                pb.redirectOutput(ProcessBuilder.Redirect.PIPE);
                pb.redirectInput(ProcessBuilder.Redirect.PIPE);
                Process p = pb.start();
                assert pb.redirectInput() == ProcessBuilder.Redirect.PIPE;
                assert p.getInputStream().read() == -1; 
                //redirect sdtout/stdin to files
                FileOutputStream STDOUTOutputStream=new FileOutputStream(new File(ARPath+"AR_sdtout.txt"));
                FileOutputStream STDERROutputStream=new FileOutputStream(new File(ARPath+"AR_sdterr.txt"));
                if (verboseAR)
                    inputStreamToOutputStream(new BufferedInputStream(p.getInputStream()), System.out);
                inputStreamToOutputStream(new BufferedInputStream(p.getInputStream()), STDOUTOutputStream);
                inputStreamToOutputStream(new BufferedInputStream(p.getErrorStream()), STDERROutputStream);
                Infos.println("External process operating reconstruction is logged in: "+new File(ARPath+"AR_sdtout.txt").getAbsolutePath());
                Infos.println("Launching reconstruction (go and take a coffee!) ...");
                try {
                    p.waitFor();
                    Thread.sleep(1000);
                } catch (InterruptedException ex) {
                    Logger.getLogger(Main_PLACEMENT.class.getName()).log(Level.SEVERE, null, ex);
                }

                STDOUTOutputStream.close();
                STDERROutputStream.close();
                Infos.println("Ancestral reconstruction finished.");
                
            }

            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////

            
            
            ////////////////////////////////////////////////////////////////////
            //LOAD THE NEW POSTERIOR PROBAS AND PAML TREE MADE FROM THE AR

            SessionNext session=new SessionNext(k, min_k, thresholdFactor, sitePPThreshold, wordPPStarThreshold/thresholdFactor);
            
            Infos.println("Loading final dataset (PAML tree and Posterior Probas ; alignment)...");
            InputManagerNext im=new InputManagerNext(InputManagerNext.SOURCE_PAML, fileRelaxedAlignmentFasta, null, statsFromRelaxedTree, s);
            session.associateStates(s);
            session.associateInputs(im);
            
            //positions for which word are built
            SequenceKnife knife=new SequenceKnife(new String(align.getCharMatrix()[0]), k, k, s, knifeMode);
            int[] refPositions=knife.getMerOrder();     
            
            
            //to raidly check that sorted probas are OK
            Infos.println("NodeId=0, 5 first PP:"+Arrays.deepToString(im.getPProbas().getPPSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(im.getPProbas().getStateSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first statesIndexes:"+ Arrays.deepToString(im.getPProbas().getStateIndexSet(0, 0, 5)));
            
            //prepare simplified hash
            SimpleHash hash=new SimpleHash();

            Infos.println("Word generator threshold will be:"+thresholdAsLog);
            //Word Explorer
            int totalWordsInHash=0;
            double startHashBuildTime=System.currentTimeMillis();
            for (int nodeId:session.tree.getInternalNodesByDFS()) {
                //Infos.println("NodeId: "+nodeId+" "+session.tree.getById(nodeId).toString() );
                    //DEBUG
//                    if(nodeId!=709)
//                        continue;
                    //DEBUG                
                
                double startMerScanTime=System.currentTimeMillis();
                WordExplorer wd =null;
                int totalWordsInNode=0;
                for (int pos:knife.getMerOrder()) {

                    if(pos+k-1>align.getLength()-1)
                        continue;
                    //DEBUG
//                    if(pos<1198 || pos>1202)
//                        continue;
                    //DEBUG
                    
                    //System.out.println("Current align pos: "+pos +" to "+(pos+(k-1)));
                    double startScanTime=System.currentTimeMillis();
                    wd =new WordExplorer(   k,
                                            pos,
                                            nodeId,
                                            im.getPProbas(),
                                            thresholdAsLog
                                        );
                    
                    for (int j = 0; j < im.getPProbas().getStateCount(); j++) {
                        wd.exploreWords(pos, j);
                    }
                    
                    //register the words in the hash
                    wd.getRetainedWords().stream().forEach((w)-> {hash.addTuple(w, w.getPpStarValue(), nodeId, w.getOriginalPosition());});
                    totalWordsInNode+=wd.getRetainedWords().size();
                    
                    //Infos.println("Words in this position:"+wd.getRetainedWords().size());
                    
                    //double endScanTime=System.currentTimeMillis();
                    //Infos.println("Word search took "+(endScanTime-startScanTime)+" ms");
                    
                    wd=null;
                }
                totalWordsInHash+=totalWordsInNode;
                
                
                //register all words in the hash
                Infos.println("Words in this node:"+totalWordsInNode);
                double endMerScanTime=System.currentTimeMillis();
                Infos.println("Word search took "+(endMerScanTime-startMerScanTime)+" ms");
                Environement.printMemoryUsage();
                
            }

            
            Infos.println("Resorting hash components...");
            hash.sortTuples();
            
            double endHashBuildTime=System.currentTimeMillis();
            Infos.println("Overall, hash built took: "+(endHashBuildTime-startHashBuildTime)+" ms");
            
            Infos.println("Words in the hash:"+totalWordsInHash);
            
            session.associateHash(hash);
            
            Infos.println("FINISHED.");
            
            
            ////////////////////////////////////////////////////////////////////
            //SAVE THE HASH BY JAVA SERIALIZATION
            
            session.store(new File(path+"PAML_session_params_k"+k+"_mk"+min_k+"_f"+thresholdFactor+"_t"+wordPPStarThreshold));
            im=null;
            session=null;  

            
            
            
            
            
            
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(WordExplorer.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(WordExplorer.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
    }
    
}
