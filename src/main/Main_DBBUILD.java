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
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;
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
        
        FileWriter fw=null;
        try {
            System.setProperty("debug.verbose", "1");
            File logProcess=new File("/media/ben/STOCK/DATA/viromeplacer/WD/log_process.csv");
            fw = new FileWriter(logProcess);
            fw.write("factor\tk\tmemory_gb\ttime_s\twords\ttuples_in_buckets\tDBsize_gb\n");
            fw.flush();
            
            int k=7;
            int maxK=7;
            int kIncrement=1;
            
            double factor=1.5;
            double maxFactor=1.5;
            double factorIncrement=0.1;
            
            for (int i = k; i < maxK+kIncrement; i+=kIncrement) {
                for (double j = factor; j < maxFactor+factorIncrement; j+=factorIncrement) {
                    Infos.println("#############################################");
                    Infos.println("#############################################");
                    Infos.println("#############################################");
                    Infos.println("#############################################");
                    Infos.println("# NEW LAUNCH; Config: k="+i+" factor="+j);
                    Infos.println("#############################################");
                    Infos.println("#############################################");

                    //DBGeneration(fw,i,(float)j);
                    fw.flush();
                    System.gc();

                }
                
            }
            
            
            fw.flush();
            
            
            fw.close();
            
            
            
        } catch (IOException ex) {
            Logger.getLogger(Main_DBBUILD.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(Main_DBBUILD.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
                
        
        
        
        
    }
    
    
    /**
     * 
     * @param processLog null if not used
     * @param k
     * @param alpha 
     * @param s 
     * @param alignmentFile 
     * @param treeFile 
     * @param workDir 
     */
    public static void DBGeneration(    FileWriter processLog, 
                                        int k, 
                                        float alpha, 
                                        int branchPerLengthAmount, 
                                        States s, 
                                        File a, 
                                        File t,
                                        File workDir,
                                        File pamlPath
                                    ) {
        
        

        
        
        try {
            


            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //PARAMETERS 
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            //logs
            String logPath=workDir+File.separator+"logs"+File.separator;
            //trees
            String extendedTreePath=workDir+File.separator+"extended_trees"+File.separator;
            //ancestral reconstruciton
            String ARPath=workDir+File.separator+"AR"+File.separator;
            
            
            //build of extended tree/////////////////////////////////////////////
            float minBranchLength=0.001f;
            String baseMLBinaries=pamlPath.getAbsolutePath();
            //String baseMLBinaries="/media/ben/STOCK/SOFTWARE/paml4.9b_hacked/bin/baseml";
            //String codeMLBinaries="/media/ben/STOCK/SOFTWARE/paml4.9b_hacked/bin/codeml";
             
            
            //build of AR///////////////////////////////////////////////////////
            boolean verboseAR=true;

            
            //build of Hash/////////////////////////////////////////////////////
            int knifeMode=SequenceKnife.SAMPLING_LINEAR;
            int min_k=k;
            
            float sitePPThreshold=Float.MIN_VALUE;
            //float wordPPStarThreshold=(float)(alpha*Math.pow(0.25,k));
            float wordPPStarThreshold=(float)Math.pow((alpha*0.25),k);
            float thresholdAsLog=(float)Math.log10(wordPPStarThreshold);
            Infos.println("k="+k);
            Infos.println("factor="+alpha);
            Infos.println("wordPPStarThreshold="+wordPPStarThreshold);
            Infos.println("wordPPStarThreshold(log10)="+thresholdAsLog);
            //site and word posterior probas thresholds
 

            
            //debug/////////////////////////////////////////////////////////////
            //skip extended tree reconstruction
            boolean buildRelaxedTree=true;
            //skip paml marginal ancestral reconstruction (made on extended tree)
            boolean launchAR=true;
            
            
            
            
            
            
            
            
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            //////////////////////
            //PREPARE DIRECTORIES
            if (!workDir.exists()) {workDir.mkdir();}
            if (!new File(logPath).exists()) {new File(logPath).mkdir();}
            if (!new File(extendedTreePath).exists()) {new File(extendedTreePath).mkdir();}
            if (!new File(ARPath).exists()) {new File(ARPath).mkdir();}
            
            //////////////////////
            //LOAD ORIGINAL ALIGNMENT
            FASTAPointer fp=new FASTAPointer(a, false);
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
            PhyloTree tree = np.parseNewickTree(t);
            
            /////////////////////
            //BUILD RELAXED TREE
            File fileRelaxedAlignmentFasta=new File(extendedTreePath+"extended_align_BrB_minbl"+minBranchLength+"_"+branchPerLengthAmount+"peredge.fasta");
            File fileRelaxedAlignmentPhylip=new File(extendedTreePath+"extended_align_BrB_minbl"+minBranchLength+"_"+branchPerLengthAmount+"peredge.phylip");
            File fileRelaxedTreewithBL=new File(extendedTreePath+"extended_tree_BrB_minbl"+minBranchLength+"_"+branchPerLengthAmount+"peredge_withBL.tree");
            File fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreePath+"extended_tree_BrB_minbl"+minBranchLength+"_"+branchPerLengthAmount+"peredge_withBL_withoutInternalLabels.tree");;
            //String extendedTreeForJplace=null;
            ExtendedTree extendedTreeOnBranches=null;
            
            if (buildRelaxedTree) {
                try {
                    System.out.println("Injecting fake nodes...");
                    extendedTreeOnBranches=new ExtendedTree(tree,minBranchLength,branchPerLengthAmount);                    
                    extendedTreeOnBranches.initIndexes(); // don't forget to reinit indexes !!!
                    ArrayList<PhyloNode> listOfNewFakeLeaves = extendedTreeOnBranches.getFakeLeaves();
                    Infos.println("RelaxedTree contains "+extendedTreeOnBranches.getLeavesCount()+ " leaves");
                    Infos.println("RelaxedTree contains "+extendedTreeOnBranches.getFakeLeaves().size()+ " FAKE_X new leaves");
                    //add new leaves to alignment
                    for (int i = 0; i < listOfNewFakeLeaves.size(); i++) {
                        PhyloNode node = listOfNewFakeLeaves.get(i);
                        char[] gapSeq=new char[align.getLength()];
                        Arrays.fill(gapSeq, '-');
                        align.addSequence(node.getLabel(), gapSeq);
                    }
                    //write alignment and tree for BrB
                    Infos.println("Write extended alignment (fasta): "+fileRelaxedAlignmentFasta.getAbsolutePath());
                    align.writeAlignmentAsFasta(fileRelaxedAlignmentFasta);
                    Infos.println("Write extended alignment (phylip): "+fileRelaxedAlignmentPhylip.getAbsolutePath());
                    align.writeAlignmentAsPhylip(fileRelaxedAlignmentPhylip);
                    //write extended trees
                    Infos.println("Write extended newick tree: "+fileRelaxedTreewithBL.getAbsolutePath());
                    NewickWriter nw=new NewickWriter(fileRelaxedTreewithBL);
                    nw.writeNewickTree(extendedTreeOnBranches, true, true, false);
                    nw.close();
                    //write version without internal nodes labels
                    Infos.println("Write extended newick tree with branch length: "+fileRelaxedTreewithBLNoInternalNodeLabels.getAbsolutePath());
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
                System.out.println("Launching PAML ancestral reconstruction...");
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
                Infos.println("Launching ancestral reconstruction (go and take a coffee, it mights take hours!) ...");
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

            SessionNext session=new SessionNext(k, min_k, alpha, sitePPThreshold, wordPPStarThreshold/alpha);
            
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
            
            System.out.println("Building hash...");
            Infos.println("Word generator threshold will be:"+thresholdAsLog);
            Infos.println("Building all words probas...");
            //Word Explorer
            int totalTuplesInHash=0;
            int nodeCounter=0;
            double[] wordsPerNode=new double[session.tree.getInternalNodesByDFS().size()];
            double startHashBuildTime=System.currentTimeMillis();
            for (int nodeId:session.tree.getInternalNodesByDFS()) {
                Infos.println("Node: "+session.tree.getById(nodeId).toString() );
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
                                            im.getPProbas(),
                                            thresholdAsLog
                                        );
                    
                    for (int j = 0; j < im.getPProbas().getStateCount(); j++) {
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
                Infos.println("Tuples in this node:"+totaTuplesInNode);
                double endMerScanTime=System.currentTimeMillis();
                Infos.println("Word generation in this node took "+(endMerScanTime-startMerScanTime)+" ms");
                //Environement.printMemoryUsageDescription();
                nodeCounter++;
                
            }

            
            Infos.println("Sorting hash components...");
            hash.sortTuples();
            
            double endHashBuildTime=System.currentTimeMillis();
            System.out.println("Hash built took: "+(endHashBuildTime-startHashBuildTime)+" ms");
            System.out.println("Words in the hash: "+hash.getKeys().size());
            System.out.println("Tuples in the hash:"+totalTuplesInHash);

            

            
                       
            ////////////////////////////////////////////////////////////////////
            //OUTPUT SOME STATS IN THE log directory
            
            double[] vals=hash.getKeys().stream().mapToDouble(w->hash.getTuples(w).size()).toArray();
            //outputWordBucketSize(vals, 40, new File(workDir+"histogram_word_buckets_size_k"+k+"_mk"+min_k+"_f"+alpha+"_t"+wordPPStarThreshold+".png"),k,alpha);
            //outputWordPerNode(wordsPerNode, 40, new File(workDir+"histogram_word_per_node_k"+k+"_mk"+min_k+"_f"+alpha+"_t"+wordPPStarThreshold+".png"), k, alpha);
            
            
            ////////////////////////////////////////////////////////////////////
            //SAVE THE HASH BY JAVA SERIALIZATION
            System.out.println("Serialization of the database...");
            session.associateHash(hash);
            File db=new File(workDir+File.separator+"PAML_session_params_k"+k+"_mk"+min_k+"_f"+alpha+"_t"+wordPPStarThreshold);
            session.store(db);
            im=null;
            session=null;  

            
            System.gc();
            double dbSize=Environement.getFileSize(db);
            System.out.println(dbSize+" Mb saved in "+db.getAbsolutePath());
            
            //double dbSize=0.0;
            //output the generation stats
            if (processLog!=null) {
                processLog.write(alpha+"\t"+k+"\t"+(Environement.getMemoryUsageAsMB()/1024)+"\t"+((endHashBuildTime-startHashBuildTime)/1000)+"\t"+hash.getKeys().size()+"\t"+(totalTuplesInHash/1e6)+"\t"+(dbSize/1024)+"\n");
            }
            
            
            Infos.println("FINISHED.");
            
            
            
            
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(WordExplorer.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(WordExplorer.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }
    
    
    
    private static void outputWordBucketSize(double[]value,int binNumber,File outputFile, int k, float factor) {
        HistogramDataset dataset = new HistogramDataset();
        dataset.setType(HistogramType.RELATIVE_FREQUENCY);
        dataset.addSeries("Histogram",value,binNumber);
        String plotTitle = "Word bucket sizes: k="+k+" fact="+factor; 
        String xaxis = "bucket_size";
        String yaxis = "proportion of words"; 
        PlotOrientation orientation = PlotOrientation.VERTICAL; 
        boolean show = false; 
        boolean toolTips = false;
        boolean urls = false; 
        JFreeChart chart = ChartFactory.createHistogram( plotTitle, xaxis, yaxis, 
                dataset, orientation, show, toolTips, urls);
        int width = 750;
        int height = 450; 
        try {
            ChartUtilities.saveChartAsPNG(outputFile, chart, width, height);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    private static void outputWordPerNode(double[]value,int binNumber,File outputFile, int k, float factor) {
        HistogramDataset dataset = new HistogramDataset();
        dataset.setType(HistogramType.RELATIVE_FREQUENCY);
        dataset.addSeries("Histogram",value,binNumber);
        String plotTitle = "Words generated per node: k="+k+" fact="+factor; 
        String xaxis = "# word";
        String yaxis = "proportion of nodes"; 
        PlotOrientation orientation = PlotOrientation.VERTICAL; 
        boolean show = false; 
        boolean toolTips = false;
        boolean urls = false; 
        JFreeChart chart = ChartFactory.createHistogram( plotTitle, xaxis, yaxis, 
                dataset, orientation, show, toolTips, urls);
        int width = 750;
        int height = 450; 
        try {
            ChartUtilities.saveChartAsPNG(outputFile, chart, width, height);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    
}
