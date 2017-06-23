/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main_v2;

import alignement.Alignment;
import core.States;
import core.Word;
import core.algos.SequenceKnife;
import core.algos.WordExplorer;
import core.hash.SimpleHash_v2;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.ARProcessResults;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
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
public class Main_DBBUILD_2 {
    
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
            Logger.getLogger(Main_DBBUILD_2.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(Main_DBBUILD_2.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
                
        
        
        
        
    }
    
    
    /**
     * 
     * @param processLog null if not used
     * @param k
     * @param alpha 
     * @param branchPerLengthAmount 
     * @param s 
     * @param a 
     * @param t 
     * @param workDir 
     * @param ARExecutablePath 
     */
    public static void DBGeneration(    FileWriter processLog, 
                                        int k, 
                                        float alpha, 
                                        int branchPerLengthAmount, 
                                        States s, 
                                        File a, 
                                        File t,
                                        File workDir,
                                        File ARExecutablePath
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
            
            
            //build of extended ARExtendedTree/////////////////////////////////////////////
            float minBranchLength=0.001f;
            
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
            //skip extended ARExtendedTree reconstruction
            boolean buildRelaxedTree=true;
            //skip paml marginal ancestral reconstruction (made on extended ARExtendedTree)
            boolean launchAR=true;
            
            boolean histogramNumberPositionsPerNode=true;
            boolean hitsogramNumberNodesPerFirstPosition=true;
            
            
            
            
            
            
            
            
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
            //read line
            String line=null;
            String tline=null;
            BufferedReader br=new BufferedReader(new FileReader(t));
            while ((line=br.readLine())!=null) {tline=line;}
            br.close();
            PhyloTree originalTree = NewickReader.parseNewickTree2(tline, false);
            
            /////////////////////
            //BUILD RELAXED TREE
            
            ExtendedTree extendedTree = null;

            File fileRelaxedAlignmentFasta=new File(extendedTreePath+"extended_align_mbl"+minBranchLength+"_pedge"+branchPerLengthAmount+".fasta");
            File fileRelaxedAlignmentPhylip=new File(extendedTreePath+"extended_align_mbl"+minBranchLength+"_pedge"+branchPerLengthAmount+".phylip");
            File fileRelaxedTreewithBL=new File(extendedTreePath+"extended_tree_mbl"+minBranchLength+"_"+branchPerLengthAmount+"_withBL.tree");
            File fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreePath+"extended_tree_mbl"+minBranchLength+"_pedge"+branchPerLengthAmount+"_withBL_withoutInterLabels.tree");;
            //String extendedTreeForJplace=null;
            
            if (buildRelaxedTree) {
                try {
                    System.out.println("Injecting fake nodes...");
                    //note, we read again the tree to build a new PhyloTree object
                    //this is necessary as its TreeModel is directly modified
                    //at instanciation of ExtendedTree
                    extendedTree=new ExtendedTree(NewickReader.parseNewickTree2(tline, false),minBranchLength,branchPerLengthAmount);                    
                    extendedTree.initIndexes(); // don't forget to reinit indexes !!!
                    ArrayList<PhyloNode> listOfNewFakeLeaves = extendedTree.getFakeLeaves();
                    Infos.println("RelaxedTree contains "+extendedTree.getLeavesCount()+ " leaves");
                    Infos.println("RelaxedTree contains "+extendedTree.getFakeLeaves().size()+ " FAKE_X new leaves");
                    //add new leaves to alignment
                    for (int i = 0; i < listOfNewFakeLeaves.size(); i++) {
                        PhyloNode node = listOfNewFakeLeaves.get(i);
                        char[] gapSeq=new char[align.getLength()];
                        Arrays.fill(gapSeq, '-');
                        align.addSequence(node.getLabel(), gapSeq);
                    }
                    //write alignment and ARExtendedTree for BrB
                    Infos.println("Write extended alignment (fasta): "+fileRelaxedAlignmentFasta.getAbsolutePath());
                    align.writeAlignmentAsFasta(fileRelaxedAlignmentFasta);
                    Infos.println("Write extended alignment (phylip): "+fileRelaxedAlignmentPhylip.getAbsolutePath());
                    align.writeAlignmentAsPhylip(fileRelaxedAlignmentPhylip);
                    //write extended trees
                    Infos.println("Write extended newick tree: "+fileRelaxedTreewithBL.getAbsolutePath());
                    NewickWriter nw=new NewickWriter(fileRelaxedTreewithBL);
                    nw.writeNewickTree(extendedTree, true, true, false);
                    nw.close();
                    //write version without internal nodes labels
                    Infos.println("Write extended newick tree with branch length: "+fileRelaxedTreewithBLNoInternalNodeLabels.getAbsolutePath());
                    nw=new NewickWriter(fileRelaxedTreewithBLNoInternalNodeLabels);
                    nw.writeNewickTree(extendedTree, true, false, false);
                    //extendedTreeForJplace=nw.getNewickTree(extendedTree, true, true, true);
                    nw.close();
                } catch (IOException ex) {
                    ex.printStackTrace();
                    System.out.println("Error raised from extended tree reconstruciton!");
                }
            }
            
            
            
            //////////////////////////////////////
            //HERE LAUNCH AR ON RELAXED TREE THROUGH EXTERNAL BINARIES
            //for now, configure AR software based executable name
            ARProcessLauncher arpl=null;
            if (ARExecutablePath.getName().contains("phyml")) {
                arpl=new ARProcessLauncher(ARProcessLauncher.AR_PHYML,ARExecutablePath,verboseAR);
            } else if (ARExecutablePath.getName().contains("baseml")){
                arpl=new ARProcessLauncher(ARProcessLauncher.AR_PAML,ARExecutablePath,verboseAR);
            } else {
                System.out.println("AR binary could not be associated to know AR configuration...");
                System.exit(1);
            }          
            if (launchAR) {
                System.out.println("Launching ancestral reconstruction...");
                File alignmentFile=null;
                File treeFile=null;
                if (buildRelaxedTree) {
                    alignmentFile=fileRelaxedAlignmentPhylip;
                    treeFile=fileRelaxedTreewithBLNoInternalNodeLabels;

                } else {
                    alignmentFile=a;
                    treeFile=t;
                }
                arpl.launchAR(new File(ARPath),alignmentFile, treeFile);
            } else {
                System.out.println("AR will be loaded from already existing files...");
            }

            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////

            
            
            ////////////////////////////////////////////////////////////////////
            //LOAD THE NEW POSTERIOR PROBAS AND PAML TREE MADE FROM THE AR

            SessionNext_v2 session=new SessionNext_v2(k, min_k, alpha, sitePPThreshold, wordPPStarThreshold/alpha);
            
            Infos.println("Loading AR modified dataset (modified tree, modified alignment, Posterior Probas...)");
            ARProcessResults arpr=null;
            //the ARProcessResults object will take in charge 2 operations:
            //1. the call of the wrapper used to parse the AR results
            //   these will return ARExtendedTree and Posterior Probas
            //2. the correspondance between the original ARExtendedTree node names
            //   and the modification operated by the AR software (which
            //   renames internal nodes/labels in its own way...)
            arpr=new ARProcessResults(    arpl,
                                        align,
                                        originalTree,
                                        extendedTree,
                                        s,
                                        new File(ARPath)
                                    );
            session.associateStates(s);
            session.associateInputs(arpr);
            
            Infos.println("#########STARTING SERIES OF TEST TO CONFIRM AR PARSING WAS FINE########");
            //to raidly check that AR ARExtendedTree was read correctly
            PhyloTree ARExtendedTree=arpr.getARTree();
            Infos.println("ARExtendedTree # nodes: "+ARExtendedTree.getNodeCount());
            Infos.println("ARExtendedTree leaves: "+ARExtendedTree.getLeavesByDFS());
            Infos.println("ARExtendedTree internal nodes: "+ARExtendedTree.getInternalNodesByDFS());
            Infos.println("ARExtendedTree nodes by DFS:      "+ARExtendedTree.getNodeIdsByDFS());
            Infos.println("ARExtendedTree node names by DFS: "+ARExtendedTree.getLabelsByDFS());
            Infos.println("Node mapping with original ExtendedTree: ("+arpr.getTreeMapping().entrySet().size()+" mappings) "+arpr.getTreeMapping().toString());
            //to raidly check that sorted probas are OK
            Infos.println("NodeId=0, 5 first PP:"+Arrays.deepToString(arpr.getPProbas().getPPSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(arpr.getPProbas().getStateSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first statesIndexes:"+ Arrays.deepToString(arpr.getPProbas().getStateIndexSet(0, 0, 5)));
            Infos.println("#######################################################################");
            
            
            //positions for which word are built
            SequenceKnife knife=new SequenceKnife(new String(align.getCharMatrix()[0]), k, k, s, knifeMode);
            int[] refPositions=knife.getMerOrder();     
            
            
            
            

            
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            
            //prepare simplified hash
            SimpleHash_v2 hash=new SimpleHash_v2();
            
            System.out.println("Building hash...");
            Infos.println("Word generator threshold will be:"+thresholdAsLog);
            Infos.println("Building all words probas...");
            //Word Explorer
            int totalTuplesInHash=0;
            int nodeCounter=0;
            double[] wordsPerNode=new double[session.ARTree.getInternalNodesByDFS().size()];
            double startHashBuildTime=System.currentTimeMillis();
            for (int nodeId:session.ARTree.getInternalNodesByDFS()) {
                Infos.println("Node: "+session.ARTree.getById(nodeId).toString() );
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
                    wd =new WordExplorer(   k,
                                            pos,
                                            nodeId,
                                            arpr.getPProbas(),
                                            thresholdAsLog
                                        );
                    
                    for (int j = 0; j < arpr.getPProbas().getStateCount(); j++) {
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
            hash.sortData();
            
            double endHashBuildTime=System.currentTimeMillis();
            System.out.println("Hash built took: "+(endHashBuildTime-startHashBuildTime)+" ms");
            System.out.println("Words in the hash: "+hash.keySet().size());
            System.out.println("Tuples in the hash:"+totalTuplesInHash);

            

            
                       
            ////////////////////////////////////////////////////////////////////
            //OUTPUT SOME STATS IN THE log directory
            
            //double[] vals=hash.keySet().stream().mapToDouble(w->hash.getPairs(w).size()).toArray();
            //outputWordBucketSize(vals, 40, new File(workDir+"histogram_word_buckets_size_k"+k+"_mk"+min_k+"_f"+alpha+"_t"+wordPPStarThreshold+".png"),k,alpha);
            //outputWordPerNode(wordsPerNode, 40, new File(workDir+"histogram_word_per_node_k"+k+"_mk"+min_k+"_f"+alpha+"_t"+wordPPStarThreshold+".png"), k, alpha);
            
            
            ////////////////////////////////////////////////////////////////////
            //SAVE THE HASH BY JAVA SERIALIZATION
            System.out.println("Serialization of the database...");
            session.associateHash(hash);
            File db=new File(workDir+File.separator+"DB_session_k"+k+"_a"+alpha+"_t"+wordPPStarThreshold);
            File dbfull=new File(db.getAbsoluteFile()+".full");
            File dbmedium=new File(db.getAbsoluteFile()+".medium");
            File dbsmall=new File(db.getAbsoluteFile()+".small");
            session.storeFullHash(dbfull);
            //System.out.println(ClassLayout.parseClass(hash.getClass()).toPrintable());
            //System.out.println(ClassLayout.parseClass(CustomNode.class).toPrintable());
            session.storeMediumHash(dbmedium);
            session.storeSmallHash(dbsmall, 10);
            arpr=null;
            session=null;  
            System.gc();
            Infos.println("DB FULL: "+Environement.getFileSize(dbfull)+" Mb saved");
            Infos.println("DB MEDIUM: "+Environement.getFileSize(dbmedium)+" Mb saved");
            Infos.println("DB SMALL: "+Environement.getFileSize(dbsmall)+" Mb saved");


            
            ////////////////////////////
            //output some stats:
            if (histogramNumberPositionsPerNode) {
                Infos.println("Building #positions_per_word histogram...");
                double[] values=new double[hash.keySet().size()];
                int i=0;
                for (Iterator<Word> iterator = hash.keySet().iterator(); iterator.hasNext();) {
                    Word next = iterator.next();
                    values[i]=new Double(hash.getPositions(next).length);
                    i++;
                }
                //jfreechat histogram construction and output as image
                HistogramDataset dataset = new HistogramDataset();
                dataset.setType(HistogramType.RELATIVE_FREQUENCY);
                int bins=25;
                dataset.addSeries("Big",values,bins);
                String plotTitle = "#positions_per_word"; 
                String xaxis = "#positions";
                String yaxis = "proportion"; 
                PlotOrientation orientation = PlotOrientation.VERTICAL; 
                boolean show = false; 
                boolean toolTips = false;
                boolean urls = false; 
                JFreeChart chart = ChartFactory.createHistogram( plotTitle, xaxis, yaxis, 
                        dataset, orientation, show, toolTips, urls);
                int width = 500;
                int height = 300; 
                try {
                    ChartUtilities.saveChartAsPNG(new File(workDir.getAbsolutePath()+File.separator+"histogram_Npositions_per_word.png"), chart, width, height);
                } catch (IOException e) {}
            }
                
            if (hitsogramNumberNodesPerFirstPosition) {
                Infos.println("Building #nodes_per_1stposition histogram...");
                double[] values=new double[hash.keySet().size()];
                int i=0;
                for (Iterator<Word> iterator = hash.keySet().iterator(); iterator.hasNext();) {
                    Word next = iterator.next();
                    values[i]=new Double(hash.getPairsOfTopPosition(next).size());
                    i++;
                }
                //jfreechat histogram construction and output as image
                HistogramDataset dataset = new HistogramDataset();
                dataset.setType(HistogramType.RELATIVE_FREQUENCY);
                int bins=25;
                dataset.addSeries("Big",values,bins,0,500);
                String plotTitle = "#nodes_per_1stposition"; 
                String xaxis = "#nodes";
                String yaxis = "proportion"; 
                PlotOrientation orientation = PlotOrientation.VERTICAL; 
                boolean show = false; 
                boolean toolTips = false;
                boolean urls = false; 
                JFreeChart chart = ChartFactory.createHistogram( plotTitle, xaxis, yaxis, 
                        dataset, orientation, show, toolTips, urls);
                int width = 500;
                int height = 300; 
                try {
                    ChartUtilities.saveChartAsPNG(new File(workDir.getAbsolutePath()+File.separator+"histogram_Nnodes_per_1stposition.png"), chart, width, height);
                } catch (IOException e) {}
            }
            
            
            Infos.println("FINISHED.");
            
            
            
            
            
        } catch (Exception ex) {
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
