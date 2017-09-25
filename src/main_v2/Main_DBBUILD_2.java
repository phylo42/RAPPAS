/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main_v2;

import outputs.ARProcessLauncher;
import alignement.Alignment;
import core.AAStates;
import core.ProbabilisticWord;
import core.QueryWord;
import core.States;
import core.Word;
import core.algos.AlignScoringProcess;
import core.algos.RandomSeqGenerator;
import core.algos.SequenceKnife;
import core.algos.WordExplorer;
import core.hash.CustomHash_v2;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.ARResults;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
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
import tree.PhyloTreeModel;

/**
 *
 * @author ben
 */
public class Main_DBBUILD_2 {
    
    public static final int TYPE_DNA=1;
    public static final int TYPE_PROT=1;
    
    
    /**
     * 
     * @param analysisType one of ArgumentsParser_v2.TYPE_*
     * @param processLog null if not used
     * @param k
     * @param alpha 
     * @param branchPerLength 
     * @param s states (DNA or Protein)
     * @param a alignment
     * @param t tree
     * @param workDir 
     * @param ARBinary binaries of external AR program
     * @param ARDirToUse if not null, search AR result in this directory instead of launching AR
     * @param exTreeDir if not null, serch extended tree and alignments in this directory instead of building them
     * @param buildDBFull the value of buildDBFull
     * @param forceRooting
     * @param dbInRAM
     * @param queries
     * @param callString
     * @param nsBound
     * @param noCalibration 
     * @param unionHash 
     * @throws java.io.FileNotFoundException 
     * @throws java.lang.ClassNotFoundException 
     */
    public static void DBGeneration(    int analysisType,
                                        FileWriter processLog,
                                        int k,
                                        float alpha,
                                        int branchPerLength,
                                        States s,
                                        File a,
                                        File t,
                                        File workDir,
                                        File ARBinary,
                                        File ARDirToUse,
                                        File exTreeDir,
                                        boolean buildDBFull,
                                        boolean forceRooting,
                                        boolean dbInRAM,
                                        File queries,
                                        String callString,
                                        Float nsBound,
                                        boolean noCalibration,
                                        boolean unionHash
                                    ) throws FileNotFoundException, IOException, ClassNotFoundException {
        
        

        
        
            


            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //PARAMETERS 
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            //logs
            String logPath=workDir.getAbsolutePath()+File.separator+"logs"+File.separator;
            //trees
            String extendedTreePath=workDir.getAbsolutePath()+File.separator+"extended_trees"+File.separator;
            //ancestral reconstruciton
            String ARPath=workDir.getAbsolutePath()+File.separator+"AR"+File.separator;
            
            
            //build of extended ARTree/////////////////////////////////////////////
            float minBranchLength=-1.0f;
            
            //build of AR///////////////////////////////////////////////////////
            boolean verboseAR=true;

            
            //build of Hash/////////////////////////////////////////////////////
            int knifeMode=SequenceKnife.SAMPLING_LINEAR;
            int min_k=k;
            
            float sitePPThreshold=Float.MIN_VALUE;
            float PPStarThreshold=(float)Math.pow((alpha*0.25),k);
            float PPStarThresholdAsLog=(float)Math.log10(PPStarThreshold);
            Infos.println("k="+k);
            Infos.println("factor="+alpha);
            Infos.println("PPStarThreshold="+PPStarThreshold);
            Infos.println("log10(PPStarThreshold)="+PPStarThresholdAsLog);
            //site and word posterior probas thresholds
 
            //score calibration/////////////////////////////////////////////////
            int meanCalibrationSequenceSize=150;
            int calibrationSampleSize=1000000;
            if (s instanceof AAStates) {
                meanCalibrationSequenceSize=meanCalibrationSequenceSize/3;
                calibrationSampleSize=calibrationSampleSize*10;
            }
            int q_quantile=100;
            int n_quantile=99;
            boolean writeTSVCalibrationLog=false;
                    
            //debug/////////////////////////////////////////////////////////////
            //skip extended ARTree reconstruction
            boolean buildRelaxedTree=true;
            //skip paml marginal ancestral reconstruction (made on extended ARTree)
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
            PhyloTree originalTree = NewickReader.parseNewickTree2(tline, forceRooting, false);
            Infos.println("Original tree read.");
            
            //test if alignment/tree labels are matching.
            //if not, exit before raising more errors later in the algo...
            List<String> alignLabels = Arrays.asList(align.getRowLabels());
            int notFoundCount=0;
            for (Iterator<String> iterator = alignLabels.iterator(); iterator.hasNext();) {
                String next = iterator.next();
                if (!originalTree.getLabelsByDFS().contains(next)) {
                    System.out.println("Alignment label \""+next+"\" not found in given tree labels.");
                    notFoundCount++;
                }
            }
            if (notFoundCount>0) {System.exit(1);}
            alignLabels=null;
            
            /////////////////////
            //BUILD RELAXED TREE
            
            ExtendedTree extendedTree = null;

            File fileRelaxedAlignmentFasta=new File(extendedTreePath+"extended_align.fasta");
            File fileRelaxedAlignmentPhylip=new File(extendedTreePath+"extended_align.phylip");
            File fileRelaxedTreewithBL=new File(extendedTreePath+"extended_tree_withBL.tree");
            File fileRelaxedTreewithBLNoInternalNodeLabels=new File(extendedTreePath+"extended_tree_withBL_withoutInterLabels.tree");
            File fileRelaxedTreeBinary=new File(extendedTreePath+"extended_tree.bin");
            File idsMappings=new File(extendedTreePath+"extended_tree_node_mapping.tsv");
            
            if (buildRelaxedTree) {
                if (exTreeDir==null) {
                    try {
                        System.out.println("Injecting fake nodes...");
                        //note, we read again the tree to buildDBFull a new PhyloTree object
                        //this is necessary as its TreeModel is directly modified
                        //at instanciation of ExtendedTree
                        //do a tree copy, to not do the extended tree extension on the original tree
                        PhyloNode rootCopy=originalTree.getRoot().copy();
                        PhyloTree treeCopy=new PhyloTree(new PhyloTreeModel(rootCopy),originalTree.isRooted(), false);
                        treeCopy.initIndexes();
                        extendedTree=new ExtendedTree(treeCopy,minBranchLength,branchPerLength);                    
                        extendedTree.initIndexes(); // don't forget to reinit indexes !!!
                        ArrayList<PhyloNode> listOfNewFakeLeaves = extendedTree.getFakeLeaves();
                        Infos.println("RelaxedTree contains "+extendedTree.getLeavesCount()+ " leaves");
                        Infos.println("RelaxedTree contains "+extendedTree.getFakeLeaves().size()+ " FAKE_X new leaves");
                        //add new leaves to alignment
                        char[] gapSeq=new char[align.getLength()];
                        Arrays.fill(gapSeq, '-');
                        ArrayList<char[]> seqs=new ArrayList<>();
                        String[] labels=new String[listOfNewFakeLeaves.size()];
                        for (int i = 0; i < listOfNewFakeLeaves.size(); i++) {
                            labels[i]=listOfNewFakeLeaves.get(i).getLabel();
                            seqs.add(gapSeq);
                        }
                        align.addAllSequences(labels,seqs);
                        //write alignment and ARTree for BrB
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
                        nw.close();
                        //save this extendedTree as a binary
                        FileOutputStream fos = new FileOutputStream(fileRelaxedTreeBinary);
                        ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,4096));
                        Infos.println("Storing binary version of Extended Tree.");
                        oos.writeObject(extendedTree);
                        oos.close();
                        fos.close();
                        //finally, for debugging, output the ids mappings
                        FileWriter fw=new FileWriter(idsMappings);
                        fw.append("original_id\toriginal_name\textended_id\textended_name");
                        LinkedHashMap<Integer, Integer> fakeNodeMapping = extendedTree.getFakeNodeMapping();
                        for (Iterator<Integer> iterator = fakeNodeMapping.keySet().iterator(); iterator.hasNext();) {
                            Integer next = iterator.next();
                            fw.append("\n");
                            fw.append(fakeNodeMapping.get(next)+"\t"+originalTree.getById(fakeNodeMapping.get(next)).getLabel()+"\t");
                            fw.append(next+"\t"+extendedTree.getById(next).getLabel());
                        }
                        fw.close();                        
                        
                    } catch (IOException ex) {
                        ex.printStackTrace();
                        System.out.println("Error raised from extended tree reconstruciton!");
                    }
                } else {
                    fileRelaxedAlignmentFasta=new File(exTreeDir.getAbsolutePath()+File.separator+"extended_align.fasta");
                    fileRelaxedAlignmentPhylip=new File(exTreeDir.getAbsolutePath()+File.separator+"extended_align.phylip");
                    fileRelaxedTreewithBL=new File(exTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL.tree");
                    fileRelaxedTreewithBLNoInternalNodeLabels=new File(exTreeDir.getAbsolutePath()+File.separator+"extended_tree_withBL_withoutInterLabels.tree");
                    fileRelaxedTreeBinary=new File(exTreeDir.getAbsolutePath()+File.separator+"extended_tree.bin");
                    FileInputStream fis = new FileInputStream(fileRelaxedTreeBinary);
                    ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(fis,4096));
                    Infos.println("Loading Extended Tree from binary file.");
                    extendedTree = (ExtendedTree)ois.readObject();
                    ois.close();
                    fis.close();
                    //simple test
                    if (extendedTree.getFakeInternalNodes().size()<1) {
                        System.out.println("Something went wrong with load of ExtendedTree from "+fileRelaxedTreeBinary.getAbsolutePath());
                        System.exit(1);
                    }

                }
            }
            
            
            
            //////////////////////////////////////
            //HERE LAUNCH AR ON RELAXED TREE THROUGH EXTERNAL BINARIES
            ARProcessLauncher arpl=new ARProcessLauncher(ARBinary,verboseAR,analysisType);
            //basique recognition of AR software is done through its ARBinary name
            //Note: even if AR is skipped (option --ardir),
            //ArgumentsParser.ARBinary value is used, which allows 
            //instanciation here. It will just not be executed.
                 
            if (launchAR) {
                File alignmentFile=null;
                File treeFile=null;
                if (buildRelaxedTree) {
                    alignmentFile=fileRelaxedAlignmentPhylip;
                    treeFile=fileRelaxedTreewithBLNoInternalNodeLabels;
                } else {
                    alignmentFile=a;
                    treeFile=t;
                } 
                if (ARDirToUse==null) {
                    System.out.println("Launching ancestral reconstruction...");
                    arpl.launchAR(new File(ARPath),alignmentFile, treeFile);
                } else {
                    System.out.println("Ancestral reconstruction loaded from directory set with --arpath.");
                    arpl.loadExistingAR(ARDirToUse, alignmentFile, treeFile);
                    ARPath=ARDirToUse.getAbsolutePath();
                }
            } 

            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////

            
            
            ////////////////////////////////////////////////////////////////////
            //LOAD THE NEW POSTERIOR PROBAS AND PAML TREE MADE FROM THE AR

            SessionNext_v2 session=new SessionNext_v2(k, min_k, alpha, branchPerLength, sitePPThreshold, PPStarThreshold,PPStarThresholdAsLog);
            
            Infos.println("Loading AR modified dataset (modified tree, modified alignment, Posterior Probas...)");
            ARResults arpr=null;
            //the ARResults object will take in charge 2 operations:
            //1. the call of the wrapper used to parse the AR results
            //   these will return ARTree and Posterior Probas
            //2. the correspondance between the original ARTree node names
            //   and the modification operated by the AR software (which
            //   renames internal nodes/labels in its own way...)
            System.out.println("Parsing Ancestral reconstruction results...");
            arpr=new ARResults(    arpl,
                                        align,
                                        originalTree,
                                        extendedTree,
                                        s
                                    );
            session.associateStates(s);
            session.associateInputs(arpr);
            
            //output in the AR directory the mapping of the nodes for debugging
            File map=new File(arpl.ARPath.getAbsolutePath()+File.separator+"ARtree_id_mapping.tsv");
            FileWriter fw=new FileWriter(map);
            fw.append("extended_id\textended_label\tARTree_id\tARtree_label");
            LinkedHashMap<Integer, Integer> map2 = session.ARTree.mapNodes(extendedTree);
            for (Iterator<Integer> iterator = map2.keySet().iterator(); iterator.hasNext();) {
                Integer ARTreeId = iterator.next();
                fw.append("\n");
                fw.append(map2.get(ARTreeId)+"\t"+session.extendedTree.getById(map2.get(ARTreeId)).getLabel()+"\t");
                fw.append(ARTreeId+"\t"+session.ARTree.getById(ARTreeId).getLabel());
            }
            fw.close();    
            
            
            
            
            Infos.println("#########STARTING SERIES OF RAPID TEST TO CONFIRM ANCESTRAL RECONSTRUCTION AND PARSING WENT FINE########");
            //to compare node mapping , output state of the original tree and extended tree
            Infos.println("OriginalTree rooted: "+originalTree.isRooted());
            Infos.println("OriginalTree # nodes: "+originalTree.getNodeCount());
            Infos.println("OriginalTree leaves: "+originalTree.getLeavesByDFS().size());
            Infos.println("OriginalTree internal nodes: "+originalTree.getInternalNodesByDFS().size());
            Infos.println("OriginalTree nodes by DFS:      "+originalTree.getNodeIdsByDFS().size());
            Infos.println("OriginalTree node names by DFS: "+originalTree.getLabelsByDFS().size());
            Infos.println("ExtendedTree rooted: "+extendedTree.isRooted());
            Infos.println("ExtendedTree # nodes: "+extendedTree.getNodeCount());
            Infos.println("ExtendedTree leaves: "+extendedTree.getLeavesByDFS().size());
            Infos.println("ExtendedTree internal nodes: "+(extendedTree.getNodeIdsByDFS().size()-extendedTree.getLeavesByDFS().size()));
            Infos.println("ExtendedTree new Fake leaves: "+extendedTree.getFakeLeaves().stream().mapToInt(n->n.getId()).toArray().length);
            Infos.println("ExtendedTree new Fake internal nodes: "+extendedTree.getFakeInternalNodes().stream().mapToInt(n->n.getId()).toArray().length);
            Infos.println("ExtendedTree nodes by DFS:      "+extendedTree.getNodeIdsByDFS().size());
            Infos.println("ExtendedTree node names by DFS: "+extendedTree.getLabelsByDFS().size());
            //Infos.println("Node mapping between ExtendedTree/OriginalTree nodes,  map(fake)=original : ("+extendedTree.getFakeNodeMapping().size()+" mappings) "+extendedTree.getFakeNodeMapping());
            Infos.println("Node mapping between ExtendedTree/OriginalTree nodes,  map(fake)=original : ("+extendedTree.getFakeNodeMapping().size()+" mappings) ");
            //to raidly check that AR ARTree was read correctly
            PhyloTree ARTree=arpr.getARTree();
            Infos.println("ARTree rooted: "+ARTree.isRooted());
            Infos.println("ARTree # nodes: "+ARTree.getNodeCount());
            Infos.println("ARTree leaves: "+ARTree.getLeavesByDFS().size());
            Infos.println("ARTree internal nodes: "+ARTree.getInternalNodesByDFS().size());
            Infos.println("ARTree nodes by DFS:      "+ARTree.getNodeIdsByDFS().size());
            Infos.println("ARTree node names by DFS: "+ARTree.getLabelsByDFS().size());
            //Infos.println("Node mapping between ARTree/ExtendedTree nodes, map(extended)=AR: ("+arpr.getTreeMapping().entrySet().size()+" mappings) "+arpr.getTreeMapping().toString());
            Infos.println("Node mapping between ARTree/ExtendedTree nodes, map(extended)=AR: ("+arpr.getTreeMapping().entrySet().size()+" mappings) ");
            //to raidly check that sorted probas are OK
            Infos.println("NodeId=0, 3 first PP:"+Arrays.deepToString(arpr.getPProbas().getPPSet(0, 0, 3)));
            Infos.println("NodeId=0, 3 first states:"+ Arrays.deepToString(arpr.getPProbas().getStateSet(0, 0, 3)));
            Infos.println("NodeId=0, 3 first statesIndexes:"+ Arrays.deepToString(arpr.getPProbas().getStateIndexSet(0, 0, 3)));
            Infos.println("#######################################################################");
            
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            // GENERATION OF ANCESTRAL WORDS
            
            //positions for which word are built
            SequenceKnife knife=new SequenceKnife(new String(align.getCharMatrix()[0]), k, k, s, knifeMode);
            int[] refPositions=knife.getMerOrder();     
            
            //prepare hash
            System.out.println("Building hash...");
            Infos.println("Word generator threshold will be:"+PPStarThresholdAsLog);
            if (unionHash) {
                    System.out.println("Union hash used.");
                    session.hash=new CustomHash_v2(k, s, CustomHash_v2.NODES_UNION);
            } else {
                    System.out.println("Positional hash used.");
                    session.hash=new CustomHash_v2(k, s, CustomHash_v2.NODES_POSITION);
            }
            
            //Word Explorer used to buildDBFull ancestral words
            //with a branch and bound approach
            Infos.println("Building all PP* probas...");
            int totalTuplesBuiltForHash=0;
            int nodeCounter=0;
            int nodeTested=ARTree.getInternalNodesByDFS().size();
            int nodeBucketSize=nodeTested/10;
            double startHashBuildTime=System.currentTimeMillis();
            for (int nodeId:session.ARTree.getInternalNodesByDFS()) {
                if (nodeCounter%nodeBucketSize==0) {
                    Infos.println("Node: "+nodeId +" ("+((0.0+nodeCounter)/nodeTested)*100.0+"%)" );
                    Infos.println("Current "+Environement.getMemoryUsage());
                }
                    //DEBUG
                    //if(nodeId!=709)
                    //    continue;
                    //DEBUG                
                
                //double startMerScanTime=System.currentTimeMillis();
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
                                            PPStarThresholdAsLog
                                        );
                    
                    for (int j = 0; j < arpr.getPProbas().getStateCount(); j++) {
                        wd.exploreWords(pos, j);
                    }
                    
                    //register the words in the hash
                    for (ProbabilisticWord w:wd.getRetainedWords()) {
                        session.hash.addTuple(w, w.getPpStarValue(), nodeId, w.getOriginalPosition());
                    }
                    totaTuplesInNode+=wd.getRetainedWords().size();
                    
                    //wd.getRetainedWords().stream().forEach((w)->System.out.println(w));
                    //Infos.println("Words in this position:"+wd.getRetainedWords().size());
                    
                    //double endScanTime=System.currentTimeMillis();
                    //Infos.println("Word search took "+(endScanTime-startScanTime)+" ms");
                    
                    wd=null;
                }
                totalTuplesBuiltForHash+=totaTuplesInNode;
                
                
                //register all words in the hash
                //Infos.println("Tuples in this node:"+totaTuplesInNode);
                //double endMerScanTime=System.currentTimeMillis();
                //Infos.println("Word generation in this node took "+(endMerScanTime-startMerScanTime)+" ms");
                //Environement.printMemoryUsageDescription();
                nodeCounter++;
                
            }

            
            Infos.println("Sorting hash components...");
            session.hash.sortData();
            
            double endHashBuildTime=System.currentTimeMillis();
            System.out.println("Hash built took: "+(endHashBuildTime-startHashBuildTime)+" ms");
            System.out.println("Words in the hash: "+session.hash.keySet().size());
            System.out.println("Tuples built for the hash:"+totalTuplesBuiltForHash);

                       
            ////////////////////////////////////////////////////////////////////
            //OUTPUT SOME STATS IN THE WORKDIR
            
            //double[] vals=hash.keySet().stream().mapToDouble(w->hash.getPairs(w).size()).toArray();
            //outputWordBucketSize(vals, 40, new File(workDir+"histogram_word_buckets_size_k"+k+"_mk"+min_k+"_f"+alpha+"_t"+PPStarThreshold+".png"),k,alpha);
            //outputWordPerNode(wordsPerNode, 40, new File(workDir+"histogram_word_per_node_k"+k+"_mk"+min_k+"_f"+alpha+"_t"+PPStarThreshold+".png"), k, alpha);
            
            //output some stats as histograms:
            if (histogramNumberPositionsPerNode && session.hash.keySet().size()>0) {
                Infos.println("Building #positions_per_word histogram...");
                double[] values=new double[session.hash.keySet().size()];
                int i=0;
                for (Iterator<Word> iterator = session.hash.keySet().iterator(); iterator.hasNext();) {
                    Word next = iterator.next();
                    values[i]=new Double(session.hash.getPositions(next).length);
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
                double[] values=new double[session.hash.keySet().size()];
                int i=0;
                for (Iterator<Word> iterator = session.hash.keySet().iterator(); iterator.hasNext();) {
                    Word next = iterator.next();
                    values[i]=new Double(session.hash.getPairsOfTopPosition(next).size());
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
            
            //keep filenames here, used even in dbInRAM mode
            File db=new File(workDir+File.separator+"DB_session_k"+k+"_a"+alpha+"_f"+branchPerLength+"_t"+session.PPStarThresholdAsLog10);
            File dbfull=new File(db.getAbsoluteFile()+".full");
            File dbmedium=new File(db.getAbsoluteFile()+".medium");
            File dbsmall=new File(db.getAbsoluteFile()+".small");
            File dbunion=new File(db.getAbsoluteFile()+".union");
            File dbsmallunion=new File(db.getAbsoluteFile()+".sunion");
            
            
            
            ////////////////////////////////////////////////////////////////////
            //OPTIONNAL: DO NOT SAVE DB TO FILES AND OPERATE DIRECTLY SOME
            //PLACEMENTS
            if (dbInRAM) {
                System.out.println("#############################");
                System.out.println("## dbInRAM mode !");
                System.out.println("## DB kept in memory and not exported in flat files");
                System.out.println("## queries: "+queries.getAbsolutePath());
                System.out.println("#############################");
                BufferedWriter bwTSVCalibration=null;
                int bufferSize=2097152; // buffer of 2mo
                //generate random sequences
                RandomSeqGenerator rs=new RandomSeqGenerator(session.states,meanCalibrationSequenceSize);
                
                if (session.hash.getHashType()==CustomHash_v2.NODES_POSITION) {
                    
                    //reduction to medium DB
                    System.out.println("Reduction to medium DB...");
                    session.hash.reduceToMediumHash();
                    System.gc();

                    Float calibrationNormScoreMedium =-1.0f;
                    if (nsBound==null) {
                        //calibration to medium DB
                        if (writeTSVCalibrationLog) {
                            bwTSVCalibration=new BufferedWriter(new FileWriter(new File(logPath+"calibration_medium.tsv")),bufferSize);
                        }
                        System.out.println("Score calibration on "+calibrationSampleSize+" random sequences (medium DB)...");
                        //generate random sequences
                        rs=new RandomSeqGenerator(session.states,meanCalibrationSequenceSize);
                        AlignScoringProcess asp=new AlignScoringProcess(session,Float.NEGATIVE_INFINITY, calibrationSampleSize);
                        //do the placement and calculate score quantiles
                        calibrationNormScoreMedium = asp.processCalibration(rs,calibrationSampleSize, null, SequenceKnife.SAMPLING_LINEAR, 0,q_quantile,n_quantile);
                        System.out.println("Score bound: "+calibrationNormScoreMedium);
                        //closes the calibration log  
                        if (writeTSVCalibrationLog){
                            bwTSVCalibration.close();
                        }
                    } else {
                        System.out.println("Using nsbound: "+nsBound.toString());
                        calibrationNormScoreMedium=nsBound;
                    }
                    //associate calibration
                    session.associateCalibrationScore(calibrationNormScoreMedium);
                    //now do placements on medium DB
                    System.out.println("Starting placement on medium DB...");
                    Main_PLACEMENT_v07 placer=new Main_PLACEMENT_v07(session,dbInRAM);
                    placer.doPlacements(queries, dbmedium, workDir, callString, nsBound);
                    //reduction to small DB
                    System.out.println("Reduction to small DB...");
                    session.hash.reducetoSmallHash_v2(100);
                    System.gc();
                    //calibration to small DB
                    //NOTE: not done, we keep medium DB calibration as the basis.
                    //now do placements on small DB
                    System.out.println("Starting placement on small DB...");
                    placer=new Main_PLACEMENT_v07(session,dbInRAM);
                    placer.doPlacements(queries, dbsmall, workDir, callString, nsBound);
                
                } else  if (session.hash.getHashType()==CustomHash_v2.NODES_UNION) {
                    
                    //calibration
                    float calibrationNormScoreUnion=Float.NEGATIVE_INFINITY;
                    if (nsBound==null) {
                        if (writeTSVCalibrationLog) {
                            bwTSVCalibration=new BufferedWriter(new FileWriter(new File(logPath+"calibration_small.tsv")),bufferSize);
                        }
                        System.out.println("Score calibration on "+calibrationSampleSize+" random sequences (union DB)...");
                        //do the placement and calculate score quantiles
                        AlignScoringProcess asp=new AlignScoringProcess(session,Float.NEGATIVE_INFINITY, calibrationSampleSize);
                        calibrationNormScoreUnion = asp.processCalibration(rs,calibrationSampleSize, null, SequenceKnife.SAMPLING_LINEAR, 0,q_quantile,n_quantile);
                        System.out.println("Score bound: "+calibrationNormScoreUnion);
                        //closes the calibration log  
                        if (writeTSVCalibrationLog){
                            bwTSVCalibration.close();
                        }
                    } else {
                        System.out.println("Using nsbound: "+nsBound.toString());
                        calibrationNormScoreUnion=nsBound;
                    }
                    //associate medium calibration
                    session.associateCalibrationScore(calibrationNormScoreUnion);
                    //now do placements on normal union DB
                    System.out.println("Starting placement on medium DB...");
                    Main_PLACEMENT_v07 placer=new Main_PLACEMENT_v07(session,dbInRAM);
                    placer.doPlacements(queries, dbsmall, workDir, callString, nsBound);
                    //reduction to small DB
                    System.out.println("Reduction to small DB...");
                    session.hash.reducetoSmallHash_v2(100);
                    System.gc();
                    //calibration to small DB
                    //NOTE: not done, we keep medium DB calibration as the basis.
                    //now do placements on small DB
                    System.out.println("Starting placement on small DB...");
                    placer=new Main_PLACEMENT_v07(session,dbInRAM);
                    placer.doPlacements(queries, dbsmall, workDir, callString, nsBound);
                    
                }
                
                System.out.println("DBINRAM OPERATIONS FINISHED.");
                return;
            }
            
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //SAVE HASH BY JAVA SERIALIZATION
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //buffer writer for optionnalTSV output
            BufferedWriter bwTSVCalibration=null;
            int bufferSize=2097152; // buffer of 2mo
            //generate random sequences for calibrations
            AlignScoringProcess asp=null;
            RandomSeqGenerator rs=new RandomSeqGenerator(session.states,meanCalibrationSequenceSize);
            
            
            
            ////////////////////////////////////////////////////////////////////
            // SAVE LARGE/MEDIUM/SMALL IF HASH BASED ON POSITION NODES
            
            if (session.hash.getHashType()==CustomHash_v2.NODES_POSITION) {
            
            
                ////////////////////////////////////////////////////////////////////
                // CALIBRATION: LARGE
                if (buildDBFull) {
                    //calibration
                    if (!noCalibration) {
                        if (writeTSVCalibrationLog) {
                            bwTSVCalibration=new BufferedWriter(new FileWriter(new File(logPath+"calibration_large.tsv")),bufferSize);
                        }
                        System.out.println("Score calibration on "+calibrationSampleSize+" random sequences (large DB)...");
                        asp=new AlignScoringProcess(session,Float.NEGATIVE_INFINITY, calibrationSampleSize);
                        //do the placement and calculate score quantiles
                        float calibrationNormScoreLarge = asp.processCalibration(rs,calibrationSampleSize, null, SequenceKnife.SAMPLING_LINEAR, 0,q_quantile,n_quantile);
                        System.out.println("Score bound: "+calibrationNormScoreLarge);
                        //closes the calibration log  
                        if (writeTSVCalibrationLog){
                            bwTSVCalibration.close();
                        }
                        //associate calibration
                        session.associateCalibrationScore(calibrationNormScoreLarge);
                    } else {
                        session.associateCalibrationScore(Float.NEGATIVE_INFINITY);
                    }
                    //store the DB
                    System.out.println("Serialization of the database (full)...");
                    session.storeHash(dbfull);
                }
                //System.out.println(ClassLayout.parseClass(hash.getClass()).toPrintable());
                //System.out.println(ClassLayout.parseClass(CustomNode.class).toPrintable());

                /////////////////////////////////////////////////:
                //SOME DEBUG TEST TO COMPARE MEDIUM/SMALL DBs
                //
    //            byte[] word={1,3,0,2,1,1,3,0};
    //            Infos.println("###########################################");
    //            Infos.println("#TEST DB FULL");
    //            Infos.println("Word: "+Arrays.toString(word));
    //            QueryWord queryWord = new QueryWord(word, 0);
    //            int[] positions=session.hash.getPositions(queryWord);
    //            Infos.println("Positions: "+Arrays.toString(positions));
    //            Infos.println("Top position: "+session.hash.getTopPosition(queryWord));
    //            Infos.println("Pairs top position: "+session.hash.getPairsOfTopPosition(queryWord));
    //            for (int i=1;i<positions.length;i++) {
    //                Infos.println("Pairs "+positions[i]+"th position: "+session.hash.getPairs(queryWord, positions[i]));
    //            }
    //            Infos.println("###########################################");



                ////////////////////////////////////////////////////////////////////
                //REDUCTION AND CALIBRATION: MEDIUM
                ////////////////////////////////////////////////////////////////////
                //1. REDUCE HASH CONTENT TO ONLY BEST POSITION ASSOCIATED TO EACH
                //KMER IN THE DATABASE
                //2. DO  PLACEMENT of N RANDOM SEQUENCES ON THE CREATED DATABASE 
                //THIS WILL BE USED TO CALCULTE QUANTILES and USE LAST QUANTILE
                //AS THE SCORE BOUND UNDER WHICH PLACEMENTS WILL NOT BE REPORTED
                //IN THE JPLACE OUPTUT.

                //reduction 
                session.hash.reduceToMediumHash();
                System.gc();
                //calibration
                float calibrationNormScoreMedium=Float.NEGATIVE_INFINITY;
                if (!noCalibration) {
                    if (writeTSVCalibrationLog) {
                        bwTSVCalibration=new BufferedWriter(new FileWriter(new File(logPath+"calibration_medium.tsv")),bufferSize);
                    }
                    System.out.println("Score calibration on "+calibrationSampleSize+" random sequences (medium DB)...");
                    asp=new AlignScoringProcess(session,Float.NEGATIVE_INFINITY, calibrationSampleSize);
                    //do the placement and calculate score quantiles
                    calibrationNormScoreMedium = asp.processCalibration(rs,calibrationSampleSize, null, SequenceKnife.SAMPLING_LINEAR, 0,q_quantile,n_quantile);
                    System.out.println("Score bound: "+calibrationNormScoreMedium);
                    //closes the calibration log  
                    if (writeTSVCalibrationLog){
                        bwTSVCalibration.close();
                    }
                }
                //associate medium calibration
                session.associateCalibrationScore(calibrationNormScoreMedium);
                //store in DB
                System.out.println("Serialization of the database (medium)...");
                session.storeHash(dbmedium);


                /////////////////////////////////////////////////:
                //SOME DEBUG TEST TO COMPARE MEDIUM/SMALL DBs
                //
    //            Infos.println("###########################################");
    //            Infos.println("#TEST DB MEDIUM");
    //            Infos.println("Word: "+Arrays.toString(word));
    //            queryWord = new QueryWord(word, 0);
    //            positions=session.hash.getPositions(queryWord);
    //            Infos.println("Positions: "+Arrays.toString(positions));
    //            Infos.println("Top position: "+session.hash.getTopPosition(queryWord));
    //            Infos.println("Pairs top position: "+session.hash.getPairsOfTopPosition(queryWord));
    //            for (int i=1;i<positions.length;i++) {
    //                Infos.println("Pairs "+positions[i]+"th position: "+session.hash.getPairs(queryWord, positions[i]));
    //            }
    //            Infos.println("###########################################");



                ////////////////////////////////////////////////////////////////////
                //REDUCTION AND CALIBRATION: SMALL
                ////////////////////////////////////////////////////////////////////
                //1. REDUCE HASH CONTENT TO ONLY BEST POSITION ASSOCIATED TO EACH
                //KMER IN THE DATABASE AND 10 NODES AT EACH POSITION
                //2. DO  PLACEMENT of N RANDOM SEQUENCES ON THE CREATED DATABASE 
                //THIS WILL BE USED TO CALCULTE QUANTILES and USE LAST QUANTILE
                //AS THE SCORE BOUND UNDER WHICH PLACEMENTS WILL NOT BE REPORTED
                //IN THE JPLACE OUPTUT.

                //reduction 
                //session.hash.reducetoSmallHash(10);
                session.hash.reducetoSmallHash_v2(100);
                System.gc();
                //calibration
                float calibrationNormScoreSmall=Float.NEGATIVE_INFINITY;
                if (!noCalibration) {
                    if (writeTSVCalibrationLog) {
                        bwTSVCalibration=new BufferedWriter(new FileWriter(new File(logPath+"calibration_small.tsv")),bufferSize);
                    }
                    System.out.println("Score calibration on "+calibrationSampleSize+" random sequences (small DB)...");
                    //do the placement and calculate score quantiles
                    asp=new AlignScoringProcess(session,Float.NEGATIVE_INFINITY, calibrationSampleSize);
                    calibrationNormScoreSmall = asp.processCalibration(rs,calibrationSampleSize, null, SequenceKnife.SAMPLING_LINEAR, 0,q_quantile,n_quantile);
                    System.out.println("Score bound: "+calibrationNormScoreSmall);
                    //closes the calibration log  
                    if (writeTSVCalibrationLog){
                        bwTSVCalibration.close();
                    }
                }
                //associate medium calibration
                session.associateCalibrationScore(calibrationNormScoreSmall);
                //store in DB
                System.out.println("Serialization of the database (small)...");
                session.storeHash(dbsmall);


                /////////////////////////////////////////////////:
                //SOME DEBUG TEST TO COMPARE MEDIUM/SMALL DBs
                //
    //            Infos.println("###########################################");
    //            Infos.println("Word: "+Arrays.toString(word));
    //            Infos.println("#TEST DB SMALL");
    //            queryWord = new QueryWord(word, 0);
    //            positions=session.hash.getPositions(queryWord);
    //            Infos.println("Positions: "+Arrays.toString(positions));
    //            Infos.println("Top position: "+session.hash.getTopPosition(queryWord));
    //            Infos.println("Pairs top position: "+session.hash.getPairsOfTopPosition(queryWord));
    //            for (int i=1;i<positions.length;i++) {
    //                Infos.println("Pairs "+positions[i]+"th position: "+session.hash.getPairs(queryWord, positions[i]));
    //            }
    //            Infos.println("###########################################");

                //serialization finished, output some log infos
                if (buildDBFull)
                    Infos.println("DB FULL: "+Environement.getFileSize(dbfull)+" Mb saved");
                Infos.println("DB MEDIUM: "+Environement.getFileSize(dbmedium)+" Mb saved");
                Infos.println("DB SMALL: "+Environement.getFileSize(dbsmall)+" Mb saved");
                System.out.println("\"Positional\" databases saved.");

                
            ////////////////////////////////////////////////////////////////////
            // CALIBRATE AND SAVE UNION HASH BASED ON UNION NODES
            
            } else  if (session.hash.getHashType()==CustomHash_v2.NODES_UNION) {
                
                //calibration
                float calibrationNormScoreUnion=Float.NEGATIVE_INFINITY;
                if (!noCalibration) {
                    if (writeTSVCalibrationLog) {
                        bwTSVCalibration=new BufferedWriter(new FileWriter(new File(logPath+"calibration_medium.tsv")),bufferSize);
                    }
                    System.out.println("Score calibration on "+calibrationSampleSize+" random sequences (normal union DB)...");
                    //do the placement and calculate score quantiles
                    asp=new AlignScoringProcess(session,Float.NEGATIVE_INFINITY, calibrationSampleSize);
                    calibrationNormScoreUnion = asp.processCalibration(rs,calibrationSampleSize, null, SequenceKnife.SAMPLING_LINEAR, 0,q_quantile,n_quantile);
                    System.out.println("Score bound: "+calibrationNormScoreUnion);
                    //closes the calibration log  
                    if (writeTSVCalibrationLog){
                        bwTSVCalibration.close();
                    }
                }
                //associate medium calibration
                session.associateCalibrationScore(calibrationNormScoreUnion);
                //store in DB
                System.out.println("Serialization of the database (normal union)...");
                session.storeHash(dbunion);
                
                //reduction 
                session.hash.reducetoSmallHash_v2(100);
                System.gc();
                //calibration
                float calibrationNormScoreSmallUnion=Float.NEGATIVE_INFINITY;
                if (!noCalibration) {
                    if (writeTSVCalibrationLog) {
                        bwTSVCalibration=new BufferedWriter(new FileWriter(new File(logPath+"calibration_small.tsv")),bufferSize);
                    }
                    System.out.println("Score calibration on "+calibrationSampleSize+" random sequences (small union DB)...");
                    //do the placement and calculate score quantiles
                    asp=new AlignScoringProcess(session,Float.NEGATIVE_INFINITY, calibrationSampleSize);
                    calibrationNormScoreSmallUnion = asp.processCalibration(rs,calibrationSampleSize, null, SequenceKnife.SAMPLING_LINEAR, 0,q_quantile,n_quantile);
                    System.out.println("Score bound: "+calibrationNormScoreSmallUnion);
                    //closes the calibration log  
                    if (writeTSVCalibrationLog){
                        bwTSVCalibration.close();
                    }
                }
                //associate medium calibration
                session.associateCalibrationScore(calibrationNormScoreSmallUnion);
                //store in DB
                System.out.println("Serialization of the database (small union)...");
                session.storeHash(dbsmallunion);

                //serialization finished, output some log infos
                Infos.println("DB UNION: "+Environement.getFileSize(dbunion)+" Mb saved");
                Infos.println("DB SMALL-UNION: "+Environement.getFileSize(dbsmallunion)+" Mb saved");
                System.out.println("\"Union\" database saved.");
            }
            
            //closing some stuff
            arpr=null;
            session=null;              
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
