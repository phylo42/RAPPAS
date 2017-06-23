/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import alignement.Alignment;
import charts.ChartsForNodes;
import charts.ChartsForReads;
import core.older.Colmer;
import core.older.ColmerSet;
import core.DNAStates;
import core.Locality;
import core.QueryWord;
import core.States;
import core.algos.SequenceKnife;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.InputManager;
import java.awt.GridLayout;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.tree.DefaultTreeCellRenderer;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;
import tree.NewickReader;
import tree.NewickWriter;
import tree.PhyloNode;
import tree.PhyloTree;
import tree.ExtendedTree;

/**
 *
 * @author ben
 */
public class Main_test_relaxedTree {
    
    public static void main(String[] args) {
        
        try {
            System.setProperty("debug.verbose", "1");
            
            
            
            //////////////////////
            //BUILD of RELAXED TREE
            double minBranchLength=0.03;
            int fakeBranchPerEdge=2;
            
            
            String wd="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/relaxed_trees/relaxed_l0.03_N2/";
//            String wd="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/relaxed_trees/relaxed_l0.03_N1/";
//            String wd="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/relaxed_trees/relaxed_l0.03_N3/";
//            String wd="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/relaxed_trees/relaxed_BrN_mode/";
            String a=wd+"fake_mod_mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta";
            
            String t="(KJ473819:1.4236764,((KF843851:0.11987693,(KJ473805:1.0000005E-6,KJ473804:0.020912504):0.14271593):0.11670331,((KF843855:0.14776544,KF843858:0.484028):0.08614692,(AB781795:0.0022537275,(AB781793:0.013681999,(AB907634:0.022757683,(KF530123:0.02997262,(AB781796:0.020279085,(AF516906:0.022711186,((AF124986:0.0054304567,(EU769559:0.0062931282,(EU769560:0.0026218926,AB907631:0.012287344):0.0020486622):0.008870568):0.008256551,((AB907632:0.012757127,AB907633:0.0020902834):0.0100815585,((DQ431016:0.0056616943,DQ431014:0.0058972114):0.011965532,((AY874541:0.011802373,(KC175339:0.014499976,(JN856008:0.035546154,AF124992:0.018812226):0.011793815):0.001528329):0.004490078,(AY878324:0.0023628832,(AB907625:0.002437988,(AB781791:0.005369946,(AB907630:0.006924321,(AB907628:0.002390499,AB781792:0.00245468):0.0028974991):0.0050557973):0.004327501):1.0000005E-6):0.0059617604):0.008440271):0.0011446134):0.0022447358):0.006862936):0.081244774):0.011016391):0.005505663):0.0055800905):0.039562356):0.41650516):0.09997534):0.75645214);";
            String t_noAF516906="(KJ473819:1.4236764,((KF843851:0.11987693,(KJ473805:1.0000005E-6,KJ473804:0.020912504):0.14271593):0.11670331,((KF843855:0.14776544,KF843858:0.484028):0.08614692,(AB781795:0.0022537275,(AB781793:0.01368199" +
"9,(AB907634:0.022757683,(KF530123:0.02997262,(AB781796:0.020279085,((AF124986:0.0054304567,(EU769559:0.0062931282,(EU769560:0.0026218926,AB907631:0.012287344):0.0020486622):0.008870568):0.008256551,((AB907" +
"632:0.012757127,AB907633:0.0020902834):0.0100815585,((DQ431016:0.0056616943,DQ431014:0.0058972114):0.011965532,((AY874541:0.011802373,(KC175339:0.014499976,(JN856008:0.035546154,AF124992:0.018812226):0.011" +
"793815):0.001528329):0.004490078,(AY878324:0.0023628832,(AB907625:0.002437988,(AB781791:0.005369946,(AB907630:0.006924321,(AB907628:0.002390499,AB781792:0.00245468):0.0028974991):0.0050557973):0.004327501)" +
":1.0000005E-6):0.0059617604):0.008440271):0.0011446134):0.0022447358):0.08810771):0.011016391):0.005505663):0.0055800905):0.039562356):0.41650516):0.09997534):0.75645214);";
            String t_basic="((L1:1,L2:2):4,((L3:2,L4:2):1,L5:4))root:0;";
            String t_basic2="((L1:0.5,L2:1)I:2,L3:4)root:0;";
            String slide_examples="(AF516906:0.07242756981647756331,(KF843851:0.15174764181859606849,KJ473804:0.11465700686342394921):0.58137165331324824891,AB781795:0.09295312119421313135):0.0;";
            
            String t_used=slide_examples;
            
            
            //very basic test of extension
            PhyloTree tree = NewickReader.parseNewickTree2(t_used, false);
            
            for (int i:tree.getNodeIdsByDFS()) {
                System.out.println(tree.getById(i));
            }
            
            System.out.println(tree.getNodeCount()); 
            ExtendedTree exTree=new ExtendedTree(tree, 0.1f, 1);
            
            System.out.println(exTree.getNodeCount());
            System.out.println(exTree.originalEdges.toString().replaceAll(",", "\n"));
            System.out.println();
            System.out.println(exTree.extendedEdges.toString().replaceAll(",", "\n"));
            
//            exTree.displayTree();
//            Thread.sleep(60000);
            
            NewickWriter nw = new NewickWriter(new File(System.getenv("HOME")+"/test.tree"));
            nw.writeNewickTree(exTree, true, true, false);
            nw.close();
            

            System.exit(1);
            
            //load alignment
            FASTAPointer fp=new FASTAPointer(new File(a), false);
            Fasta f=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((f=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(f);
            }
            Alignment align=new Alignment(fastas);
            
            
            
            
            //build relaxed tree
////            RelaxedTree relaxedTreeOnNodes=new RelaxedTree(tree, RelaxedTree.BRANCHING_ON_NODE);
////            np.writeNewickTree(relaxedTreeOnNodes, new File("test_relaxed_on_nodes.tree"),true);
//            RelaxedTree relaxedTreeOnBranches=new RelaxedTree(tree,minBranchLength,fakeBranchPerEdge);
//            ArrayList<PhyloNode> listOfNewFakeLeaves = relaxedTreeOnBranches.getListOfNewFakeLeaves();
//            //relaxedTreeOnBranches.displayTree();
//            
//            //add new leaves to alignment
//            for (int i = 0; i < listOfNewFakeLeaves.size(); i++) {
//                PhyloNode node = listOfNewFakeLeaves.get(i);
//                char[] gapSeq=new char[align.getLength()];
//                Arrays.fill(gapSeq, '-');
//                align.addSequence(node.getLabel(), gapSeq);
//            }
//            
//            //write alignment and tree for BrB
//            align.writeAlignmentAsFasta(new File(wd+"relaxed_align_BrB_minbl"+minBranchLength+"_"+fakeBranchPerEdge+"peredge.fasta"));
//            align.writeAlignmentAsPhylip(new File(wd+"relaxed_align_BrB_minbl"+minBranchLength+"_"+fakeBranchPerEdge+"peredge.phylip"));
//            np.writeNewickTree(relaxedTreeOnBranches, new File(wd+"relaxed_tree_BrB_minbl"+minBranchLength+"_"+fakeBranchPerEdge+"peredge.tree"),true,false);
//            np.writeNewickTree(relaxedTreeOnBranches, new File(wd+"relaxed_tree_BrB_minbl"+minBranchLength+"_"+fakeBranchPerEdge+"peredge_withBL.tree"),true,true);
//
//            //version BrN
////            align.writeAlignmentAsFasta(new File(wd+"relaxed_align_BrN.fasta"));
////            align.writeAlignmentAsPhylip(new File(wd+"relaxed_align_BrN.phylip"));
////            np.writeNewickTree(relaxedTreeOnBranches, new File(wd+"relaxed_tree_BrN.tree"),true,false);
//            
//            
//            
//            System.exit(0);
            
            
        //////////////////////
        //TEST PLACEMENT ON RELAXED TREE
        //baseDirs
        ArrayList<File> fakeBranchDirectories=new ArrayList<>();
        String baseDir="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/relaxed_trees/";
        fakeBranchDirectories.add(new File(baseDir+"relaxed_l0.03_N1/"));
        fakeBranchDirectories.add(new File(baseDir+"relaxed_l0.03_N2/"));
        fakeBranchDirectories.add(new File(baseDir+"relaxed_l0.03_N3/"));
        
        //build complete file lists from the base
        ArrayList<File> fakeAlign=new ArrayList<>();
        ArrayList<File> fakeTree=new ArrayList<>();
        ArrayList<File> stats=new ArrayList<>();
        ArrayList<File> reads=new ArrayList<>();
        int count=0;
        for (File file : fakeBranchDirectories) {
            count++;
            fakeAlign.add(new File(file.getAbsolutePath()+File.separator+"relaxed_align_BrB_minbl0.03_"+count+"peredge.fasta"));
            fakeTree.add(new File(file.getAbsolutePath()+File.separator+"relaxed_tree_BrB_minbl0.03_"+count+"peredge.tree"));
            stats.add(new File(file.getAbsolutePath()+File.separator+"rst"));
            reads.add(new File(file.getAbsolutePath()+File.separator+"alphaTest1"));
        }
        //BrN test
//        fakeBranchDirectories.add(new File(baseDir+"relaxed_BrN_mode/"));
//        fakeAlign.add(new File(fakeBranchDirectories.get(3).getAbsolutePath()+File.separator+"relaxed_align_BrN.fasta"));
//        fakeTree.add(new File(fakeBranchDirectories.get(3).getAbsolutePath()+File.separator+"relaxed_align_BrN.ftree"));
//        stats.add(new File(fakeBranchDirectories.get(3).getAbsolutePath()+File.separator+"rst"));
//        reads.add(new File(fakeBranchDirectories.get(3).getAbsolutePath()+File.separator+"alphaTest1"));
        
        System.out.println("###################################################");
            
        int k=8;
        
        //////////////////////
        //TEST PLACEMENT ON RELAXED TREE
        File alignmentFile=null;
        File treeFile=null;
        File probaFile=null;
        File readsFile=null;
        
        ////////////////
        // !!!!!!!!!!!!!!!!
        int intervalStart=1920;
        int intervalEnd=2020;
        // !!!!!!!!!!!!!!!!
        ////////////////
        
        
//        for (int i = 0; i < fakeBranchDirectories.size(); i++) {
//            
//            //debug
//            i=1;
//            
//            Infos.println("Opening : "+fakeAlign.get(i));
//            alignmentFile=fakeAlign.get(i);
//            Infos.println("Opening : "+fakeTree.get(i));
//            treeFile=fakeTree.get(i);
//            Infos.println("Opening : "+stats.get(i));
//            probaFile=stats.get(i);
//            Infos.println("Opening : "+reads.get(i));
//            readsFile=reads.get(i);
//            
//            States s=new DNAStates();
//            Session c=new Session(k,k,1e-6,1e-6);
//            c.associateStates(new DNAStates());
//            
//            //!!!!!!!
//            // PAML remove the branch length in its tree repreentation
//            // so session.tree give a branchlengthless tree
//            //!!!!!!!
//            
//            InputManager im=new InputManager(InputManager.SOURCE_PAML, alignmentFile, treeFile, probaFile, s);
//            c.associateInputs(im);
//            c.generateColmerSetDataForInterval(false,ColmerSet.SAMPLING_LINEAR,intervalStart,intervalEnd);
//
//            c.store(new File(fakeBranchDirectories.get(i)+File.separator+"PAML_session_params_"+k+"_"+k+"_1e-6_1e-6_coronavirus_alpha"));
//            im=null;
//            c=null;
//            System.gc();
//            
//            //debug
//            break;
//            
//        }
//        
//
//        
//        
//        System.exit(0);


        
        
        //load tests
        k=8;

        //Order of pruning tests
        /*
        0:relaxed_l0.03_N1
        1:relaxed_l0.03_N2
        2:relaxed_l0.03_N3
        3:ake_to_47_L0.5_better
        */
        int testIndex=1;
        readsFile=reads.get(testIndex);


        final Session session=Session.load(new File(fakeBranchDirectories.get(testIndex ).getAbsolutePath()+"/PAML_session_params_"+k+"_"+k+"_1e-6_1e-6_coronavirus_alpha"));

        System.out.println("k: "+session.k);
        System.out.println("minK: "+session.minK);
        System.out.println("probaThreshold: "+session.stateThreshold);
        System.out.println("wordThreshold: "+session.wordThreshold);
        System.out.println("states: "+session.states.getStateCount());
        System.out.println("Alignment test: "+session.align.getLength()+" bp long");
        System.out.println("Tree test: node 1, "+session.tree.getById(1));
        System.out.print("PPStats test (site 1): ");
        for (int i = 0; i < 4; i++) {
            System.out.print(session.parsedProbas.getPP(session.tree.getById(1).getId(),1, i)+";");
        }
        System.out.println("");
        System.out.println("ColmerSet test: colmer count="+session.cs.getColmerCount());

        //@TODO pre-read the reads file to determine the maximum read size (here put manually to 101).
        //prepare datasets to plot the reads matches
        //XYSeriesCollection matchCollection = new XYSeriesCollection();
        //MatrixSeriesCollection matchCollection = new MatrixSeriesCollection();

        //map registering all XYZDatasets (1 per node)
        //in each dataset,
        //series 0 are the read matches [x,y,z]=[reference,read,PP*]
        //series 1 ... n can be reerved to display more information"
        HashMap<Integer,DefaultXYZDataset> datasetPerNode=new HashMap<>();


        HashMap<Integer,Integer> nodeIdToTableIndex =new HashMap<>();
        int seriesCounter=0;
        for (int i : session.tree.getInternalNodesByDFS()) {
            PhyloNode n=session.tree.getById(i);

            datasetPerNode.put(n.getId(), new DefaultXYZDataset());

            nodeIdToTableIndex.put(n.getId(), seriesCounter);
            System.out.println("NodeId,Series association: "+n.getId()+"("+n.getLabel()+"),"+ seriesCounter);
//            //XYSeries matchForNode = new XYSeries(seriesCounter);
//            //matchCollection.addSeries(matchForNode);
            seriesCounter++;
        }
//        System.out.println("# of XYSeries: "+seriesCounter);


        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
        int readLength=101;


        //load coronavirus read and register matches as some graphs
        //we will have a table with table[1]=X table[2]=Y and table[3]=Z
        //to fill the whole graph, we need X*Y lines in this 3 columns table
        double[][][] graphData =new double[session.tree.getInternalNodesByDFS().size()][3][session.align.getLength()*readLength];

        fp=new FASTAPointer(readsFile,true);
        Fasta fasta=null;
        QueryWord qw=null;
        int readCounter=0;
        int merCounter=0;
        FileWriter fw=null;
        
        // FOR EACH READ
        while ((fasta=fp.nextSequenceAsFastaObject())!=null) {


            System.out.println("################################################");
            System.out.println("Fasta: "+fasta.getHeader()+";"+fasta.getSequence());
            SequenceKnife knife=new SequenceKnife(fasta, session.k, session.minK, session.states,SequenceKnife.SAMPLING_LINEAR);
            //System.out.println("merOrder: "+Arrays.toString(knife.getMerOrder()));
            long startTime=System.currentTimeMillis();
            //init the Z axis (PP*) to very small values for all possible (X,Y)
            for (int i = 0; i < session.tree.getInternalNodesByDFS().size(); i++) {
                for (int line = 0; line < readLength; line++) {
                    for (int col = 0; col < session.align.getLength(); col++) {
                        //System.out.println(col+" "+line+" "+(col+(line*session.align.getLength())));
                        graphData[i][0][col+(line*session.align.getLength())]=col;
                        graphData[i][1][col+(line*session.align.getLength())]=line;
                        graphData[i][2][col+(line*session.align.getLength())]=Math.log10(1e-6);
                    }
                }
            }

//                        for (int i = 0; i < graphData[0][0].length; i++) {
//                            System.out.println(graphData[0][0][i]+","+graphData[0][1][i]+","+graphData[0][2][i]);
//                        }
            //diagonal sum of the log(PP*)
            LinkedHashMap<Integer,double[]> diagPPSums=new LinkedHashMap<>();
            //diagonal sum of the log(PP*)
            LinkedHashMap<Integer,int[]> diagSums=new LinkedHashMap<>();
            for (Iterator<Integer> it=session.tree.getInternalNodesByDFS().iterator();it.hasNext();) {
                int nodeId=it.next();
                double[] temp=new double[session.align.getLength()];
                Arrays.fill(temp, 0.0);
                diagPPSums.put(nodeId, temp);
                int[] t2=new int[session.align.getLength()];
                Arrays.fill(t2, 0);
                diagSums.put(nodeId, t2);
            }

            merCounter=0;
            
            try {
                
                //prepare ouput data
                fw = new FileWriter(new File ("merStructure.csv"));
                BufferedWriter bw = new BufferedWriter(fw);
                bw.append("node_label;node_id;PPCounter;log10(PP*)");
                bw.newLine();
                                
                while ((qw=knife.getNextWord())!=null) {//take next mer from the query
                    for (Integer nodeId : session.cs.getRegister().get(qw).keySet()) {//take next nodeId in which the word exists
                        int nodeIndex = nodeIdToTableIndex.get(nodeId);
                        Locality nextLocality = session.cs.getRegister().get(qw).get(nodeId);
                        if (nodeId==23 || nodeId==36) {
                        System.out.println(" WORD: "+qw.toString());
                        System.out.println("   Word found in node: "+session.tree.getById(nodeId));
                        System.out.println("   Node translated as series: "+nodeIndex);                    //for know, only plot best proba with circle size as log10 of the proba.
                        System.out.println("    Best tuple in this node:"+nextLocality.getTuplesSortedByProbas().get(0));
                        System.out.println("    Original Mer position: "+qw.getOriginalPosition());
                        }
                        List<Locality.Tuple> tupleList=nextLocality.getTuplesSortedByProbas();
                        if (nodeId==23 || nodeId==36) {
                        System.out.println("    Total # tuples:"+nextLocality.getTuplesSortedByProbas().size());
                        System.out.println("    Tuples: "+nextLocality.getTuplesSortedByProbas());
                        }
                        Colmer currentColmer=null;
                        //System.out.println("seriesId: "+seriesId);
                        //System.out.println("items in the series: "+matchCollection.getSeries(seriesId).getItemCount());
                        //matchCollection.getSeries(seriesId).add(currentColmer.getStartSite(),qw.getOriginalPosition());
                        //Math.abs(Math.log10(t.getPpStar())));
                        for (Locality.Tuple tuple:tupleList) {
                            
                            //GRAPH DATA
                            
                            //                        System.out.println("      Current Tuple: "+t);
                            currentColmer=session.cs.getColmerById(tuple.getColmerId());
                            //                        System.out.println("      Plot:"+currentColmer.getStartSite()+","+ qw.getOriginalPosition()+","+ Math.log10(t.getPpStar()));
                            for (int i = 0; i < 1; i++) { //set i<session.k here for full mer dots
                                //the following mers can overlap previous coordinates which had a score
                                //this is allowed only if the PP* is superior to previous cases
                                int xIndex=(currentColmer.getStartSite()+i)+((qw.getOriginalPosition()+i)*session.align.getLength());
                            //                            System.out.println("currentTestedIndex: "+xIndex+"    value:"+graphData[nodeIndex][2][xIndex]);
                            //                            System.out.println("log10(tuple) "+Math.log10(t.getPpStar())+"   log10(currentValue) "+graphData[nodeIndex][2][xIndex]);
                                //if (Math.log10(t.getPpStar())>graphData[nodeIndex][2][xIndex]) {
                                    graphData[nodeIndex][2][xIndex]= Math.log10(tuple.getPpStar());
                                //}
                            }
                            //System.out.println("operation: "+currentColmer.getStartSite()+"-"+qw.getOriginalPosition());
                            
                            //SCORING GLOBAL
//                            if (nodeId==23 || nodeId==36) {
//                            System.out.println("       Best scoring: bestPP="+t.getPpStar()+" : bestPosition="+currentColmer.getStartSite()+" : queryPosition="+qw.getOriginalPosition());
//                            }
//                            if (((currentColmer.getStartSite()-qw.getOriginalPosition())>=0) && (Math.log10(t.getPpStar())!=0)) {
//                                //System.out.println("       Set!: diagSums("+nodeId+")["+(bestColmer.getStartSite()-qw.getOriginalPosition())+"]="+(diagSums.get(nodeId)[bestColmer.getStartSite()-qw.getOriginalPosition()]+Math.log10(bestPP)));
//                                bw.append(session.tree.getById(nodeId).getLabel()+";"+nodeId+";"+merCounter+";"+Math.log10(t.getPpStar()));
//                                bw.newLine();
//                                diagPPSums.get(nodeId)[currentColmer.getStartSite()-qw.getOriginalPosition()]+= Math.log10(t.getPpStar());
//                                diagSums.get(nodeId)[currentColmer.getStartSite()-qw.getOriginalPosition()]+=1;
//                            }
                            
                        }

                        //for now, this doesn't take into account reads which are partially mapped at the beginning of the sequence
                        //SCORING BEST ONLY
                        double bestPP=nextLocality.getBestProbas(1)[0];
                        int bestColmerId=nextLocality.getBestLocalities(1)[0];
                        Colmer bestColmer=session.cs.getColmerById(bestColmerId);
                        if (nodeId==23 || nodeId==36) {
                        System.out.println("       Best scoring: bestPP="+bestPP+" : bestPosition="+bestColmer.getStartSite()+" : queryPosition="+qw.getOriginalPosition());
                        }
                        if (((bestColmer.getStartSite()-qw.getOriginalPosition())>=0) && (Math.log10(bestPP)!=0)) {
                            //System.out.println("       Set!: diagSums("+nodeId+")["+(bestColmer.getStartSite()-qw.getOriginalPosition())+"]="+(diagSums.get(nodeId)[bestColmer.getStartSite()-qw.getOriginalPosition()]+Math.log10(bestPP)));
                            bw.append(session.tree.getById(nodeId).getLabel()+";"+nodeId+";"+merCounter+";"+Math.log10(bestPP));
                            bw.newLine();
                            diagPPSums.get(nodeId)[bestColmer.getStartSite()-qw.getOriginalPosition()]+= Math.log10(bestPP);
                            diagSums.get(nodeId)[bestColmer.getStartSite()-qw.getOriginalPosition()]+=1;
                        }
                    }


                    merCounter++;
                    if (merCounter>101) { 
                        break;
                    }

                }
                
                bw.close();
                fw.close();
                
            } catch (Exception ex) {
                Logger.getLogger(Main_fakebranch_tests.class.getName()).log(Level.SEVERE, null, ex);
            }
                
                
                
                

            long endTime=System.currentTimeMillis();
            Infos.println("Matching and generating graphs for this read took: "+(endTime-startTime)+"ms");

            XYSeriesCollection dataset = new XYSeriesCollection();
            
            try {
                //prepare graph data
                File fdiag=new File ("diagSums.csv");
                fw = new FileWriter(fdiag);
                System.out.println("Diagsums details written in : "+fdiag.getAbsolutePath());
                BufferedWriter bw = new BufferedWriter(fw);
                bw.append("pos;node_label;node_id;sum(log10(PP*));sum(mers);score");
                bw.newLine();

                for (Iterator<Integer> iterator = diagPPSums.keySet().iterator(); iterator.hasNext();) {
                    Integer nodeId = iterator.next();
                    PhyloNode n=session.tree.getById(nodeId);

                    //FOR NOW, SCORING IS DONE HERE, WILL BE MOVED AS ITS OWN FUNCTIONS LATER
                    //////////////////////////////////////////////////////////////////////////////
                    System.out.println("DiagSums at node "+nodeId+" :"+Arrays.toString(diagPPSums.get(nodeId)));
                    XYSeries series1 = new XYSeries("id="+n.getId()+":label="+n.getLabel());
                    for (int i = 0; i < diagPPSums.get(nodeId).length; i++) {
                        double d=diagPPSums.get(nodeId)[i];
                        int involvedMers=diagSums.get(nodeId)[i];
                        if (d!=0.0) {

                            //all mer that are not considered are attributed a PP* of minK (ex: 1e-6)
                            //indeed min(score)= #mer . wordThreshold
                            //and    max(score)= #mer . 0 = 0
                            double score=d+(merCounter-involvedMers)*Math.log10(session.wordThreshold); 


                            bw.append(i+";"+n.getLabel()+";"+n.getId()+";"+d+";"+involvedMers+";"+score);
                            bw.newLine();

                            series1.add(i, score);

                            //series1.add(i, (-d+involvedMers*6));                 //here do a normalisation operation between 0 [worst proba] and 1[best probas]
                            //System.out.println("pos="+i+", sum(log10(PP*))="+d+" with sum(mers)= "+involvedMers+" -> score="+score);
                        } else {
                            //bw.append(0.0+";"+0+";"+0.0);
                            //bw.newLine();
                            series1.add(i, 0.0);
                        }
                    }


                    dataset.addSeries(series1);
                }
                

                bw.close();
                fw.close();

            } catch (Exception ex) {
                Logger.getLogger(Main_fakebranch_tests.class.getName()).log(Level.SEVERE, null, ex);
            }
                

            JPanel buildDiagScoreGraph = ChartsForReads.buildDiagScoreGraph(session.tree, dataset,fasta.getHeader());
            JFrame infos=new JFrame();
            infos.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            infos.setLayout(new GridLayout(1, 1));
            infos.add(buildDiagScoreGraph);
            infos.pack();
            RefineryUtilities.centerFrameOnScreen(infos);
            infos.setVisible(true);
            readCounter++;

            if (readCounter>2) {
                break;
            }

                
                
        }
                




        //after all values are registered, set all the XYZdatasets
        for (int i : session.tree.getInternalNodesByDFS()) {
            datasetPerNode.get(i).addSeries(0, graphData[nodeIdToTableIndex.get(i)]);
        }

        //add node listener

        //some moficiation to the default tree to display the probas when a node is clisked
        DefaultTreeCellRenderer renderer = (DefaultTreeCellRenderer) session.tree.getCellRenderer();
        renderer.setLeafIcon(null);
        renderer.setOpenIcon(null);
        renderer.setClosedIcon(null);
        session.tree.addTreeSelectionListener((TreeSelectionEvent e) -> {

            PhyloNode n=(PhyloNode)e.getNewLeadSelectionPath().getLastPathComponent();
            if (!n.isLeaf()) {
                JFrame infos=new JFrame();
                infos.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                infos.setLayout(new GridLayout(1, 1));
                System.out.println("NodeId:"+n.getId()+" seriesId:"+nodeIdToTableIndex.get(n.getId()));

                //XYSeries matrix = matchCollection.getSeries(nodeIdToMatriSeriesCOllection.get(n.getId()));
                //MatrixSeries matrix = queryMatch.getSeries(nodeIdToMatriSeriesCOllection.get(n.getId()));
//                for (Object o : matrix.getItems()) {
//                    System.out.println(((XYDataItem)o));
//                }

                //MatrixSeriesCollection matrix = matchCollection.getSeries(nodeIdToMatriSeriesCollection.get(n.getId()));


                //infos.add(ChartsForNodes.buildReadMatchForANode3(n,session.align,101, new XYSeriesCollection(matrix)));

                infos.add(ChartsForNodes.buildReadMatchForANode3(n,session.align,readLength, datasetPerNode.get(n.getId()),Math.log10(session.wordThreshold),0.0000001));
                infos.setSize(1024, 250);
                infos.pack();
                RefineryUtilities.centerFrameOnScreen(infos);
                infos.setVisible(true);

            }
        });

        session.tree.displayTree();
        
        NewickWriter np2=new NewickWriter(new File("phylotree.tree"));
        np2.writeNewickTree(session.tree, true, true, false);
        np2=null;
        
        
            
        } catch (Exception ex) {
            Logger.getLogger(Main_test_relaxedTree.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
    
    
}
