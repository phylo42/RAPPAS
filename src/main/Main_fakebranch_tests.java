/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import au.com.bytecode.opencsv.CSVWriter;
import charts.ChartsForNodes;
import charts.ChartsForReads;
import core.Colmer;
import core.ColmerSet;
import core.DNAStates;
import core.Locality;
import core.Locality.Tuple;
import core.QueryWord;
import core.States;
import core.algos.ReadKnife;
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
import tree.PhyloNode;

/**
 *
 * @author ben
 */
public class Main_fakebranch_tests {

    public static void main(String[] args) {
        System.setProperty("debug.verbose", "1");
        
        //////////////////////
        //TEST PRUNING
        //baseDirs
        ArrayList<File> fakeBranchDirectories=new ArrayList<>();
        String baseDir="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/fake_branch_versions/";
        fakeBranchDirectories.add(new File(baseDir+"fake_to_47_L0.01/"));
        fakeBranchDirectories.add(new File(baseDir+"fake_to_47_L0.1/"));            
        fakeBranchDirectories.add(new File(baseDir+"fake_to_47_L0.5/"));
        fakeBranchDirectories.add(new File(baseDir+"fake_to_47_L0.5_better/"));
        //build complete file lists from the base
        ArrayList<File> fakeAlign=new ArrayList<>();
        ArrayList<File> fakeTree=new ArrayList<>();
        ArrayList<File> stats=new ArrayList<>();
        ArrayList<File> reads=new ArrayList<>();
        for (File file : fakeBranchDirectories) {
            System.out.println(file.getAbsolutePath()+File.separator+"fake_mod_mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta");
            fakeAlign.add(new File(file.getAbsolutePath()+File.separator+"fake_mod_mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta"));
            System.out.println(file.getAbsolutePath()+File.separator+"fake_RAxML_bipartitionsBranchLabels.result_alpha_RNApol_REROOTED.tree");
            fakeTree.add(new File(file.getAbsolutePath()+File.separator+"fake_RAxML_bipartitionsBranchLabels.result_alpha_RNApol_REROOTED.tree"));
            System.out.println(new File(file.getAbsolutePath()+File.separator+"rst"));
            stats.add(new File(file.getAbsolutePath()+File.separator+"rst"));
            System.out.println(new File(file.getAbsolutePath()+File.separator+"alphaTest1"));
            reads.add(new File(file.getAbsolutePath()+File.separator+"alphaTest1"));
        }
        System.out.println("###################################################");
            
        int k=8;
        
        //////////////////////
        //TEST PRUNING
        File alignmentFile=null;
        File treeFile=null;
        File probaFile=null;
        File readsFile=null;
        
//        for (int i = 0; i < fakeBranchDirectories.size(); i++) {
//            
//            //debug
//            i=3;
//            
//            alignmentFile=fakeAlign.get(i);
//            treeFile=fakeTree.get(i);
//            probaFile=stats.get(i);
//            readsFile=reads.get(i);
//            
//            States s=new DNAStates();
//            Session c=new Session(k,k,1e-6,1e-6);
//            c.associateStates(new DNAStates());
//            InputManager im=new InputManager(InputManager.SOURCE_PAML, alignmentFile, treeFile, probaFile, s);
//            c.associateInputs(im);
//            c.generateColmerSetData(false,ColmerSet.SAMPLING_LINEAR);
//            
//            //test
//            System.out.println("TUPLES########\n");
//            byte[] test0= {3, 0, 1, 1, 2, 1, 1, 3};
//            System.out.println("Node 41, word 68 "+c.cs.getRegister().get(new QueryWord(test0, 68)).get(23).getTuples());
//            byte[] test1= {0, 1, 1, 2, 1, 1, 3, 0};
//            System.out.println("Node 41, word 69 "+c.cs.getRegister().get(new QueryWord(test1, 69)).get(23).getTuples());
//            byte[] test2 = {1, 1, 2, 1, 1, 3, 0, 3};
//            System.out.println("Node 41, word 70 "+c.cs.getRegister().get(new QueryWord(test2, 70)).get(23).getTuples());
//            byte[] test3 = {1, 2, 1, 1, 3, 0, 3, 1};
//            System.out.println("Node 41, word 71 "+c.cs.getRegister().get(new QueryWord(test3, 71)).get(23).getTuples());
//            
//            c.store(new File(fakeBranchDirectories.get(i)+File.separator+"PAML_session_params_"+k+"_"+k+"_1e-6_1e-6_coronavirus_alpha"));
//            im=null;
//            c=null;
//            System.gc();
//            
//            //debug
//            //break;
//            
//        }
        

        
        
//        System.exit(0);


        
        
        //load tests
        k=8;

        //Order of pruning tests
        /*
        0:ake_to_47_L0.01
        1:ake_to_47_L0.1
        2:ake_to_47_L0.5
        3:ake_to_47_L0.5_better
        */
        int testIndex=3;
        readsFile=reads.get(testIndex);

        //final Session session=Session.load(new File("/media/ben/STOCK/DATA/viromeplacer/PAML_session_params_"+k+"_"+k+"_1e-6_1e-6_coronavirus_all"));

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

        FASTAPointer fp=new FASTAPointer(readsFile,true);
        Fasta fasta=null;
        QueryWord qw=null;
        int readCounter=0;
        int merCounter=0;
        FileWriter fw=null;
        
        // FOR EACH READ
        while ((fasta=fp.nextSequenceAsFastaObject())!=null) {


            System.out.println("################################################");
            System.out.println("Fasta: "+fasta.getHeader()+";"+fasta.getSequence());
            ReadKnife knife=new ReadKnife(fasta, session.k, session.minK, session.states,ReadKnife.SAMPLING_LINEAR);
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

            //            for (int i = 0; i < graphData[0][0].length; i++) {
            //                System.out.println(graphData[0][0][i]+","+graphData[0][1][i]+","+graphData[0][2][i]);
            //            }
            //diagonal sum of the log(PP*)
            LinkedHashMap<Integer,double[]> diagPPSums=new LinkedHashMap<>();
            //diagonal sum of the log(PP*)
            LinkedHashMap<Integer,int[]> diagSums=new LinkedHashMap<>();
            for (Iterator<Integer> it=session.tree.getInternalNodesByDFS().iterator();it.hasNext();) {
                int nodeId=it.next();
                double[] t=new double[session.align.getLength()];
                Arrays.fill(t, 0.0);
                diagPPSums.put(nodeId, t);
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
                        List<Tuple> tupleList=nextLocality.getTuplesSortedByProbas();
                        if (nodeId==23 || nodeId==36) {
                        System.out.println("    Total # tuples:"+nextLocality.getTuplesSortedByProbas().size());
                        System.out.println("    Tuples: "+nextLocality.getTuplesSortedByProbas());
                        }
                        Colmer currentColmer=null;
                        //System.out.println("seriesId: "+seriesId);
                        //System.out.println("items in the series: "+matchCollection.getSeries(seriesId).getItemCount());
                        //matchCollection.getSeries(seriesId).add(currentColmer.getStartSite(),qw.getOriginalPosition());
                        //Math.abs(Math.log10(t.getPpStar())));
                        for (Tuple t:tupleList) {
                            
                            //GRAPH DATA
                            
                            //                        System.out.println("      Current Tuple: "+t);
                            currentColmer=session.cs.getColmerById(t.getColmerId());
                            //                        System.out.println("      Plot:"+currentColmer.getStartSite()+","+ qw.getOriginalPosition()+","+ Math.log10(t.getPpStar()));
                            for (int i = 0; i < 1; i++) { //set i<session.k here for full mer dots
                                //the following mers can overlap previous coordinates which had a score
                                //this is allowed only if the PP* is superior to previous cases
                                int xIndex=(currentColmer.getStartSite()+i)+((qw.getOriginalPosition()+i)*session.align.getLength());
                            //                            System.out.println("currentTestedIndex: "+xIndex+"    value:"+graphData[nodeIndex][2][xIndex]);
                            //                            System.out.println("log10(tuple) "+Math.log10(t.getPpStar())+"   log10(currentValue) "+graphData[nodeIndex][2][xIndex]);
                                //if (Math.log10(t.getPpStar())>graphData[nodeIndex][2][xIndex]) {
                                    graphData[nodeIndex][2][xIndex]= Math.log10(t.getPpStar());
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
                fw = new FileWriter(new File ("diagSums.csv"));
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


                            //
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

  

        
        

        
        //here put the action to open the graph corresponding to the node


    }
    

}
