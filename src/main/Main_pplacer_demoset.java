/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import charts.ChartsForNodes;
import charts.ChartsForReads;
import core.DNAStates;
import core.PProbas;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
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
public class Main_pplacer_demoset {

    public static void main(String[] args) {
        System.setProperty("debug.verbose", "1");

        int k=8;

        
        ////////////////////////
        //TEST BASIC PPLACER DEMO
        File alignmentFile=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/bv_refs_aln.fasta");
        File treeFile=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/rst");
        File probaFile=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/rst");
        File readsFile=new File("/media/ben/STOCK/SOFTWARE/pplacer-Linux-v1.1.alpha18-2-gcb55169/fhcrc-microbiome-demo-730d268/src/p4z1r36.fasta");
        //File merDBFile=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s/pplacer_demoset_params_"+k+"_"+k+"_1e-6_1e-6");
        System.out.println("###################################################");
            

        
        States s=new DNAStates();
        Session S=new Session(k,k,1e-50,1e-50);
        S.associateStates(new DNAStates());
        System.out.println("Align:"+alignmentFile.getAbsolutePath());
        System.out.println("Tree:"+treeFile.getAbsolutePath());
        System.out.println("Probas:"+probaFile.getAbsolutePath());
        InputManager im=new InputManager(InputManager.SOURCE_PAML, alignmentFile, treeFile, probaFile, s);
        S.associateInputs(im);
        System.out.println("k: "+S.k);
        System.out.println("minK: "+S.minK);
        System.out.println("probaThreshold: "+S.stateThreshold);
        System.out.println("wordThreshold: "+S.wordThreshold);
        System.out.println("states: "+S.states.getStateCount());
        System.out.println("Alignment test: "+S.align.getLength()+" bp long");
        System.out.println("Tree test: node 1, "+S.tree.getById(1));
        System.out.print("PPStats test (site 1): ");
        for (int i = 0; i < 4; i++) {
            System.out.print(S.parsedProbas.getPP(S.tree.getById(1).getId(),1, i)+";");
        }
        System.out.println("");
        //System.out.println("ColmerSet test: colmer count="+S.cs.getColmerCount());

        //@TODO pre-read the reads file to determine the maximum read size (here put manually to 101).
        //prepare datasets to plot the reads matches
        //map registering all XYZDatasets (1 per node)
        //in each dataset,
        //series 0 are the read matches [x,y,z]=[reference,read,PP*]
        //series 1 ... n can be reerved to display more information"
//        HashMap<Integer,DefaultXYZDataset> datasetPerNode=new HashMap<>();
//        HashMap<Integer,Integer> nodeIdToTableIndex =new HashMap<>();
//        int seriesCounter=0;
//        for (int i : S.tree.getInternalNodesByDFS()) {
//            PhyloNode n=S.tree.getById(i);
//            datasetPerNode.put(n.getId(), new DefaultXYZDataset());
//            nodeIdToTableIndex.put(n.getId(), seriesCounter);
//            //System.out.println("NodeId,Series association: "+n.getId()+"("+n.getLabel()+"),"+ seriesCounter);
//            seriesCounter++;
//        }


        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
        
        int maxReadLength=1500;


        //load read and register matches as some graphs
        //we will have a table with table[1]=X table[2]=Y and table[3]=Z
        //to fill the whole graph, we need X*Y lines in this 3 columns table
        //double[][][] graphData =new double[S.tree.getInternalNodesByDFS().size()][3][S.align.getLength()*maxReadLength];

        FASTAPointer fp=new FASTAPointer(readsFile,true);
        Fasta fasta=null;
        QueryWord qw=null;
        int readCounter=0;
        int merCounter=0;
        FileWriter fw=null;
        
        // FOR EACH READ
        while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
            
            //prepare graphical stuff for the debugging
            long startTime=System.currentTimeMillis();
            
//            //init the Z axis (PP*) to very small values for all possible (X,Y)
//            for (int i = 0; i < S.tree.getInternalNodesByDFS().size(); i++) {
//                for (int line = 0; line < maxReadLength; line++) {
//                    for (int col = 0; col < S.align.getLength(); col++) {
//                        //System.out.println(col+" "+line+" "+(col+(line*session.align.getLength())));
//                        graphData[i][0][col+(line*S.align.getLength())]=col;
//                        graphData[i][1][col+(line*S.align.getLength())]=line;
//                        graphData[i][2][col+(line*S.align.getLength())]=Math.log10(1e-6);
//                    }
//                }
//            }
//            //diagonal sum of the log(PP*), fill with zeros
//            LinkedHashMap<Integer,double[]> diagPPSums=new LinkedHashMap<>();
//            LinkedHashMap<Integer,int[]> diagSums=new LinkedHashMap<>();
//            for (Iterator<Integer> it=S.tree.getInternalNodesByDFS().iterator();it.hasNext();) {
//                int nodeId=it.next();
//                double[] t=new double[S.align.getLength()];
//                Arrays.fill(t, 0.0);
//                diagPPSums.put(nodeId, t);
//                int[] t2=new int[S.align.getLength()];
//                Arrays.fill(t2, 0);
//                diagSums.put(nodeId, t2);
//            }

            try {
                
                //prepare ouput data
                fw = new FileWriter(new File ("merStructure.csv"));
                BufferedWriter bw = new BufferedWriter(fw);
                bw.append("node_label;node_id;PPCounter;log10(PP*)");
                bw.newLine();
                
                
                /////////////////////////////////////////////////////////////
                //step 1
                /////////////////////////////////////////////////////////////
                //generate SAMPLING_NON_OVERLAPPING mers from ReadKnife
                System.out.println("################################################");
                System.out.println("Fasta: "+fasta.getHeader()+";"+fasta.getSequence());
                SequenceKnife knife=new SequenceKnife(fasta, k, S.minK, S.states,SequenceKnife.SAMPLING_NON_OVERLAPPING);
                System.out.println("merOrder: "+Arrays.toString(knife.getMerOrder()));

                System.exit(1);
                
                
                merCounter=0;               
                while ((qw=knife.getNextWord())!=null) {//take next mer from the query
                    
                    /////////////////////////////////////////////////////////////
                    //step 2
                    /////////////////////////////////////////////////////////////
                    //use PPstats in one single pass to build diagsums along the reference from this basic mer sampling
                    //do that for how many nodes ? let's every "nodeStep", from nodes as returned by DFS
                    //set the best diagsums as a list of potential localities taht will be further investgated
                    //
                    
                    //word matches are calculted for every nodeStep nodes 
                    int nodeStep=3;
                    byte[] word = qw.getWord();
                    PProbas P=S.parsedProbas;
                    ArrayList<Integer> nodesByDFS = S.tree.getNodeIdsByDFS();
                    //ArrayBlockingQueue<Double> queue =new ArrayBlockingQueue(k, false);
                    
                    HashMap<Integer,ArrayList<Double>> sums=new HashMap<>(S.parsedProbas.getNodeCount());
                    
                    for (int i = 0; i < nodesByDFS.size(); i=i+nodeStep) {
                        int nodeId = nodesByDFS.get(i);
                        if (!sums.containsKey(nodeId)) {
                            sums.put(nodeId, new ArrayList<>(S.parsedProbas.getSiteCount()));
                            for (int j = 0; j < S.parsedProbas.getSiteCount(); j++) {
                                sums.get(nodeId).add(0.0);
                            }
                        }
                        
                        //this loop works only for the sampling currenlty used (ReadKnife.SAMPLING_NON_OVERLAPPING)
                        double PPStar=1.0;
                        for (int site = 0; site < P.getSiteCount(); site++) {
                            PPStar=PPStar*P.getPP(nodeId, site,word[site%k]);
                            if (site>0 && (site+1)%k==0) {
                                //PP* are registered only every k sites
                                sums.get(nodeId).set(site,PPStar);
                                PPStar=1.0;
                                
                            } else {
                                PPStar=PPStar*P.getPP(nodeId, site,word[site%k]);
                            }
                            
                        }
                         
                    }
                    
                    for (Iterator<Integer> iterator = sums.keySet().iterator(); iterator.hasNext();) {
                        Integer next = iterator.next();
                        System.out.println(sums.get(next));
                        
                    }
                    
                    
                    
                    
                    
                    
                    
//                    for (Integer nodeId : S.cs.getRegister().get(qw).keySet()) {//take next nodeId in which the word exists
//                        int nodeIndex = nodeIdToTableIndex.get(nodeId);
//                        Locality nextLocality = S.cs.getRegister().get(qw).get(nodeId);
//                        if (nodeId==23 || nodeId==36) {
//                        System.out.println(" WORD: "+qw.toString());
//                        System.out.println("   Word found in node: "+S.tree.getById(nodeId));
//                        System.out.println("   Node translated as series: "+nodeIndex);                    //for know, only plot best proba with circle size as log10 of the proba.
//                        System.out.println("    Best tuple in this node:"+nextLocality.getTuplesSortedByProbas().get(0));
//                        System.out.println("    Original Mer position: "+qw.getOriginalPosition());
//                        }
//                        List<Tuple> tupleList=nextLocality.getTuplesSortedByProbas();
//                        if (nodeId==23 || nodeId==36) {
//                        System.out.println("    Total # tuples:"+nextLocality.getTuplesSortedByProbas().size());
//                        System.out.println("    Tuples: "+nextLocality.getTuplesSortedByProbas());
//                        }
//                        Colmer currentColmer=null;
//                        //System.out.println("seriesId: "+seriesId);
//                        //System.out.println("items in the series: "+matchCollection.getSeries(seriesId).getItemCount());
//                        //matchCollection.getSeries(seriesId).add(currentColmer.getStartSite(),qw.getOriginalPosition());
//                        //Math.abs(Math.log10(t.getPpStar())));
//                        for (Tuple t:tupleList) {
//                            //System.out.println("      Current Tuple: "+t);
//                            currentColmer=S.cs.getColmerById(t.getColmerId());
//                            //System.out.println("      Plot:"+currentColmer.getStartSite()+","+ qw.getOriginalPosition()+","+ Math.log10(t.getPpStar()));
//                            for (int i = 0; i < 1; i++) { //set i<session.k here for full mer dots
//                                //the following mers can overlap previous coordinates which had a score
//                                //this is allowed only if the PP* is superior to previous cases
//                                int xIndex=(currentColmer.getStartSite()+i)+((qw.getOriginalPosition()+i)*S.align.getLength());
//                            //  System.out.println("currentTestedIndex: "+xIndex+"    value:"+graphData[nodeIndex][2][xIndex]);
//                            //  System.out.println("log10(tuple) "+Math.log10(t.getPpStar())+"   log10(currentValue) "+graphData[nodeIndex][2][xIndex]);
//                                //if (Math.log10(t.getPpStar())>graphData[nodeIndex][2][xIndex]) {
//                                    graphData[nodeIndex][2][xIndex]= Math.log10(t.getPpStar());
//                                //}
//                            }
//
//                            
//                        }
//
//                        //for now, this doesn't take into account reads which are partially mapped at the beginning of the sequence
//                        //SCORING BEST ONLY
//                        double bestPP=nextLocality.getBestProbas(1)[0];
//                        int bestColmerId=nextLocality.getBestLocalities(1)[0];
//                        Colmer bestColmer=S.cs.getColmerById(bestColmerId);
//                        if (nodeId==23 || nodeId==36) {
//                        System.out.println("       Best scoring: bestPP="+bestPP+" : bestPosition="+bestColmer.getStartSite()+" : queryPosition="+qw.getOriginalPosition());
//                        }
//                        if (((bestColmer.getStartSite()-qw.getOriginalPosition())>=0) && (Math.log10(bestPP)!=0)) {
//                            //System.out.println("       Set!: diagSums("+nodeId+")["+(bestColmer.getStartSite()-qw.getOriginalPosition())+"]="+(diagSums.get(nodeId)[bestColmer.getStartSite()-qw.getOriginalPosition()]+Math.log10(bestPP)));
//                            bw.append(S.tree.getById(nodeId).getLabel()+";"+nodeId+";"+merCounter+";"+Math.log10(bestPP));
//                            bw.newLine();
//                            diagPPSums.get(nodeId)[bestColmer.getStartSite()-qw.getOriginalPosition()]+= Math.log10(bestPP);
//                            diagSums.get(nodeId)[bestColmer.getStartSite()-qw.getOriginalPosition()]+=1;
//                        }
                        
                        
                    }
                    
                    /////////////////////////////////////////////////////////////
                    //step 3
                    /////////////////////////////////////////////////////////////
                    //focus on the potential localities
                    //build full diagsum for a region including the locality
                    //do that for all nodes, to define which is the best node.
                    //
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    


                merCounter++;
                if (merCounter>maxReadLength) { 
                    break;
                }

                
                
                bw.close();
                fw.close();
                
            } catch (Exception ex) {
                Logger.getLogger(Main_pplacer_demoset.class.getName()).log(Level.SEVERE, null, ex);
            }
                
             
            //#######################################################
            System.exit(1);
                
                

            long endTime=System.currentTimeMillis();
            Infos.println("Matching and generating graphs for this read took: "+(endTime-startTime)+"ms");

            XYSeriesCollection dataset = new XYSeriesCollection();
            
            try {
                //prepare graph data
                fw = new FileWriter(new File ("diagSums.csv"));
                BufferedWriter bw = new BufferedWriter(fw);
                bw.append("pos;node_label;node_id;sum(log10(PP*));sum(mers);score");
                bw.newLine();

//                for (Iterator<Integer> iterator = diagPPSums.keySet().iterator(); iterator.hasNext();) {
//                    Integer nodeId = iterator.next();
//                    PhyloNode n=S.tree.getById(nodeId);
//
//                    //FOR NOW, SCORING IS DONE HERE, WILL BE MOVED AS ITS OWN FUNCTIONS LATER
//                    //////////////////////////////////////////////////////////////////////////////
//                    System.out.println("DiagSums at node "+nodeId+" :"+Arrays.toString(diagPPSums.get(nodeId)));
//                    XYSeries series1 = new XYSeries("id="+n.getId()+":label="+n.getLabel());
//                    for (int i = 0; i < diagPPSums.get(nodeId).length; i++) {
//                        double d=diagPPSums.get(nodeId)[i];
//                        int involvedMers=diagSums.get(nodeId)[i];
//                        if (d!=0.0) {
//
//                            //all mer that are not considered are attributed a PP* of minK (ex: 1e-6)
//                            //indeed min(score)= #mer . wordThreshold
//                            //and    max(score)= #mer . 0 = 0
//                            double score=d+(merCounter-involvedMers)*Math.log10(S.wordThreshold); 
//
//
//                            bw.append(i+";"+n.getLabel()+";"+n.getId()+";"+d+";"+involvedMers+";"+score);
//                            bw.newLine();
//
//
//                            //
//                            series1.add(i, score);
//
//                            //series1.add(i, (-d+involvedMers*6));                 //here do a normalisation operation between 0 [worst proba] and 1[best probas]
//                            //System.out.println("pos="+i+", sum(log10(PP*))="+d+" with sum(mers)= "+involvedMers+" -> score="+score);
//                        } else {
//                            //bw.append(0.0+";"+0+";"+0.0);
//                            //bw.newLine();
//                            series1.add(i, 0.0);
//                        }
//                    }
//
//
//                    dataset.addSeries(series1);
//                }

                bw.close();
                fw.close();

            } catch (Exception ex) {
                Logger.getLogger(Main_pplacer_demoset.class.getName()).log(Level.SEVERE, null, ex);
            }
                

            JPanel buildDiagScoreGraph = ChartsForReads.buildDiagScoreGraph(S.tree, dataset,fasta.getHeader());
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
//        for (int i : S.tree.getInternalNodesByDFS()) {
//            datasetPerNode.get(i).addSeries(0, graphData[nodeIdToTableIndex.get(i)]);
//        }

        //add node listener

        //some moficiation to the default tree to display the probas when a node is clisked
        DefaultTreeCellRenderer renderer = (DefaultTreeCellRenderer) S.tree.getCellRenderer();
        renderer.setLeafIcon(null);
        renderer.setOpenIcon(null);
        renderer.setClosedIcon(null);
        S.tree.addTreeSelectionListener((TreeSelectionEvent e) -> {

            PhyloNode n=(PhyloNode)e.getNewLeadSelectionPath().getLastPathComponent();
            if (!n.isLeaf()) {
                JFrame infos=new JFrame();
                infos.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                infos.setLayout(new GridLayout(1, 1));
                //System.out.println("NodeId:"+n.getId()+" seriesId:"+nodeIdToTableIndex.get(n.getId()));


                //infos.add(ChartsForNodes.buildReadMatchForANode3(n,S.align,maxReadLength, datasetPerNode.get(n.getId()),Math.log10(S.wordThreshold),0.0000001));
                infos.setSize(1024, 250);
                infos.pack();
                RefineryUtilities.centerFrameOnScreen(infos);
                infos.setVisible(true);

            }
        });

        S.tree.displayTree();

        
        //System.exit(0);

        
//        String word="ATCCTGAC";
//        k=word.length();
//        DNAStates states=new DNAStates();
//        QueryWord qw=new QueryWord(states.getBytes(word.toCharArray()), 0);
//        System.out.println(" - matches of : "+word+" "+Arrays.toString(qw.getWord()));
//        AtomicInteger ii=new AtomicInteger(0);
//        session.cs.getRegister().get(qw).forEach(
//                (node,locality) -> 
//                {
//                    if (ii.getAndIncrement()<50) {
//                        System.out.println(session.tree.getById(node));
//                        System.out.println("   "+locality);
//                    }
//                }
//        );
//
//        
//        System.out.println("");
        
        //some moficiation to the default tree to display the probas when a node is clisked
//        DefaultTreeCellRenderer renderer2 = (DefaultTreeCellRenderer) session.tree.getCellRenderer();
//        renderer.setLeafIcon(null);
//        renderer.setOpenIcon(null);
//        renderer.setClosedIcon(null);
//        session.tree.addTreeSelectionListener((TreeSelectionEvent e) -> {
//            
//            PhyloNode n=(PhyloNode)e.getNewLeadSelectionPath().getLastPathComponent();
//            if (!n.isLeaf()) {
//                double[] array = session.cs.getRegister().get(null).get(n.getId()).getAllProbas().toArray();
//                JFrame infos=new JFrame();
//                infos.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
//                infos.setLayout(new GridLayout(1, 1));
//                infos.add(ChartsForNodes.buildPPHistogramForANode(n, array));
//                infos.pack();
//                RefineryUtilities.centerFrameOnScreen(infos);
//                infos.setVisible(true);
//                
//            }
//        });
        
        //here put the action to open the graph corresponding to the node
        
        
        //session.tree.displayTree();
        //c.testMethod(); 

    }
    

}
