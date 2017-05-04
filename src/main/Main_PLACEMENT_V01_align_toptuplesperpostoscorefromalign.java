/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import alignement.Alignment;
import au.com.bytecode.opencsv.CSVWriter;
import charts.ChartsForNodes;
import core.AAStates;
import core.DNAStates;
import core.DiagSum;
import core.PProbasSorted;
import core.QueryWord;
import core.hash.SimpleHash;
import core.States;
import core.algos.SequenceKnife;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.awt.GridLayout;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.ui.RefineryUtilities;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import tree.NewickReader;
import tree.NewickWriter;
import tree.Tree;

/**
 * Algo in this version:
 * 
 * 1. all words of the query are searched in the hash, ONLY in a preselection of
 *    a few nodes (param nodeShift)
 *    --> int hese nodes best PP* hit corresponds to a node/position
 *    --> position is used for the alignment
 * 2. using the alignment, for all nodes we search query words associated to the
 *    positions memorized in the alignment
 *    --> we build diagsums for all nodes (1 per node)
 *    --> we keep as peeks only diagsum position holding more than X words
 *    --> the best score in these peeks define the placements
 * 
 * 
 * 
 * @author ben
 */
public class Main_PLACEMENT_V01_align_toptuplesperpostoscorefromalign {
    
    //sequence type
    public static int TYPE_DNA=1;
    public static int TYPE_PROTEIN=2;
    //memory mode
    public static int MEMORY_LOW=1;
    public static int MEMORY_LARGE=2;
    
    public static void main(String[] args) {

        try {
            System.setProperty("debug.verbose","1");
            
            ////////////////////////////////////////////////////////////////////
            // INPUT QUERY
            //pplacer benchmark, stat not based on relaxed tree
            String inputsPath="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/";
            //String q=inputsPath+"mod_p4z1r36_query_only2.fasta";
            //String q=inputsPath+"mod_p4z1r36_query_1st_seq_expanded.fasta";
            //String q=inputsPath+"mod_p4z1r36_query_ancestrals.fasta";
            String q=inputsPath+"mod_p4z1r36_query_ancestrals_rearranged.fasta";
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //PARAMETERS
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            
            //preplacement parameters//////////////////////////////////////////
            //for a number of nodes =total_nodes/nodeShift , build diagsum vectors
            //think to keep sum of the diagsums to do mean at the end and highlight positions > to mean
            int nodeShift=1; // carefull, brings an error when nodeShift<2
            if (nodeShift<2) {nodeShift=2;}
            // pour stocker ou non dans hash
            //minimum read/ref overlap,in bp. When not respected, read not reported
            int minOverlap=100;
            //word sampling method
            int queryWordSampling=SequenceKnife.SAMPLING_LINEAR;

            
            //debug/////////////////////////////////////////////////////////////
            //max number of queries treated 
            int queryLimit=3;
            //which log to write, !!!
            //more logs= much slower placement because of disk access latency
            boolean logDiagsums=false;
            boolean logPrePlacementDetailedDiagsums=false;
            //graph of words alignment
            boolean graphAlignment=true;
            

            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //LOADING ALL DATA FROM THE INPUT DATABASE
            long startLoadTime=System.currentTimeMillis();
            
            //LOAD SESSION//////////////////////////////////////////////////////
            //base path for outputs
            String path="/media/ben/STOCK/DATA/viromeplacer/WD/";
            //logs
            String logPath=path+"logs/";
            //trees
            String relaxedTreePath=path+"relaxed_trees/";
            //ancestral reconstruciton
            String ARPath=path+"AR/";
            //session itself
            boolean loadHash=true;
            SessionNext session= SessionNext.load(new File(path+"PAML_session_params_k8_mk8_f1.5_t3.9106607E-4"),loadHash);
            
            //type of Analysis//////////////////////////////////////////////////
            States s=session.states; 
            int analysisType=-1;
            //States: DNA or AA
            if (s instanceof DNAStates)
                analysisType=Main_PLACEMENT_V01_align_toptuplesperpostoscorefromalign.TYPE_DNA;
            else if (s instanceof AAStates)
                analysisType=Main_PLACEMENT_V01_align_toptuplesperpostoscorefromalign.TYPE_PROTEIN;
            
            //posterior probas analysis/////////////////////////////////////////
            //mers size
            int k=session.k;
            int min_k=session.minK;
            Infos.println("k="+k);
            Infos.println("min_k="+min_k);
            //site and word posterior probas thresholds
            float sitePPThreshold=session.stateThreshold;
            float wordPPStarThreshold=session.wordThreshold;
            float factor=session.factor;
            float thresholdAsLog=(float)Math.log10(wordPPStarThreshold);

            Infos.println("factor="+factor);
            Infos.println("sitePPThreshold="+sitePPThreshold+" (for info, unused)");
            Infos.println("wordPPStarThreshold="+wordPPStarThreshold);
            Infos.println("wordPPStarThreshold(log10)="+thresholdAsLog);
            



            //PREPARE DIRECTORIES///////////////////////////////////////////////
            if (!new File(path).exists()) {new File(path).mkdir();}
            if (!new File(logPath).exists()) {new File(logPath).mkdir();}
            if (!new File(relaxedTreePath).exists()) {new File(relaxedTreePath).mkdir();}
            if (!new File(ARPath).exists()) {new File(ARPath).mkdir();}
            
            
            //LOAD ORIGINAL ALIGNMENT///////////////////////////////////////////
            Alignment align=session.align;
            Infos.println(align.describeAlignment(false));
            
            //LOAD TREES////////////////////////////////////////////////////////
            NewickReader np=new NewickReader();
            Tree tree = session.tree;
            NewickWriter nw=new NewickWriter(new File("null"));
            Infos.println("# nodes in the tree: "+tree.getNodeCount());
            Infos.println("# leaves in the tree: "+tree.getLeavesCount());
            Infos.println("# internal nodes in the tree: "+tree.getInternalNodesByDFS().size());
            String relaxedTreeForJplace=nw.getNewickTree(tree, true, true, true);
            //Infos.println(relaxedTreeForJplace);
            nw.close();
            
            //LOAD THE POSTERIOR PROBAS/////////////////////////////////////////
            PProbasSorted pprobas = session.parsedProbas;
            //to raidly check that sorted probas are OK
            Infos.println("NodeId=0, 5 first PP:"+Arrays.deepToString(pprobas.getPPSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateIndexSet(0, 0, 5)));
            
            SimpleHash hash=session.hash;
            Infos.println(Environement.getMemoryUsage());
            long endLoadTime=System.currentTimeMillis();
            Infos.println("Loading the database took "+(endLoadTime-startLoadTime)+" ms");

            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //LOAD FINISHED, TIME TO START THE PLACEMENT PROCESS
            

            long startTotalTime=System.currentTimeMillis();

            
            
            ///////////////////////////////////////////////
            //SELECTION OF NODES SEARCHED DURING PREPLACEMENT
            int internalNodesCount=tree.getInternalNodesByDFS().size();
            int nodesToScan= internalNodesCount/nodeShift ;
            //nodeIdsTested[nodeId]=true -> the node will be used in preplamcement
            boolean[] nodeIdsTested=new boolean[tree.getNodeCount()]; 
            for (int i = 0; i < tree.getInternalNodesByDFS().size(); i++) {
                if (i%nodeShift==0) {
                    nodeIdsTested[tree.getInternalNodesByDFS().get(i)]=true;
                } else {
                    nodeIdsTested[tree.getInternalNodesByDFS().get(i)]=false;
                }
            }

            Infos.println("Number of internal nodes to scan: "+nodesToScan);
            

            ////////////////////////////////////////////////////////////////////
            //MAP TO REGISTER PLACEMENTS: map(query)=(map(nodeId)=positions_in_align)  !position in align, not diagsum
            HashMap<Fasta,HashMap<Integer,ArrayList<ArrayList<Object>>>> placementsPerQuery=new HashMap();
            
            LinkedHashMap<Fasta,ArrayList<Fasta>> identicalQueries=new LinkedHashMap();

            /////////////////////
            //PLACEMENT OF THE QUERIES
            FASTAPointer fp=new FASTAPointer(new File(q), false);
            FileWriter fw =new FileWriter(new File("Queries.fasta"));
            
            int queryCounter=0;
            int queryMatchingRefCounter=0;
            Fasta fasta=null;
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                if (queryCounter>=queryLimit)
                    break;
                
                Infos.println("#######################################################################");
                Infos.println("### PLACEMENT FOR QUERY #"+queryCounter+" : "+fasta.getHeader());
                Infos.println("#######################################################################");
                fw.append(fasta.getFormatedFasta()+"\n");
                long startScanTime=System.currentTimeMillis();
                int queryLength=fasta.getSequence().length();
                Infos.println("Query length: "+queryLength);
                
                //check if another query didn't had exactly the same sequence
                //fasta is a Comparable based on the sequence
                if (identicalQueries.containsKey(fasta)) {
                    identicalQueries.get(fasta).add(fasta);
                    Infos.println("Sequence identicial to "+identicalQueries.get(fasta).get(0).getHeader()+". Placement skipped.");
                    queryCounter++;
                    continue;
                }
                identicalQueries.put(fasta, new ArrayList<>());
                identicalQueries.get(fasta).add(fasta);
                
                //register query in placement table
                if (!placementsPerQuery.containsKey(fasta))
                    placementsPerQuery.put(fasta, new HashMap<>());
                



                
                ///////////////////////////////////
                // LOOP ON QUERY WORDS
                QueryWord qw=null;
                SequenceKnife sk=new SequenceKnife(fasta, k, min_k, s, queryWordSampling);
                //Infos.println("Mer order: "+Arrays.toString(rk.getMerOrder()));
                int queryWordCounter=0;
                int queryWordFoundCounter=0;
                
               
                //a diagsum based only on a the single top best PP* associated to a query word.
                //accumulating all nodes, just taing higher PP* and corresponding positions
                DiagSum bestPPStarsDiagsum=new DiagSum(queryLength, align.getLength(), minOverlap, k, sk.getStep());
                bestPPStarsDiagsum.init(thresholdAsLog);
                Infos.println("Diagsum vector size: "+bestPPStarsDiagsum.getSize());
                //a diagsum based only on a the top best PP* for for each alignment poistion associated to a query word.
                //accumulating all nodes, just taing higher PP* and corresponding positions
                DiagSum bestPPStarsPerPositionDiagsum=new DiagSum(queryLength, align.getLength(), minOverlap, k, sk.getStep());
                bestPPStarsPerPositionDiagsum.init(thresholdAsLog);
                

                //preparing preplacement graph data :
                double[][] graphDataForTopTuplesPerPos =null;
                double[][] graphDataForTopTuples =null;
                int xSize=bestPPStarsDiagsum.getSize();
                int ySize=queryLength;
                int xSize2=bestPPStarsPerPositionDiagsum.getSize();
                int ySize2=queryLength;
                
                if (graphAlignment) {
                    graphDataForTopTuples =new double[3][xSize*ySize];
                    //init the Z axis (PP*) to very small values for all possible (X,Y)
                        for (int y = 0; y < ySize; y++) {
                            for (int x = 0; x < xSize; x++) {
                                //System.out.println(col+" "+line+" "+(col+(line*session.align.getLength())));
                                graphDataForTopTuples[0][x+(y*xSize)]=x;
                                graphDataForTopTuples[1][x+(y*xSize)]=y;
                                graphDataForTopTuples[2][x+(y*xSize)]=thresholdAsLog;
                            }
                        }
                    graphDataForTopTuplesPerPos =new double[3][xSize2*ySize2];
                    //init the Z axis (PP*) to very small values for all possible (X,Y)
                        for (int y = 0; y < ySize2; y++) {
                            for (int x = 0; x < xSize2; x++) {
                                //System.out.println(col+" "+line+" "+(col+(line*session.align.getLength())));
                                graphDataForTopTuplesPerPos[0][x+(y*xSize2)]=x;
                                graphDataForTopTuplesPerPos[1][x+(y*xSize2)]=y;
                                graphDataForTopTuplesPerPos[2][x+(y*xSize2)]=thresholdAsLog;
                            }
                        }
                }
                    
                
                
                
                /////////////////////////////
                // LOG OUTPUT 1:words detais
                CSVWriter writerDetails=null;
                if (logPrePlacementDetailedDiagsums) {
                    Infos.println("Write logs: preplacement_diagSums_details.tsv");
                    writerDetails=new CSVWriter(new FileWriter(new File(logPath+fasta.getHeader()+"_preplacement_diagSums_details.tsv")), '\t');
                    String[] headerData=new String[2+bestPPStarsDiagsum.getSize()];
                    headerData[0]="best_node";
                    headerData[1]="word";
                    for (int i = 0; i < bestPPStarsDiagsum.getSize(); i++) {
                        headerData[i+2]="pos="+i;
                    }
                    writerDetails.writeNext(headerData);
                }
                /////////////////
                
                
                
                
                
                
                // MATRIX DESCRIBING ALIGNMENT
                //////////////////////////////////////////
                
                //[x][y] = [queryPosition][n_ieme_hit_in_alignment]
                ArrayList<Integer>[] alignmentMatrix=new ArrayList[queryLength-k+1];
                for (int i = 0; i < alignmentMatrix.length; i++) {
                    alignmentMatrix[i]=new ArrayList<>(5);
                }
                //Arrays.fill(alignmentMatrix,new ArrayList<Integer>()); //not working, ask why in internet
                
                
                // ITERTION ON ALL WORDS
                //////////////////////////////////////////
                
                while ((qw=sk.getNextWord())!=null) {
                    //Infos.println("Query mer: "+qw.toString());


                    SimpleHash.Tuple topTuple = hash.getTopTuple(qw);
                    //if this word is not registered in the hash
                    if (topTuple==null) {
                        queryWordCounter++;
                        continue;
                    }
                    //word is in the hash
                    queryWordFoundCounter++;
                    
                    /////////////////
                    // LOG OUTPUT 1
                    String[] data=null;
                    if (logPrePlacementDetailedDiagsums) {
                        data=new String[2+bestPPStarsDiagsum.getSize()];
                        data[0]=String.valueOf(topTuple.getNodeId());
                        data[1]=String.valueOf(qw.getOriginalPosition());
                    }
                    /////////////////    
                    
                    //get best PP* of each query word (whatever node/position)
                    //////////////////////////////////////////////////////////
                    
                    //a diagsum based only on a the single top best PP* associate dto a word.
                    int topDiagSumPos=topTuple.getRefPos()-qw.getOriginalPosition()+(queryLength-minOverlap);
                    if (topDiagSumPos>-1 && topDiagSumPos<bestPPStarsDiagsum.getSize()) {
                        bestPPStarsDiagsum.sum(topDiagSumPos, -thresholdAsLog+topTuple.getPPStar());
                        if (logPrePlacementDetailedDiagsums) {
                            data[2+topDiagSumPos]=String.valueOf(topTuple.getPPStar()); 
                        }
                        if(graphAlignment) {
                            graphDataForTopTuples[2][queryWordCounter*xSize+topTuple.getRefPos()]= topTuple.getPPStar();                        
                        }
                    }
                    
                    //get best PP* of each position --> nodes under nodeShift !!!
                    /////////////////////////////////////////////////
                    List<SimpleHash.Tuple> allTuples = hash.getTopTuplesUnderNodeShift(qw, thresholdAsLog, nodeIdsTested);
                    List<SimpleHash.Tuple> bestTuplePerPosition=new ArrayList<SimpleHash.Tuple>();
                    HashMap<Integer,Boolean> map=new HashMap<>();
                    for (int i = 0; i < allTuples.size(); i++) {
                        SimpleHash.Tuple tuple = allTuples.get(i);
                        if (!map.containsKey(tuple.getRefPos())) {
                            map.put(tuple.getRefPos(), true);
                            bestTuplePerPosition.add(tuple);
                        }
                    }
                    
                    for (int i = 0; i < bestTuplePerPosition.size(); i++) {
                        SimpleHash.Tuple tuple = bestTuplePerPosition.get(i);
                        int diagSumPos=tuple.getRefPos()-qw.getOriginalPosition()+(queryLength-minOverlap);
                        if (diagSumPos>-1 && diagSumPos<bestPPStarsPerPositionDiagsum.getSize()) {
                            //System.out.println("diagsumPos="+diagSumPos+" tuplePos="+tuple.getRefPos());
                            bestPPStarsPerPositionDiagsum.sum(diagSumPos, -thresholdAsLog+tuple.getPPStar());
                            if (graphAlignment) {
                                graphDataForTopTuplesPerPos[2][queryWordCounter*xSize2+tuple.getRefPos()]= tuple.getPPStar();   
                            }
                            //store in  alignment matrix
                            alignmentMatrix[qw.getOriginalPosition()].add(tuple.getRefPos());
                        }

                        
                    }
//                    if (bestTuplePerPosition.size()>1)
//                        System.out.println(bestTuplePerPosition.size()+" possible positions at query pos="+qw.getOriginalPosition());
                    
                    if(logPrePlacementDetailedDiagsums)
                        writerDetails.writeNext(data);  
                    
                    queryWordCounter++;
                    
                    //DEBUG
                    //if (queryWordCounter>1000)
                    //        break;
                    //DEBUG
                    
                }
                
                //check alignment vs diagsum stats
                for (int i = 0; i < alignmentMatrix.length; i++) {
                    System.out.println("alignment i="+i+": "+alignmentMatrix[i]);
                }
                //check diagsum stats
//                for (int i = 0; i < bestPPStarsPerPositionDiagsum.getSize(); i++) {
//                    System.out.println("pre-diagsum i="+i+": #words="+bestPPStarsPerPositionDiagsum.getEffectiveWordCount(i)+" PP*="+bestPPStarsPerPositionDiagsum.getSum(i));
//                }
                
                
                if(logPrePlacementDetailedDiagsums) {
                    writerDetails.flush();
                    writerDetails.close();
                }

                
                if (graphAlignment) {
                    //graph for topTuplePerPos
                    DefaultXYZDataset datasetForGraph=new DefaultXYZDataset();
                    //datasetForGraph.addSeries(0, graphDataForTopTuples);
                    datasetForGraph.addSeries(0, graphDataForTopTuplesPerPos);
                    JFrame infos=new JFrame();
                    infos.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                    infos.setLayout(new GridLayout(1, 1));
                    infos.add(ChartsForNodes.buildReadMatchForANode4("Top_tuple_per_pos of read="+fasta.getHeader(),bestPPStarsDiagsum.getSize(),queryLength, datasetForGraph,thresholdAsLog,0));
                    infos.setSize(1024, 250);
                    infos.pack();
                    RefineryUtilities.centerFrameOnScreen(infos);
                    infos.setVisible(true);
                }
               

                Infos.println("Proportion of query words retrieved in the hash: "+queryWordFoundCounter+"/"+queryWordCounter);
                long endScanTime=System.currentTimeMillis();
                Infos.println("Scan on "+nodesToScan+" nodes took "+(endScanTime-startScanTime)+" ms");
                //Infos.println("Diagsum of top positions:\n"+bestPPStarsPerPositionDiagsum);
           

              
                
                
               
                
                
                // HERE reuse the alignment vector to determine which positions 
                //should be tested for all nodes
                
                //load the PP* of these positions nodes per nodes
                //build the diagsums per nodes 
                //select peeks including #words > threshold
                //select the node with the best peeks as the best placement
                
                startScanTime=System.currentTimeMillis();
                //to register which node is assined to which line of preplacementAllDiagsum
                Infos.println("Launching search on all nodes...");
                int[] nodeIdToDiagsumListIndex = new int[tree.getNodeCount()];  //index[nodeId]=preplacementDiagSumIndex
                int[] diagsumListIndexToNodeId =new int[internalNodesCount];
                Arrays.fill(nodeIdToDiagsumListIndex,-1);
                Arrays.fill(diagsumListIndexToNodeId,-1);
                //the list of diagSums vectors one per internal node.
                ArrayList<DiagSum> diagsumsPerNode = new ArrayList<>(tree.getInternalNodesByDFS().size());
                for (int i=0;i<tree.getInternalNodesByDFS().size();i++) {
                    int nodeId=tree.getInternalNodesByDFS().get(i);
                    diagsumsPerNode.add(new DiagSum(queryLength, align.getLength(), minOverlap, k, sk.getStep()));
                    diagsumListIndexToNodeId[i]=nodeId;
                    nodeIdToDiagsumListIndex[nodeId]=i;
                    diagsumsPerNode.get(i).init(thresholdAsLog); 
                }
                
                
                
                //now, pass again on all algnment positions, but search them in 
                //all internal nodes.
                //only alignment positions whose coordinates matches a
                //diagsumPos associated to more than (20-k+1)/s words
                //(20 consecutive represented by consecutive words)
                //will be tested
                
                
                //V1
//                boolean[] positionSurvivingTest=new boolean[diagsumsPerNode.get(0).getSize()];
//                for (int i = 0; i < tree.getInternalNodesByDFS().size(); i++) {
//                    int nodeId=tree.getInternalNodesByDFS().get(i);
//                    for (int alignPos = 0; alignPos < alignmentMatrix.length; alignPos++) {
//                        QueryWord word=sk.getWordAt(alignPos);
//                        ArrayList<Integer> allRefPos=alignmentMatrix[alignPos];
//                        for (int refPos:allRefPos) {
//                            int diagSumPos=refPos-word.getOriginalPosition()+(queryLength-minOverlap);
//                            if (bestPPStarsPerPositionDiagsum.getEffectiveWordCount(diagSumPos)>=((20-k+1)/sk.getStep())) {
//                                SimpleHash.Tuple tuple = hash.getTuplePerNodeAndRefPosition(word, nodeId, refPos);
//                                positionSurvivingTest[diagSumPos]=true;
//                                if (tuple!=null) {
//                                    diagsumsPerNode.get(nodeIdToDiagsumListIndex[nodeId]).sum(diagSumPos, -thresholdAsLog+tuple.getPPStar());
//                                }
//                            }
//                        }
//                    }
//                }
                
                //V2, pickup nodes directly from tuples corrsponding to position, avoiding to check all nodes
                boolean[] positionSurvivingTest=new boolean[diagsumsPerNode.get(0).getSize()];
                boolean[] nodesAttained=new boolean[tree.getNodeCount()];
                for (int alignPos = 0; alignPos < alignmentMatrix.length; alignPos++) {
                    QueryWord word=sk.getWordAt(alignPos);
                    ArrayList<Integer> allRefPos=alignmentMatrix[alignPos];
                    for (int refPos:allRefPos) {
                        int diagSumPos=refPos-word.getOriginalPosition()+(queryLength-minOverlap);
                        if (bestPPStarsPerPositionDiagsum.getEffectiveWordCount(diagSumPos)>=((20-k+1)/sk.getStep())) {
                            positionSurvivingTest[diagSumPos]=true;
                            List<SimpleHash.Tuple> tuple = hash.getTuplePerRefPosition(word, refPos);
                            for (int i = 0; i < tuple.size(); i++) {
                                SimpleHash.Tuple get = tuple.get(i);
                                nodesAttained[get.getNodeId()]=true;
                                diagsumsPerNode.get(nodeIdToDiagsumListIndex[get.getNodeId()]).sum(diagSumPos, -thresholdAsLog+get.getPPStar());
                            }
                        }
                    }
                }            
                
                
                endScanTime=System.currentTimeMillis();
                Infos.println("Targeted scan on "+tree.getInternalNodesByDFS().size()+" nodes took "+(endScanTime-startScanTime)+" ms");
                //check how many peeks have susrvived
                int countSurviving=0;
                for (int i = 0; i < positionSurvivingTest.length; i++) {
                    if(positionSurvivingTest[i]==true)
                        countSurviving++;
                }
                Infos.println("Peeks surviving tests: "+countSurviving);
                //check how many peeks have susrvived
                int countNodesAttained=0;
                for (int i = 0; i < nodesAttained.length; i++) {
                    if(nodesAttained[i]==true)
                        countNodesAttained++;
                }
                Infos.println("Nodes attained: "+countNodesAttained);    
                
                /////////////////
                // LOG OUTPUT 2
                CSVWriter writerDiagsum=null;
                if (logDiagsums) {
                    Infos.println("Write logs: placement_diagSums.tsv");
                    writerDiagsum=new CSVWriter(new FileWriter(new File(logPath+fasta.getHeader()+"_placement_diagSums.tsv")), '\t');
                    String[] headerData=new String[3+diagsumsPerNode.size()];
                    headerData[0]="diagsum_position";
                    headerData[1]="reference_align_site";
                    headerData[2]="BPP";
                    for (int i = 0; i < diagsumsPerNode.size(); i++) {
                        headerData[i+3]="nodeid="+diagsumListIndexToNodeId[i];
                    }
                    writerDiagsum.writeNext(headerData);
                    

                    for (int i = 0; i < diagsumsPerNode.get(0).getSize(); i++) { //for each site
                        String[] data=new String[3+diagsumsPerNode.size()];
                        data[0]=String.valueOf(i);
                        data[1]=String.valueOf(i-(queryLength-minOverlap));
                        data[2]=String.valueOf(bestPPStarsPerPositionDiagsum.getSum(i));
                        for (int j = 0; j < diagsumsPerNode.size(); j++) {
                            data[j+3]=String.valueOf(diagsumsPerNode.get(j).getSum(i));
                        }
                        writerDiagsum.writeNext(data);
                    }
                    writerDiagsum.flush();
                    writerDiagsum.close();
                }
                /////////////////
                
                
                
                //as potential peeks are now detrmined, it's time to find the max values
                //we go through the diagsums and retain the max
                placementsPerQuery.put(fasta, new HashMap<>());
                for (int i = 0; i < positionSurvivingTest.length; i++) {
                    if(positionSurvivingTest[i]==true) {
                        float peekValue=Float.NEGATIVE_INFINITY;
                        int peekIndex=-1;
                        for (int j = 0; j < diagsumsPerNode.size(); j++) {
                            float val=diagsumsPerNode.get(j).getSum(i);
                            if (val>peekValue) {
                                peekValue=val;
                                peekIndex=j;
                            }
                        }
                        Infos.println("For peek at diagSumPos="+i+" nodeId="+diagsumListIndexToNodeId[peekIndex]+" appears to be the best placement.");
                        Infos.println("Node details: "+tree.getById(diagsumListIndexToNodeId[peekIndex]));
                        placementsPerQuery.get(fasta).put(diagsumListIndexToNodeId[peekIndex], new ArrayList<>());
                        ArrayList<Object> vals=new ArrayList<>();
                        vals.add(diagsumsPerNode.get(0).getPositionInReference(i)); //position
                        vals.add(peekValue); //PP*
                        placementsPerQuery.get(fasta).get(diagsumListIndexToNodeId[peekIndex]).add(vals);
                    }
                }                
                
                
                
                
                
                
                

                
                queryCounter++;
            }
            
            fw.close();
            fp.closePointer();
            
            

            
            Infos.println("#######################################################################");
            Infos.println("### "+queryCounter+" READS WERE ANALYZED, NOW GENERATING OUTPUTS...");
            Infos.println("#######################################################################");
            
            
            
            ////////////////////////////////////////////////////////////////////
            //OUTPUT THE PLACEMENTS IN TSV FORMAT
            ////////////////////////////////////////////////////////////////////
            Infos.println("Generating placement output in .tsv format...");
            CSVWriter fwPlacement=new CSVWriter(new FileWriter(new File(logPath+"placements.tsv")), '\t');
            String[] header=new String[5];
            header[0]="Query";
            header[1]="NodeId";
            header[2]="ExtNodeId";
            header[3]="Align_start";
            header[4]="PP*";
            fwPlacement.writeNext(header);
            String[] data=new String[5];
            for (Iterator<Fasta> iterator = placementsPerQuery.keySet().iterator(); iterator.hasNext();) {
                Fasta fastaMaster = iterator.next();
                for (Iterator<Fasta> iteratorSlave = identicalQueries.get(fastaMaster).iterator(); iteratorSlave.hasNext();) {
                    Fasta fastaSlave = iteratorSlave.next();
                    data[0]=fastaSlave.getHeader();
                    for (Iterator<Integer> iterator1 = placementsPerQuery.get(fastaSlave).keySet().iterator(); iterator1.hasNext();) {
                        Integer nextNode = iterator1.next();
                        data[1]=String.valueOf(nextNode);
                        data[2]=String.valueOf(tree.getById(nextNode).getExternalId());
                        for (Iterator<ArrayList<Object>> iterator2 = placementsPerQuery.get(fastaSlave).get(nextNode).iterator(); iterator2.hasNext();) {
                            ArrayList<Object> nextPosition = iterator2.next();
                            data[3]=String.valueOf((Integer)nextPosition.get(0));
                            data[4]=String.valueOf((Float)nextPosition.get(1));
                            fwPlacement.writeNext(data);
                        }
                    }
                }
            }
            fwPlacement.close();
            
            
            
            
            ////////////////////////////////////////////////////////////////////
            //OUTPUT THE PLACEMENTS IN JPLACE (JSON) FORMAT
            ////////////////////////////////////////////////////////////////////
            Infos.println("Generating placement output in .jplace format...");
            //with library json.simple
            int level=0;
            JSONObject top=new JSONObject(); 
            LinkedHashMap topMap=new LinkedHashMap();
            //object tree (mandatory)
            topMap.put("tree",relaxedTreeForJplace);
            //we do one placement object for each group of identical reads
            JSONArray placements=new JSONArray();
            for (Iterator<Fasta> iterator = identicalQueries.keySet().iterator(); iterator.hasNext();) {
                Fasta currentFasta = iterator.next();
                
                //object representing placements (mandatory), it contains 2 key/value couples
                JSONObject placement=new JSONObject();
                
                //first we build the "p" object, containing the position/scores of all reads
                JSONArray allPlaces=new JSONArray();
                //map(query)=(map(nodeId)=positions_in_align)
                for (Iterator<Integer> iterator1 = placementsPerQuery.get(currentFasta).keySet().iterator(); iterator1.hasNext();) {
                    //this is the node to which the read was placed
                    Integer nodeId = iterator1.next();
                    for (Iterator<ArrayList<Object>> iterator2 = placementsPerQuery.get(currentFasta).get(nodeId).iterator(); iterator2.hasNext();) {
                        //this is the position of the peek
                        ArrayList<Object> next2 = iterator2.next();
                        JSONArray placeColumns=new JSONArray();
                        placeColumns.add(nodeId);  //nodeId
                        placeColumns.add((Integer)next2.get(0)); //ref_position
                        placeColumns.add((Float)next2.get(1)); //PP*
                        allPlaces.add(placeColumns);
                    }
                }
                placement.put("p", allPlaces);
                
                //second we build the nm object, containing the reads identifiers
                JSONArray allIdentifiers=new JSONArray();
                //these are simply all the identical sequences
                for (Iterator<Fasta> iterator1 = identicalQueries.get(currentFasta).iterator(); iterator1.hasNext();) {
                    Fasta next = iterator1.next();
                    JSONArray readMultiplicity=new JSONArray();
                    readMultiplicity.add(next.getHeader());
                    readMultiplicity.add(1);
                    allIdentifiers.add(readMultiplicity);
                }
                placement.put("nm", allIdentifiers);
                
                //store the placement in the list of placements
                placements.add(placement);
            }
            //associate the list of placements to the top level
            top.put("placements",placements);
            //object version (mandatory
            top.put("version",3);
            //object metadata (mandatory): sub-objet "version is mandatory"
            JSONObject invoc=new JSONObject();
            invoc.put("invocation", "metaplacer xxxxxxxxxxxxxxx");
            top.put("metadata", invoc);
            //object fields
            JSONArray fList=new JSONArray();
            fList.add("edge_id");
            fList.add("son_node_id");
            fList.add("son_node_externa_id");
            fList.add("reference_localization");
            top.put("fields", fList);
            
            //put all the elements in the top JSON object
            top.putAll(topMap);
            String out=top.toJSONString();

            out=out.replaceAll("\\},\\{", "\n\\},\\{\n\t"); //},{
            out=out.replaceAll("\\],\"","\\],\n\t\"");   //],"
            out=out.replaceAll("\\]\\}\\],", "\\]\n\\}\n\\},\n"); //]}]
            //out=out.replace("]},", "]},"); //]}
            
            FileWriter fwJSON =new FileWriter(new File(logPath+"placements.jplace"));
            fwJSON.append(out);
            fwJSON.close();
                    
            
            

            long endTotalTime=System.currentTimeMillis();
            Infos.println("#######################################################################");
            Infos.println("### DONE, total execution (excluding DB load): "+(endTotalTime-startTotalTime)+" ms");
            Infos.println("#######################################################################");
            
            
            
            
            //System.exit(0);
            
            
        } catch (IOException ex) {
            Logger.getLogger(Main_PLACEMENT_V01_align_toptuplesperpostoscorefromalign.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        
        
        
    }
    

    
    
    public static void inputStreamToOutputStream(final InputStream inputStream, final OutputStream out) {
        Thread t = new Thread(new Runnable() {

            @Override
            public void run() {
                try {
                    int d;
                    while ((d = inputStream.read()) != -1) {
                        out.write(d);
                    }
                } catch (IOException ex) {
                    Infos.println(ex.getCause());
                }
            }
        });
        t.setDaemon(true);
        t.start();
    }
    
    
    
    
    
}
