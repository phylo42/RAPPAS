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
import core.older.PProbas;
import core.PProbasSorted;
import core.QueryWord;
import core.SimpleHash;
import core.States;
import core.algos.SequenceKnife;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.PAMLWrapper;
import java.awt.GridLayout;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.ProcessBuilder.Redirect;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import static main.Main_DBBUILD.TYPE_DNA;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.ui.RefineryUtilities;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import prog.ProgRunner;
import tree.NewickReader;
import tree.NewickWriter;
import tree.PhyloNode;
import tree.PhyloTree;
import tree.ExtendedTree;
import tree.Tree;

/**
 *
 * @author ben
 */
public class Main_PLACEMENT_3 {
    
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
            String q=inputsPath+"mod_p4z1r36_query_ancestrals.fasta";
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //PARAMETERS
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            
            //preplacement parameters//////////////////////////////////////////
            //for a number of nodes =total_nodes/nodeShift , build diagsum vectors
            //think to keep sum of the diagsums to do mean at the end and highlight positions > to mean
            int nodeShift=100; // carefull, brings an error when nodeShift<2
            // pour stocker ou non dans hash
            //minimum read/ref overlap,in bp. When not respected, read not reported
            int minOverlap=100;
            //word sampling method
            int preplacementWordSampling=SequenceKnife.SAMPLING_LINEAR;
            //peek detection parameters
            int w=31; //sliding window size
            int w_min=17;//minimum window size (for diagSums borders, basically, the w_win/2 1st position will not show a ratio, i.e. output value is -1.0)
            if (w_min>w) {w_min=w;}//!!! must be w_min<=w. if not, the program force w_min=w
            double r_win=0.8; //peek/windows_mean ratio
            double r_nodes=0.20; //node_where_peek_detected/scanned_nodes ratio
            int maxAllowedPeeks=50; //maximum number of localizations allowed for a read, if greater read discarded because clearly a repeat
            //fine placement parameters//////////////////////////////////////////
            //word sampling method
            int wordSampling=SequenceKnife.SAMPLING_NON_OVERLAPPING;
            //searchMethod
            boolean deepSearch=false;
            //when doing deepsearch, bp on left and right around the (peek;peek+query_length) interval
            int peekBoundary=10; 
            //the 'top' nodes (of highest score) are reported for a peek
            //int topNodesSize=5;
            
            //debug/////////////////////////////////////////////////////////////
            //max number of queries treated 
            int queryLimit=5;
            //which log to write, !!!
            //more logs= much slower placement because of disk access latency
            boolean logDetailedDiagsums=false;
            boolean logDiagsums=false;
            boolean logPeekRatios=false;
            boolean logPreplacementDiagsums=true;
            boolean logPreplacementDetailedDiagsums=true;

            

            
            
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
                analysisType=Main_PLACEMENT_3.TYPE_DNA;
            else if (s instanceof AAStates)
                analysisType=Main_PLACEMENT_3.TYPE_PROTEIN;
            
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
            

          

            
            

            ////////////////////////////////////////////////////////////////////
            //MAP TO REGISTER PLACEMENTS: map(query)=(map(nodeId)=positions_in_align)  !position in align, not diagsum
            HashMap<Fasta,HashMap<Integer,ArrayList<Integer>>> placementsPerQuery=new HashMap();
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
                
                //infos about nodes to scan
                int internalNodesCount=tree.getInternalNodesByDFS().size();
                int nodesToScan= internalNodesCount/nodeShift ;
                //alternatively, force node for tests
                
                Infos.println("Number of internal nodes to scan: "+nodesToScan);
                int queryLength=fasta.getSequence().length();
                Infos.println("Query length: "+queryLength);
                //to register which node is assined to which line of preplacementAllDiagsum
                int[] preplaceDiagIndex = new int[tree.getNodeCount()];  //index[nodeId]=preplacementDiagSumIndex
                int[] preplaceDiagNodeIndex =new int[internalNodesCount];
                Arrays.fill(preplaceDiagIndex,-1);
                Arrays.fill(preplaceDiagNodeIndex,-1);
                
                ///////////////////////////////////
                // LOOP ON QUERY WORDS
                QueryWord qw=null;
                SequenceKnife sk=new SequenceKnife(fasta, k, min_k, s, preplacementWordSampling);
                //Infos.println("Mer order: "+Arrays.toString(rk.getMerOrder()));
                int queryWordCounter=0;
                int queryWordFoundCounter=0;
                
                
                //the list of diagSums vectors one per internal node.
                ArrayList<DiagSum> preplacementAllDiagsums = new ArrayList<>(nodesToScan);
                boolean[] nodeScannedInPreplacement=new boolean[tree.getNodeCount()];
                Arrays.fill(nodeScannedInPreplacement,false);
                for (int i=0;i<nodesToScan;i++) {
                    int nodeId=tree.getInternalNodesByDFS().get(i*nodeShift);
                    preplacementAllDiagsums.add(new DiagSum(queryLength, align.getLength(), minOverlap, k, sk.getStep()));
                    nodeScannedInPreplacement[nodeId]=true;
                    //System.out.println("i:"+nodeId);
                    preplaceDiagNodeIndex[i]=nodeId;
                    preplaceDiagIndex[nodeId]=i;
                    preplacementAllDiagsums.get(i).init(thresholdAsLog); 
                }
                Infos.println("Diagsum vector size: "+preplacementAllDiagsums.get(0).getSize());
                
                //a diagsum based only on a the single top best PP* associated to a query word.
                //accumulating all nodes, just taing higher PP* and corresponding positions
                DiagSum bestPPStarsDiagsum=new DiagSum(queryLength, align.getLength(), minOverlap, k, sk.getStep());
                bestPPStarsDiagsum.init(thresholdAsLog);
                
                //a diagsum based only on a the top best PP* for for each alignment poistion associated to a query word.
                //accumulating all nodes, just taing higher PP* and corresponding positions
                DiagSum bestPPStarsPerPositionDiagsum=new DiagSum(queryLength, align.getLength(), minOverlap, k, sk.getStep());
                bestPPStarsPerPositionDiagsum.init(thresholdAsLog);
                
                
                
                //preparing graph data:
                int xSize=bestPPStarsDiagsum.getSize();
                int ySize=queryLength;
                double[][] graphDataForTopTuples =new double[3][xSize*ySize];
                //init the Z axis (PP*) to very small values for all possible (X,Y)
                    for (int y = 0; y < ySize; y++) {
                        for (int x = 0; x < xSize; x++) {
                            //System.out.println(col+" "+line+" "+(col+(line*session.align.getLength())));
                            graphDataForTopTuples[0][x+(y*xSize)]=x;
                            graphDataForTopTuples[1][x+(y*xSize)]=y;
                            graphDataForTopTuples[2][x+(y*xSize)]=thresholdAsLog;
                        }
                    }
                int xSize2=bestPPStarsPerPositionDiagsum.getSize();
                int ySize2=queryLength;
                double[][] graphDataForTopTuplesPerPos =new double[3][xSize2*ySize2];
                //init the Z axis (PP*) to very small values for all possible (X,Y)
                    for (int y = 0; y < ySize2; y++) {
                        for (int x = 0; x < xSize2; x++) {
                            //System.out.println(col+" "+line+" "+(col+(line*session.align.getLength())));
                            graphDataForTopTuplesPerPos[0][x+(y*xSize2)]=x;
                            graphDataForTopTuplesPerPos[1][x+(y*xSize2)]=y;
                            graphDataForTopTuplesPerPos[2][x+(y*xSize2)]=thresholdAsLog;
                        }
                    }
                    
                
                //DEBUG
                boolean testWord=false;
                FileWriter hashBucketLoadWriter=null;
                if (testWord) {
                    hashBucketLoadWriter=new FileWriter(new File(path+"hash_bucket_load_test_k"+k+"_f"+factor));
                    hashBucketLoadWriter.write("word,nodeId,position,PPStar\n");
                }
                //DEBUG
                
                
                /////////////////////////////
                // LOG OUTPUT 1:words detais
                CSVWriter writerDetails=null;
                if (logPreplacementDetailedDiagsums) {
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
                    if (logPreplacementDetailedDiagsums) {
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
                        if (logPreplacementDetailedDiagsums) {
                            data[2+topDiagSumPos]=String.valueOf(topTuple.getPPStar()); 
                        }
                        graphDataForTopTuples[2][queryWordCounter*xSize+topTuple.getRefPos()]= topTuple.getPPStar();                        
                    }
                    
                    //get best PP* of each position (whatever nodes)
                    /////////////////////////////////////////////////
                    List<SimpleHash.Tuple> allTuples = hash.getTopTuplesUnderNodeShift(qw, -1.0f, nodeScannedInPreplacement);
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
                            bestPPStarsPerPositionDiagsum.sum(diagSumPos, -thresholdAsLog+tuple.getPPStar());
                            graphDataForTopTuplesPerPos[2][queryWordCounter*xSize2+tuple.getRefPos()]= tuple.getPPStar();   
                            //store in  alignment matrix
                            alignmentMatrix[qw.getOriginalPosition()].add(tuple.getRefPos());
                        }

                        
                    }
                    if (bestTuplePerPosition.size()>1)
                        System.out.println(bestTuplePerPosition.size()+" possible positions at query pos="+qw.getOriginalPosition());
                    
                    
                    
                    
                    
//                    //System.out.println("Tuple size: "+allTuples.size());
//                    for (int i = 0; i < allTuples.size(); i++) {
//                        SimpleHash.Tuple tuple = allTuples.get(i);
//                        if (testWord)
//                            hashBucketLoadWriter.write(queryWordCounter+","+tuple.toStringCSV()+"\n");
//                        //System.out.println(tuple);
//                        int diagSumPos=tuple.getRefPos()-qw.getOriginalPosition()+(queryLength-minOverlap);
//                        //System.out.println("Match diagsumPos:"+diagSumPos);
//                        if (diagSumPos>-1 && diagSumPos<preplacementAllDiagsums.get(0).getSize()) {
//                            //System.out.println("Modified: "+(-thresholdAsLog+tuple.getPPStar()));
//                            preplacementAllDiagsums.get(preplaceDiagIndex[tuple.getNodeId()]).sum(diagSumPos, -thresholdAsLog+tuple.getPPStar());
//                            float val=preplacementAllDiagsums.get(preplaceDiagIndex[tuple.getNodeId()]).getSum(diagSumPos);
//
//                        }
//                        
//                    }


                    










                    
                    if(logPreplacementDetailedDiagsums)
                        writerDetails.writeNext(data);  
                    
                    queryWordCounter++;
                    
                    //DEBUG
                    //if (queryWordCounter>1000)
                    //        break;
                    //DEBUG
                    
                }
                
                
                for (int i = 0; i < alignmentMatrix.length; i++) {
                    System.out.println("alignment i="+i+": "+alignmentMatrix[i]);
                }
                
                
                
                if(logPreplacementDetailedDiagsums) {
                    writerDetails.flush();
                    writerDetails.close();
                }

                if (hashBucketLoadWriter!=null)
                    hashBucketLoadWriter.close(); 
                
                
                //graph for topTuples/topTuplePerPos
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
                

                Infos.println("Proportion of query words found: "+queryWordFoundCounter+"/"+queryWordCounter);
                long endScanTime=System.currentTimeMillis();
                Infos.println("Scan on "+nodesToScan+" nodes took "+(endScanTime-startScanTime)+" ms");
                
           
                
                
                
                
                
                /////////////////
                // LOG OUTPUT 2
                CSVWriter writerDiagsum=null;
                if (logPreplacementDiagsums) {
                    Infos.println("Write logs: preplacement_diagSums.tsv");
                    writerDiagsum=new CSVWriter(new FileWriter(new File(logPath+fasta.getHeader()+"_preplacement_diagSums.tsv")), '\t');
                    String[] headerData=new String[3+preplacementAllDiagsums.size()];
                    headerData[0]="diagsum_position";
                    headerData[1]="reference_align_site";
                    headerData[2]="BEST_MIX";
                    for (int i = 0; i < preplacementAllDiagsums.size(); i++) {
                        headerData[i+3]="nodeid="+preplaceDiagNodeIndex[i];
                    }
                    writerDiagsum.writeNext(headerData);
                    

                    for (int i = 0; i < preplacementAllDiagsums.get(0).getSize(); i++) { //for each site
                        String[] data=new String[3+preplacementAllDiagsums.get(0).getSize()];
                        data[0]=String.valueOf(i);
                        data[1]=String.valueOf(i-(queryLength-minOverlap));
                        data[2]=String.valueOf(bestPPStarsDiagsum.getSum(i));
                        for (int j = 0; j < preplacementAllDiagsums.size(); j++) {
                            data[j+3]=String.valueOf(preplacementAllDiagsums.get(j).getSum(i));
                        }
                        writerDiagsum.writeNext(data);
                    }
                    writerDiagsum.flush();
                    writerDiagsum.close();
                }
                /////////////////
              
                
                
                queryCounter++;
            }
            
            fw.close();
            fp.closePointer();
            
            

            
            
            
            
            
            
            ////////////////////////////////////////////////////////////////////
            //OUTPUT THE PLACEMENTS IN TSV FORMAT
            ////////////////////////////////////////////////////////////////////
            
            CSVWriter fwPlacement=new CSVWriter(new FileWriter(new File(logPath+"placements.tsv")), '\t');
            String[] header=new String[4];
            header[0]="Query";
            header[1]="NodeId";
            header[2]="ExtNodeId";
            header[3]="Loc_start";
            fwPlacement.writeNext(header);
            String[] data=new String[4];
            for (Iterator<Fasta> iterator = placementsPerQuery.keySet().iterator(); iterator.hasNext();) {
                Fasta fastaMaster = iterator.next();
                for (Iterator<Fasta> iteratorSlave = identicalQueries.get(fastaMaster).iterator(); iteratorSlave.hasNext();) {
                    Fasta fastaSlave = iteratorSlave.next();
                    data[0]=fastaSlave.getHeader();
                    for (Iterator<Integer> iterator1 = placementsPerQuery.get(fastaSlave).keySet().iterator(); iterator1.hasNext();) {
                        Integer nextNode = iterator1.next();
                        data[1]=String.valueOf(nextNode);
                        data[2]=String.valueOf(tree.getById(nextNode).getExternalId());
                        for (Iterator<Integer> iterator2 = placementsPerQuery.get(fastaSlave).get(nextNode).iterator(); iterator2.hasNext();) {
                            Integer nextPosition = iterator2.next();
                            data[3]=String.valueOf(nextPosition);
                            fwPlacement.writeNext(data);
                        }
                    }
                }
            }
            fwPlacement.close();
            
            
            
            
            ////////////////////////////////////////////////////////////////////
            //OUTPUT THE PLACEMENTS IN JPLACE (JSON) FORMAT
            ////////////////////////////////////////////////////////////////////
            //with library json.simple
            int level=0;
            JSONObject top=new JSONObject(); 
            LinkedHashMap topMap=new LinkedHashMap();
            //object tree (mandatory)
            topMap.put("tree",relaxedTreeForJplace);
            //object placements (mandatory)
            topMap.put("placements", new ArrayList());
            for (Iterator<Fasta> iterator = identicalQueries.keySet().iterator(); iterator.hasNext();) {
                Fasta currentFasta = iterator.next();
                JSONObject p=new JSONObject();
                //will add section "p" and section "nm" in the object p
                
                
                for (Iterator<Integer> iterator1 = placementsPerQuery.get(currentFasta).keySet().iterator(); iterator1.hasNext();) {
                    Integer next1 = iterator1.next();

                    for (Iterator<Integer> iterator2 = placementsPerQuery.get(currentFasta).get(next1).iterator(); iterator2.hasNext();) {
                        Integer next2 = iterator2.next();

                        
                    }
                }
                //topMap.put("p", );
                //topMap.put("nm", );
                
            }
            top.put("placements","XXX");
            
            //object version (mandatory
            top.put("version",3);
            //object metadata (mandatory): sub-objet "version is mandatory"
            JSONObject invoc=new JSONObject();
            invoc.put("invocation", "viromeplacer xxxxxxxxxxxxxxx");
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
            System.out.println(out);
            
            
            

            System.out.println("Result: "+queryMatchingRefCounter+"/"+queryCounter+" needed to be localized.");
            
            
            
            
            
            
            //System.exit(0);
            
            
        } catch (IOException ex) {
            Logger.getLogger(Main_PLACEMENT_3.class.getName()).log(Level.SEVERE, null, ex);
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
